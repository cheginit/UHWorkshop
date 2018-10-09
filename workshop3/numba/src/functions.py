import numpy as np
from numba import njit, stencil, prange


class Grid2D(object):
    """ Generates a mesh assuming equal spacing in x and y direction.
        Also it allocates required arrays for velocity and pressure.

        Args:
            ngx (int): Number of grid points in x-direction
            ngy (int): Number of grid points in y-direction
            l_lid (float): Dimension of the box

        Returns:
            Grid class
    """
    def __init__(self, ngx, ngy, l_lid):
        self.l_lid = l_lid
        self.ngx = ngx
        self.ngy = ngy

        self.xlo = 1
        self.ylo = 1
        self.xhi = ngx - 1
        self.yhi = ngy - 1
        self.xtot = ngx + 1
        self.ytot = ngy + 1

        self.dx = l_lid / np.float(ngx - 1)
        self.dy = self.dx

        # Variables on staggered grid points
        self.u = np.zeros([self.ngx, self.ytot], dtype=np.float64)
        self.v = np.zeros([self.xtot, self.ngy], dtype=np.float64)
        self.p = np.zeros([self.xtot, self.ytot], dtype=np.float64)

        # Variabls on original grid points
        self.ug = np.zeros_like(self.u)
        self.vg = np.zeros_like(self.v)
        self.pg = np.zeros_like(self.p)

        # Velocities on cell centers
        self.uc = np.zeros_like(self.u)
        self.vc = np.zeros_like(self.v)

        # Set of boundary conditions
        self.BC = {'u': {'t': 1.0, 'b': 0.0, 'r': 0.0, 'l': 0.0},
                   'v': {'t': 0.0, 'b': 0.0, 'r': 0.0, 'l': 0.0},
                   'p': {'t': 0.0, 'b': 0.0, 'r': 0.0, 'l': 0.0}}

    # Dirichlet BC for velocity field
    def BC_u(self):
        self.u[self.xlo:self.xhi, 0] = \
            self.BC['u']['b'] - self.u[self.xlo:self.xhi, 1]
        self.u[self.xlo:self.xhi, self.ytot - 1] = \
            2.0 * self.BC['u']['t'] - self.u[self.xlo:self.xhi, self.ytot - 2]
        self.u[0, :self.ytot] = self.BC['u']['l']
        self.u[self.xhi, :self.ytot] = self.BC['u']['r']

    def BC_v(self):
        self.v[self.xlo:self.xtot - 1, 0] = self.BC['v']['b']
        self.v[self.xlo:self.xtot - 1, self.yhi] = self.BC['v']['t']
        self.v[0, :self.ngy] = \
            self.BC['v']['l'] - self.v[1, :self.ngy]
        self.v[self.xtot - 1, :self.ngy] = \
            self.BC['v']['r'] - self.v[self.xtot - 2, :self.ngy]

    # Neumann BC for pressure field
    def BC_p(self):
        self.p[self.xlo:self.xtot - 1, 0] = \
            self.p[self.xlo:self.xtot - 1, 1] - self.dy * self.BC['p']['b']
        self.p[self.xlo:self.xtot - 1, self.ytot - 1] = \
            self.p[self.xlo:self.xtot - 1, self.ytot - 2] - \
            self.dy * self.BC['p']['t']
        self.p[0, :self.ytot] = \
            self.p[1, :self.ytot] - self.dx * self.BC['p']['l']
        self.p[self.xtot - 1, :self.ytot] = \
            self.p[self.xtot - 2, :self.ytot] - \
            self.dx * self.BC['p']['r']


class Simulation(object):
    """ Generates simulation parameters and solve momentum and continuity
        equations.

        Args:
            grid (class): 2D grid class
            cfl (float): Courant number
            c2 (float): Squared artificial sound speed (m^2/s^2)
            Re (float): Reynolds number
            tol (float): Desired tolerence (default = 1e-7)
            itr_max (int): Maximum number of iterations

        Returns:
            Simulation Class
    """
    def __init__(self, grid, cfl, c2, Re, tol=1e-7, itr_max=1e6):
        self.grid = grid
        self.dt = cfl * min(grid.dx, grid.dy) / grid.BC['u']['t']
        self.nu = grid.BC['u']['t'] * grid.l_lid / Re
        self.c2 = c2
        self.re = Re
        self.tol = tol
        self.itr_max = itr_max

        self.dtxy = self.dt * grid.dx * grid.dy

        self.init_cond()
        grid.BC_u()
        grid.BC_v()
        grid.BC_p()

    # Impose initial condition
    def init_cond(self):
        g = self.grid
        g.u[g.xlo:g.xhi, g.ytot - 1] = g.BC['u']['t']
        g.u[g.xlo:g.xhi, g.ytot - 2] = g.BC['u']['t']

    # Compute L2-norm based on new and pold time steps data
    def l2norm(self, phi_n, phi_o):
        return np.sqrt(self.dtxy * np.sum((phi_n - phi_o)**2))


# Compute average in x-direction
@stencil
def ave_x(arr):
    return 0.5 * (arr[-1, 0] + arr[0, 0])


# Compute average in y-direction
@stencil
def ave_y(arr):
    return 0.5 * (arr[0, -1] + arr[0, 0])


# Compute average in 2D
@stencil
def ave_xy(arr):
    return 0.25 * (arr[0, 0] + arr[0, -1] + arr[-1, 0] + arr[-1, -1])


# Solve momentum equation in x and y direction
@njit(parallel=True)
def momentum(uo, vo, po, dt, dx, dy, nu):
    un = uo.copy()
    vn = vo.copy()

    dtdx = dt / dx
    dtdy = dt / dy
    dtdxx = dt / (dx * dx)
    dtdyy = dt / (dy * dy)

    ug = ave_y(uo)
    uc = ave_x(uo)
    vg = ave_x(vo)
    vc = ave_y(vo)

    for i in prange(1, uo.shape[0] - 1):
        for j in prange(1, uo.shape[1] - 1):
            un[i, j] = uo[i, j] \
                - dtdx * (uc[i+1, j] * uc[i+1, j] - uc[i, j] * uc[i, j]) \
                - dtdy * (ug[i, j+1] * vg[i+1, j] - ug[i, j] * vg[i+1, j-1]) \
                - dtdx * (po[i+1, j] - po[i, j]) \
                + nu * (dtdxx * (uo[i-1, j] - 2 * uo[i, j] + uo[i+1, j])
                        + dtdyy * (uo[i, j-1] - 2 * uo[i, j] + uo[i, j+1]))

    for i in prange(1, vo.shape[0] - 1):
        for j in prange(1, vo.shape[1] - 1):
            vn[i, j] = vo[i, j] \
                - dtdx * (ug[i, j+1] * vg[i+1, j] - ug[i-1, j+1] * vg[i, j]) \
                - dtdy * (vc[i, j+1] * vc[i, j+1] - vc[i, j] * vc[i, j]) \
                - dtdy * (po[i, j+1] - po[i, j]) \
                + nu * (dtdxx * (vo[i-1, j] - 2 * vo[i, j] + vo[i+1, j])
                        + dtdyy * (vo[i, j-1] - 2 * vo[i, j] + vo[i, j+1]))
    return un, vn


# Solve continuity equation
@njit(parallel=True)
def continuity(uo, vo, po, dt, dx, dy, c2):
    dtdx = dt / dx
    dtdy = dt / dy

    cn_err = np.zeros_like(po)
    pn = po.copy()

    err = 0.0
    for i in prange(1, vo.shape[0] - 1):
        for j in prange(1, uo.shape[1] - 1):
            cn_err[i, j] = dtdx * (uo[i, j] - uo[i-1, j]) \
                + dtdy * (vo[i, j] - vo[i, j-1])
            pn[i, j] = po[i, j] - c2 * cn_err[i, j]
            err += cn_err[i, j]

    return pn, err
