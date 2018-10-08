import numpy as np


class Grid2D(object):
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

        self.u = np.zeros([self.ngx, self.ytot], dtype=np.float64)
        self.v = np.zeros([self.xtot, self.ngy], dtype=np.float64)
        self.p = np.zeros([self.xtot, self.ytot], dtype=np.float64)

        self.ug = np.zeros_like(self.u)
        self.vg = np.zeros_like(self.v)
        self.pg = np.zeros_like(self.p)

        self.uc = np.zeros_like(self.u)
        self.vc = np.zeros_like(self.v)

        # Set of boundary conditions
        self.BC = {'u': {'t': 1.0, 'b': 0.0, 'r': 0.0, 'l': 0.0},
                   'v': {'t': 0.0, 'b': 0.0, 'r': 0.0, 'l': 0.0},
                   'p': {'t': 0.0, 'b': 0.0, 'r': 0.0, 'l': 0.0}}

    # Compute velocity field on grid points and cell centers
    def phi_g(self):
        self.ug[:, self.ylo:] = 0.5 * (self.u[:, self.ylo:]
                                       + self.u[:, self.ylo - 1:-1])
        self.vg[self.xlo:, :] = 0.5 * (self.v[self.xlo:, :]
                                       + self.v[self.xlo - 1:-1, :])

    def phi_c(self):
        self.uc[self.xlo:, :] = 0.5 * (self.u[self.xlo:, :]
                                       + self.u[self.xlo - 1:-1, :])
        self.vc[:, self.ylo:] = 0.5 * (self.v[:, self.ylo:]
                                       + self.v[:, self.ylo - 1:-1])

    # Compute pressure field on grid points
    def p_g(self):
        self.pg[self.xlo:, self.ylo:] = 0.25 \
            * (self.p[self.xlo - 1:-1, self.ylo:]
               + self.p[self.xlo - 1:-1, self.ylo - 1:-1]
               + self.p[self.xlo:, self.ylo - 1:-1]
               + self.p[self.xlo:, self.ylo:])

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

    # Shift a given range slice to left, right, top and bottom by one
    def slice(self, xl, xh, yl, yh):
        self.s_in = np.s_[xl:xh, yl:yh]
        self.s_xr = np.s_[xl + 1:xh + 1]
        self.s_xl = np.s_[xl - 1:xh - 1]
        self.s_yt = np.s_[yl + 1:yh + 1]
        self.s_yb = np.s_[yl - 1:yh - 1]

    # Impose initial condition
    def init_cond(self):
        g = self.grid
        g.u[g.xlo:g.xhi, g.ytot - 1] = g.BC['u']['t']
        g.u[g.xlo:g.xhi, g.ytot - 2] = g.BC['u']['t']

    # @profile
    def momentum(self):
        g = self.grid
        dtdx = self.dt / g.dx
        dtdxx = self.dt / (g.dx * g.dx)
        dtdy = self.dt / g.dy
        dtdyy = self.dt / (g.dy * g.dy)

        self.slice(g.xlo, g.xhi, g.ylo, g.ytot - 1)
        s_in = self.s_in
        s_xr = self.s_xr
        s_xl = self.s_xl
        s_yt = self.s_yt
        s_yb = self.s_yb

        un = g.u.copy()
        un[s_in] = g.u[s_in] \
            - dtdx * (g.uc[s_xr, s_in[1]] * g.uc[s_xr, s_in[1]]
                      - g.uc[s_in] * g.uc[s_in]) \
            - dtdy * (g.ug[s_in[0], s_yt] * g.vg[s_xr, s_in[1]]
                      - g.ug[s_in] * g.vg[s_xr, s_yb]) \
            - dtdx * (g.p[s_xr, s_in[1]] - g.p[s_in]) \
            + self.nu * (dtdxx
                         * (g.u[s_xl, s_in[1]] - 2 * g.u[s_in]
                            + g.u[s_xr, s_in[1]])
                         + dtdyy
                         * (g.u[s_in[0], s_yb] - 2 * g.u[s_in]
                            + g.u[s_in[0], s_yt]))

        self.slice(g.xlo, g.xtot - 1, g.ylo, g.yhi)
        s_in = self.s_in
        s_xr = self.s_xr
        s_xl = self.s_xl
        s_yt = self.s_yt
        s_yb = self.s_yb

        vn = g.v.copy()
        vn[s_in] = g.v[s_in] \
            - dtdx * (g.ug[s_in[0], s_yt] * g.vg[s_xr, s_in[1]]
                      - g.ug[s_xl, s_yt] * g.vg[s_in]) \
            - dtdy * (g.vc[s_in[0], s_yt] * g.vc[s_in[0], s_yt]
                      - g.vc[s_in] * g.vc[s_in]) \
            - dtdy * (g.p[s_in[0], s_yt] - g.p[s_in]) \
            + self.nu * (dtdxx
                         * (g.v[s_xl, s_in[1]] - 2 * g.v[s_in]
                            + g.v[s_xr, s_in[1]])
                         + dtdyy
                         * (g.v[s_in[0], s_yb] - 2 * g.v[s_in]
                            + g.v[s_in[0], s_yt]))
        return un, vn

    def continuity(self):
        g = self.grid
        dtdx = self.dt / g.dx
        dtdy = self.dt / g.dy

        self.slice(g.xlo, g.xtot - 1, g.ylo, g.ytot - 1)
        s_in = self.s_in
        s_xl = self.s_xl
        s_yb = self.s_yb
        dxy = self.grid.dx * self.grid.dy

        cn_err = np.zeros_like(g.p)
        cn_err[s_in] = dtdx * (g.u[s_in] - g.u[s_xl, s_in[1]]) \
            + dtdy * (g.v[s_in] - g.v[s_in[0], s_yb])

        pn = g.p.copy()
        pn[s_in] = g.p[s_in] - self.c2 * cn_err[s_in]
        return pn, np.sum(cn_err)

    # Compute L2-norm based on new and pold time steps data
    def l2norm(self, phi_n, phi_o):
        return np.sqrt(self.dtxy * np.sum((phi_n - phi_o)**2))
