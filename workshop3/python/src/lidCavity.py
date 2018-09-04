import numpy as np
from sys import argv

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

        self.u = np.zeros([self.ngx, self.ytot])
        self.v = np.zeros([self.xtot, self.ngy])
        self.p = np.zeros([self.xtot, self.ytot])

        self.ug = np.zeros_like(self.u)
        self.vg = np.zeros_like(self.v)
        self.pg = np.zeros_like(self.p)

        self.uc = np.zeros_like(self.u)
        self.vc = np.zeros_like(self.v)

        self.BC = {'u': {'t': 1.0, 'b': 0.0, 'r': 0.0, 'l': 0.0},
                   'v': {'t': 0.0, 'b': 0.0, 'r': 0.0, 'l': 0.0},
                   'p': {'t': 0.0, 'b': 0.0, 'r': 0.0, 'l': 0.0}}

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

    def p_g(self):
        self.pg[self.xlo:, self.ylo:] = 0.25 \
            * (self.p[self.xlo - 1:-1, self.ylo:]
               + self.p[self.xlo - 1:-1, self.ylo - 1:-1]
               + self.p[self.xlo:, self.ylo - 1:-1]
               + self.p[self.xlo:, self.ylo:])

    def BC_u(self):
        self.u[self.xlo:self.xhi, 0] = \
            self.BC['u']['b'] - self.u[self.xlo:self.xhi, 1]
        self.u[self.xlo:self.xhi, self.ytot - 1] = \
            2 * self.BC['u']['t'] - self.u[self.xlo:self.xhi, self.ytot - 1]
        self.u[0, :self.ytot] = self.BC['u']['l']
        self.u[self.xhi, :self.ytot] = self.BC['u']['r']

    def BC_v(self):
        self.v[self.xlo:self.xtot - 1, 0] = self.BC['v']['b']
        self.v[self.xlo:self.xtot - 1, self.yhi] = self.BC['v']['t']
        self.v[0, :self.ngy] = \
            self.BC['v']['l'] - self.v[1, :self.ngy]
        self.v[self.xtot - 1, :self.ngy] = \
            self.BC['v']['r'] - self.v[self.xtot - 2, :self.ngy]

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

    def phi_midy(self, phi):
        return 0.5 * np.sum(phi[:,
                                np.int(phi.shape[1] / 2) - 1:
                                np.int(phi.shape[1] / 2) + 1], axis=1)

    def phi_midx(self, phi):
        return 0.5 * np.sum(phi[np.int(phi.shape[0] / 2) - 1:
                                np.int(phi.shape[0] / 2) + 1, :], axis=0)


class Simulation(object):
    def __init__(self, grid, cfl, c2, Re, tol=1e-8, itr_max=1e6):
        self.grid = grid
        self.dt = cfl * min(grid.dx, grid.dy) / grid.BC['u']['t']
        self.nu = grid.BC['u']['t'] * grid.l_lid / Re
        self.c2 = c2
        self.re = Re
        self.tol = tol
        self.itr_max = itr_max
        self.init_cond()
        grid.BC_u()
        grid.BC_v()
        grid.BC_p()

    def slice(self, xl, xh, yl, yh):
        self.s_in = np.s_[xl:xh, yl:yh]
        self.s_xr = np.s_[xl + 1:xh + 1]
        self.s_xl = np.s_[xl - 1:xh - 1]
        self.s_yt = np.s_[yl + 1:yh + 1]
        self.s_yb = np.s_[yl - 1:yh - 1]

    def init_cond(self):
        g = self.grid
        g.u[g.xlo:g.xhi, g.ytot - 1] = g.BC['u']['t']
        g.u[g.xlo:g.xhi, g.ytot - 2] = g.BC['u']['t']

    def momentum(self):
        g = self.grid

        self.slice(g.xlo, g.xhi, g.ylo, g.ytot - 1)
        s_in = self.s_in
        s_xr = self.s_xr
        s_xl = self.s_xl
        s_yt = self.s_yt
        s_yb = self.s_yb

        un = g.u.copy()
        un[s_in] = g.u[s_in] \
            - self.dt / g.dx * (g.uc[s_xr, s_in[1]] * g.uc[s_xr, s_in[1]]
                                - g.uc[s_in] * g.uc[s_in]) \
            - self.dt / g.dy * (g.ug[s_in[0], s_yt] * g.vc[s_xr, s_in[1]]
                                - g.ug[s_in] * g.vc[s_xr, s_yb]) \
            - self.dt / g.dx * (g.p[s_xr, s_in[1]] - g.p[s_in]) \
            + self.nu * (self.dt / (g.dx * g.dx)
                         * (g.u[s_xl, s_in[1]] - 2 * g.u[s_in]
                            + g.u[s_xr, s_in[1]])
                         + self.dt / (g.dy * g.dy)
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
            - self.dt / g.dx * (g.ug[s_in[0], s_yt] * g.vc[s_xr, s_in[1]]
                                - g.ug[s_xl, s_yt] * g.vc[s_in]) \
            - self.dt / g.dy * (g.vg[s_in[0], s_yt] * g.vg[s_in[0], s_yt]
                                - g.vg[s_in] * g.vg[s_in]) \
            - self.dt / g.dy * (g.p[s_in[0], s_yt] - g.p[s_in]) \
            + self.nu * (self.dt / (g.dx * g.dx)
                         * (g.v[s_xl, s_in[1]] - 2 * g.v[s_in]
                            + g.v[s_xr, s_in[1]])
                         + self.dt / (g.dy * g.dy)
                         * (g.v[s_in[0], s_yb] - 2 * g.v[s_in]
                            + g.v[s_in[0], s_yt]))
        return un, vn

    def continuity(self):
        g = self.grid
        self.slice(g.xlo, g.xtot - 1, g.ylo, g.ytot - 1)
        s_in = self.s_in
        s_xl = self.s_xl
        s_yb = self.s_yb
        dtxy = self.grid.dx * self.grid.dy * self.dt

        cn_err = np.zeros_like(g.p)
        cn_err[s_in] = self.dt / g.dx * (g.u[s_in] - g.u[s_xl, s_in[1]]) \
            + self.dt / g.dy * (g.v[s_in] - g.v[s_in[0], s_yb])

        pn = g.p.copy()
        pn[s_in] = g.p[s_in] - self.c2 * cn_err[s_in]
        return pn, dtxy * np.sum(cn_err)

    def l2norm(self, phi_n, phi_o):
        dtxy = self.grid.dx * self.grid.dy * self.dt
        return np.sqrt(dtxy * np.sum((phi_n - phi_o)**2))

Re = np.float(argv[1]) if len(argv) == 2 else 100.0
print("Re number is set to {:d}".format(int(Re)))
g = Grid2D(128, 128, 1.0)
s = Simulation(g, cfl=0.15, c2=5.0, Re=Re)

flog=open('data/residual','ab')

itr = 1
while True:
    g.phi_g()
    g.phi_c()

    u_n, v_n = s.momentum()
    u_err = s.l2norm(u_n, g.u)
    v_err = s.l2norm(v_n, g.v)

    if np.isnan(np.sum(u_n)) or np.isnan(np.sum(v_n)):
        print("Diverged after {:d} iterations.".format(itr))
        break

    g.u = u_n
    g.v = v_n
    g.BC_u()
    g.BC_v()

    p_n, c_err = s.continuity()
    p_err = s.l2norm(p_n, g.p)
    err = np.max([u_err, v_err, p_err, c_err])

    np.savetxt(flog, np.c_[itr, err, u_err, v_err, p_err, c_err], fmt='%.8f')

    g.p = p_n
    g.BC_p()

    if err < s.tol or itr > s.itr_max:
        print("Converged after {:d} iterations".format(itr))
        break
    itr += 1

if itr > s.itr_max:
    printf("Maximum number of iterations, {:d}, exceeded".format(itr))
flog.close()

g.phi_g()
g.p_g()

v_x = g.ug[:, 1:]
v_y = g.vg[1:, :]
p_g = g.pg[1:, 1:]

xp = np.linspace(0, g.l_lid, g.ngx)
yp = np.linspace(0, g.l_lid, g.ngy)
xpos, ypos = np.meshgrid(xp, yp)

np.savetxt('data/Central_U', np.c_[g.phi_midx(v_x), yp], fmt='%.8f')
np.savetxt('data/Central_V', np.c_[g.phi_midy(v_y), xp], fmt='%.8f')
np.savetxt('data/xyuvp', np.c_[xpos.reshape(-1), ypos.reshape(-1),
                               v_x.reshape(-1), v_y.reshape(-1),
                               p_g.reshape(-1)], fmt='%.8f')
