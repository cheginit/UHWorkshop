import numpy as np
from sys import argv, exit
from functions import Grid2D, Simulation


Re = np.float64(argv[1]) if len(argv) == 2 else 100.0
print("Re number is set to {:d}".format(int(Re)))
g = Grid2D(128, 128, 1.0)
# s = Simulation(g, cfl=0.15, c2=5.0, Re=Re)
s = Simulation(g, cfl=0.20, c2=5.8, Re=Re)

flog = open('data/residual', 'ab')

itr = 1
while True:
    g.phi_g()
    g.phi_c()

    u_n, v_n = s.momentum()

    if np.isnan(np.sum(u_n)) or np.isinf(np.sum(u_n)):
        print("Diverged after {:d} iterations.".format(itr))
        flog.close()
        exit(1)

    u_err = s.l2norm(u_n, g.u)
    v_err = s.l2norm(v_n, g.v)

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
        flog.close()
        break
    itr += 1

if itr > s.itr_max:
    print("Maximum number of iterations, {:d}, exceeded".format(itr))
    exit(2)

g.phi_g()
g.p_g()

v_x = g.ug[:, 1:]
v_y = g.vg[1:, :]
p_g = g.pg[1:, 1:]

xp = np.linspace(0.0, g.l_lid, g.ngx)
yp = np.linspace(0.0, g.l_lid, g.ngy)
xpos, ypos = np.meshgrid(xp, yp)

v_mid = 0.5 * np.sum(v_y[:,
                         np.int(v_x.shape[1] / 2) - 1:
                         np.int(v_x.shape[1] / 2) + 1], axis=1)
u_mid = 0.5 * np.sum(v_x[np.int(v_y.shape[0] / 2) - 1:
                         np.int(v_y.shape[0] / 2) + 1, :], axis=0)

np.savetxt('data/Central_U', np.c_[u_mid, yp], fmt='%.8f')
np.savetxt('data/Central_V', np.c_[v_mid, xp], fmt='%.8f')
np.savetxt('data/xyuvp', np.c_[xpos.reshape(-1), ypos.reshape(-1),
                               v_x.reshape(-1), v_y.reshape(-1),
                               p_g.reshape(-1)], fmt='%.8f')
