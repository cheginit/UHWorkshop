# ========================================================================== #
# Solves Navier-Stokes equations for incompressible, laminar, steady flow
# using artificial compressibility method on staggered grid.
#
# The governing equations are as follows:
# P_t + c^2 div[u] = 0
# u_t + u . grad[u] = - grad[P] + nu div[grad[u]]
#
# where P is p/rho and c represents artificial sound's speed.
#
# Lid-Driven Cavity case:
# Dimensions : 1x1 m
# Grid size  : 128 x 128
# Re number  : 100 / 1000 / 5000 / 10000
# Grid type  : Staggered Arakawa C
# Boundary Conditions: u, v -> Dirichlet (as shown below)
#                      p    -> Neumann (grad[p] = 0)
#
#                                 u=1, v=0
#                              ---------------
#                             |               |
#                         u=0 |               | u=0
#                         v=0 |               | v=0
#                             |               |
#                             |               |
#                              ---------------
#                                   u=0, v=0
# ========================================================================== #
import numpy as np
from sys import argv, exit
import functions as fn


Re = np.float64(argv[1]) if len(argv) == 2 else 100.0
print("Re number is set to {:d}".format(int(Re)))
g = fn.Grid2D(128, 128, 1.0)

if Re < 500:
    s = fn.Simulation(g, cfl=0.15, c2=5.0, Re=Re)
if Re < 2000:
    s = fn.Simulation(g, cfl=0.20, c2=5.8, Re=Re)
else:
    s = fn.Simulation(g, cfl=0.05, c2=5.8, Re=Re)

flog = open('data/residual', 'ab')

itr = 1
while True:
    # Compute velocity field on grid points and cell centers
    g.phi_g()
    g.phi_c()

    # Compute velocity field using momentum equations
    u_n, v_n = s.momentum()

    # Check if the results diverged (NAN or INF)
    if np.isnan(np.sum(u_n)) or np.isinf(np.sum(u_n)):
        print("Diverged after {:d} iterations.".format(itr))
        flog.close()
        exit(1)

    # Compute L2-norm error based on the new and old time steps
    u_err = s.l2norm(u_n, g.u)
    v_err = s.l2norm(v_n, g.v)

    # Update velocity field and apply Dirichlet boundary condition
    g.u = u_n
    g.v = v_n
    g.BC_u()
    g.BC_v()

    # Compute pressure using continuity equation and update velocity field,
    # as well as continuity error (div[v]).
    p_n, c_err = s.continuity()
    p_err = s.l2norm(p_n, g.p)
    err = np.max([u_err, v_err, p_err, c_err])

    # Write obtained errors
    np.savetxt(flog, np.c_[itr, err, u_err, v_err, p_err, c_err], fmt='%.8f')

    # Update pressure field and apply Neumman boundary condition
    g.p = p_n
    g.BC_p()

    # Check if convergence achieved
    if err < s.tol or itr > s.itr_max:
        print("Converged after {:d} iterations".format(itr))
        flog.close()
        break
    itr += 1

if itr > s.itr_max:
    print("Maximum number of iterations, {:d}, exceeded".format(itr))
    exit(2)

# Compute velocity and pressure fields on grid points for visualization
g.phi_g()
g.p_g()

# Slice the fields to exclude the ghost cells
v_x = g.ug[:, 1:]
v_y = g.vg[1:, :]
p_g = g.pg[1:, 1:]

xp = np.linspace(0.0, g.l_lid, g.ngx)
yp = np.linspace(0.0, g.l_lid, g.ngy)
xpos, ypos = np.meshgrid(xp, yp)

# Compute v at the middle of domain along x-axis
v_mid = 0.5 * np.sum(v_y[:,
                         np.int(v_x.shape[1] / 2) - 1:
                         np.int(v_x.shape[1] / 2) + 1], axis=1)
# Compute u at the middle of domain along y-axis
u_mid = 0.5 * np.sum(v_x[np.int(v_y.shape[0] / 2) - 1:
                         np.int(v_y.shape[0] / 2) + 1, :], axis=0)

# Write the field data visulization
np.savetxt('data/Central_U', np.c_[u_mid, yp], fmt='%.8f')
np.savetxt('data/Central_V', np.c_[v_mid, xp], fmt='%.8f')
np.savetxt('data/xyuvp', np.c_[xpos.reshape(-1), ypos.reshape(-1),
                               v_x.reshape(-1), v_y.reshape(-1),
                               p_g.reshape(-1)], fmt='%.8f')
