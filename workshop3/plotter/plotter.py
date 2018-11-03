from setup import savepgf, newfig
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sys import argv


re = argv[1]

# =========================================================================== #
uvp = np.loadtxt("data/xyuvp", dtype=np.float64)
x = np.unique(uvp[:, 0])
y = np.unique(uvp[:, 1])
u, v, p = [uvp[:, i].reshape(np.shape(x)[0], -1).T for i in [2, 3, 4]]

# Plot velocity at the middle of the domain along the x and y axes for lidCavity
if len(argv) > 2:
    rc = int(argv[2])
    um = u[:, int(u.shape[1] / 2)]
    vm = v[int(v.shape[0] / 2), :]

    ue = np.loadtxt("../../plotter/data/yu", dtype=np.float64)
    ve = np.loadtxt("../../plotter/data/xv", dtype=np.float64)

    fig, ax_u = newfig(0.8)
    ax_v = ax_u.twinx()
    ax_v2 = ax_v.twiny()

    ax_u.plot(um, y, color='g', label="Numerical, $U_x$", linewidth=0.6)
    ax_u.plot(
        ue[:, rc],
        ue[:, 0],
        '*',
        color='r',
        label="Ghia, $U_x$",
        linewidth=0.6)

    ax_v2.plot(x, vm, label="Numerical, $U_y$", linewidth=0.6)
    ax_v2.plot(
        ve[:, 0],
        ve[:, rc],
        'x',
        label="Ghia, $U_y$",
        linewidth=0.6)

    ax_u.set_xlim([-1, 1])
    ax_u.set_ylim([0, 1])
    ax_v.set_ylim([-1, 1])
    ax_v2.set_xlim([0, 1])

    ax_u.grid(alpha=0.5)
    ax_u.tick_params(direction='out', top=False, right=False)
    ax_u.set_title("Velocity profile along the middle of axes for Re=" + re)

    ax_u.set_ylabel('$y$ (m)')
    ax_v2.set_xlabel('$x$ (m)')
    ax_u.set_xlabel('$U_y$ (m.s$^{-1}$)')
    ax_v.set_ylabel('$U_x$ (m.s$^{-1}$)')

    ax_u.legend(loc='center left', bbox_to_anchor=(1.25, 0.1))
    ax_v2.legend(loc='center left', bbox_to_anchor=(1.25, 0.27))

    ax_u.set_aspect('auto')
    ax_v.set_aspect('auto')

    plt.tight_layout()
    savepgf("velocity")

    plt.clf()

# =========================================================================== #
# Extract fields data from input including X, Y, u, v and p

vel = np.sqrt(u * u + v * v)
X, Y = np.meshgrid(x, y)

# Plot pressure contours and streamlines
fig, ax = newfig(0.7)
divider = make_axes_locatable(ax)
lw = 3 * vel / vel.max()
cax = divider.append_axes('right', size='5%', pad=0.05)

cf = ax.contourf(X, Y, p, alpha=0.5, cmap=cm.viridis)
fig.colorbar(cf, cax=cax, orientation='vertical')
ax.contour(X, Y, p, cmap=cm.viridis)
ax.streamplot(X, Y, u, v, color='b', linewidth=lw, cmap='autumn')
ax.set_aspect('equal')
ax.set_xlim([X.min(), X.max()])
ax.set_ylim([Y.min(), Y.max()])
ax.tick_params(direction='out', top=False, right=False)
ax.set_title("Pressure contours and streamlines")
ax.set_xlabel('$x$ (m)')
ax.set_ylabel('$y$ (m)')
plt.tight_layout()
savepgf("pressure_stream")
plt.clf()

# =========================================================================== #
# Plot residuals
data = np.loadtxt("data/residual", dtype=np.float)
c = 10

fig, ax = newfig(0.8)
ax.semilogy(data[c:, 0], data[c:, 1], label="tot", linewidth=0.6)
ax.semilogy(data[c:, 0], data[c:, 2], label="u", linewidth=0.6)
ax.semilogy(data[c:, 0], data[c:, 3], label="v", linewidth=0.6)
ax.semilogy(data[c:, 0], data[c:, 4], label="p", linewidth=0.6)
ax.semilogy(data[c:, 0], data[c:, 5], label="div", linewidth=0.6)

ax.tick_params(direction='out', top=False, right=False)
ax.set_title("Residual")
ax.set_xlabel('Time')
ax.set_ylabel('Residual')
plt.tight_layout()
ax.legend(loc="best")
savepgf("residual")
