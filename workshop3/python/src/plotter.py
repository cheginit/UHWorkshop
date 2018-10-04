from setup import savepgf, newfig
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sys import argv


# Read field data
re = argv[1]
rc = int(argv[2])

files = {'uvp': "data/xyuvp", 'uc': "data/Central_U", 'vc': "data/Central_V",
         'ue': "data/yu", 've': "data/xv"}

data = {k: np.loadtxt(f, dtype=np.float64) for k, f in files.items()}

# =========================================================================== #
# Plot velocity at the middle of the domain along the x and y axes
fig, ax_u = newfig(0.8)
ax_v = ax_u.twiny().twinx()

ax_u.plot(data['uc'][:, 0], data['uc'][:, 1], color= 'g',
        label="Numerical, $U_x$", linewidth=0.6)
ax_u.plot(data['ue'][:, rc], data['ue'][:, 0], '*', color= 'r',
        label="Ghia, $U_x$", linewidth=0.6)

ax_v.plot(data['vc'][:, 1], data['vc'][:, 0],
        label="Numerical, $U_y$", linewidth=0.6)
ax_v.plot(data['ve'][:, 0], data['ve'][:, rc], 'x',
        label="Ghia, $U_y$", linewidth=0.6)

ax_u.grid(alpha=0.5)
ax_u.tick_params(direction='out', top=False, right=False)
ax_u.set_title("Velocity profile along the middle of axes for Re=" + re)
ax_v.set_xlabel('$x$ (m)')
ax_v.set_ylabel('$U_y$ (m s$^{-1}$)')
ax_u.set_ylabel('$y$ (m)')
ax_u.set_xlabel('$U_x$ (m s$^{-1}$)')
ax_u.legend(loc='center left', bbox_to_anchor=(1.2, 0.1))
ax_v.legend(loc='center left', bbox_to_anchor=(1.2, 0.27))

x0, x1 = ax_u.get_xlim()
y0, y1 = ax_u.get_ylim()
ax_u.set_aspect((x1-x0)/(y1-y0))

x0, x1 = ax_v.get_xlim()
y0, y1 = ax_v.get_ylim()
ax_v.set_aspect((x1-x0)/(y1-y0))

plt.tight_layout()
savepgf("velocity")

plt.clf()

# =========================================================================== #
# Extract fields data from input including X, Y, u, v and p
x = np.unique(data['uvp'][:, 0])
y = np.unique(data['uvp'][:, 1])

u, v, p = [data['uvp'][:, i].reshape(np.shape(x)[0], -1).T for i in [2, 3, 4]]
vel = np.sqrt(u * u + v * v)
X, Y = np.meshgrid(x, y)

# Plot pressure contours and streamlines
fig, ax = newfig(0.7)
divider = make_axes_locatable(ax)
lw = 3*vel / vel.max()
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
