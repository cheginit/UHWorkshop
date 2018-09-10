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
# Plot v component  of velocit at the middle of the domain along the x-dir
fig, ax = newfig(0.7)
ax.plot(data['uc'][:, 0], data['uc'][:, 1],
        label="Numerical", linewidth=0.6)
ax.plot(data['ue'][:, rc], data['ue'][:, 0], 'x',
        label="Experimental", linewidth=0.6)
ax.grid(alpha=0.5)
ax.tick_params(direction='out', top=False, right=False)
ax.set_title("Velocity along the middle of $x$-axis at Re=" + re)
ax.set_xlabel(r'$U$ (m s$^{-1}$)')
ax.set_ylabel('$y$ (m)')
ax.legend(loc="best")
plt.tight_layout()
savepgf("uvelocity")

plt.clf()

# =========================================================================== #
# Plot v component  of velocit at the middle of the domain along the x-dir
fig, ax = newfig(0.7)
ax.plot(data['vc'][:, 1], data['vc'][:, 0],
        label="Numerical", linewidth=0.6)
ax.plot(data['ve'][:, 0], data['ve'][:, rc], 'x',
        label="Experimental", linewidth=0.6)
ax.grid(alpha=0.5)
ax.tick_params(direction='out', top=False, right=False)
ax.set_title("Velocity along the middle of $y$-axis at Re=" + re)
ax.set_xlabel('$x$ (m)')
ax.set_ylabel(r'$V$ (m s$^{-1}$)')
ax.legend(loc="best")
plt.tight_layout()
savepgf("vvelocity")

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
ax.set_title("Velocity vector field")
ax.set_xlabel('$x$ (m)')
ax.set_ylabel('$y$ (m)')
plt.tight_layout()
savepgf("U_field")
