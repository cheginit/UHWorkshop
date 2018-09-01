from setup import savepgf, newfig
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

re = argv[1]
rc = int(argv[2])

files = {'uvp': "data/xyuvp", 'uc': "data/Central_U", 'vc': "data/Central_V",
         'ue': "data/yu", 've': "data/xv"}

data = {k: np.loadtxt(f, dtype=np.float64) for k, f in files.iteritems()}

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

x = np.unique(data['uvp'][:, 0])
y = np.unique(data['uvp'][:, 1])

u, v, p = [data['uvp'][:, i].reshape(np.shape(x)[0], -1).T for i in [2, 3, 4]]

vel = np.sqrt(u * u + v * v)

X, Y = np.meshgrid(x, y)

fig, ax = newfig(0.7)
levels = np.arange(0, 1, 0.1)
cs = ax.contour(X, Y, vel, levels=levels)
ax.clabel(cs, inline=1, fontsize=8)
ax.set_aspect('equal')
ax.tick_params(direction='out', top=False, right=False)
ax.set_title("Velocity contours for Re=" + re)
ax.set_xlabel('$x$ (m)')
ax.set_ylabel('$y$ (m)')
plt.tight_layout()
savepgf("vel_contour")

plt.clf()

fig, ax = newfig(0.7)
levels = np.arange(np.amin(p), np.amax(p), 0.05)
cs = ax.contour(X, Y, p, levels=levels)
ax.clabel(cs, inline=1, fontsize=8)
ax.set_aspect('equal')
ax.tick_params(direction='out', top=False, right=False)
ax.set_title("Pressure contours for Re=" + re)
ax.set_xlabel('$x$ (m)')
ax.set_ylabel('$y$ (m)')
plt.tight_layout()
savepgf("p_contour")

plt.clf()

skip = (slice(None, None, 5), slice(None, None, 5))
fig, ax = newfig(0.7)
ax.quiver(X[skip], Y[skip], u[skip], v[skip], vel[skip])
ax.set_aspect('equal')
ax.tick_params(direction='out', top=False, right=False)
ax.set_title("Velocity vector field")
ax.set_xlabel('$x$ (m)')
ax.set_ylabel('$y$ (m)')
plt.tight_layout()
savepgf("U_field")
