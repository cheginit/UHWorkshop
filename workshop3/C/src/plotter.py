from setup import savepgf, newfig
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

re = argv[1]
rc = int(argv[2])

files = ["data/xyuvp", "data/Central_U", "data/Central_V",
         "data/yu", "data/xv"]

data = []
[data.append(np.loadtxt(f, dtype=np.float64)) for f in files]

fig, ax = newfig(0.7)
ax.plot(data[1][:, 0], data[1][:, 1],
        label="Numerical", linewidth=0.6)
ax.plot(data[3][:, rc], data[3][:, 0], 'x',
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
ax.plot(data[2][:, 0], data[2][:, 1],
        label="Numerical", linewidth=0.6)
ax.plot(data[4][:, rc], data[4][:, 0], 'x',
        label="Experimental", linewidth=0.6)
ax.grid(alpha=0.5)
ax.tick_params(direction='out', top=False, right=False)
ax.set_title("Velocity along the middle of $y$-axis at Re=" + re)
ax.set_xlabel(r'$V$ (m s$^{-1}$)')
ax.set_ylabel('$x$ (m)')
ax.legend(loc="best")
plt.tight_layout()
savepgf("vvelocity")

plt.clf()

x = np.unique(data[0][:, 0])
y = np.unique(data[0][:, 1])

u, v, p = [data[0][:, i].reshape(np.shape(x)[0], -1) for i in [2, 3, 4]]

vel = np.zeros_like(u)
vel[:] = np.sqrt(np.power(u[:], 2) + np.power(v[:], 2))

X, Y = np.meshgrid(x, y)

fig, ax = newfig(0.7)
levels = np.arange(0, 1, 0.1)
cs = ax.contour(Y, X, vel, levels=levels)
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
cs = ax.contour(Y, X, p, levels=levels)
ax.clabel(cs, inline=1, fontsize=8)
ax.set_aspect('equal')
ax.tick_params(direction='out', top=False, right=False)
ax.set_title("Pressure contours for Re=" + re)
ax.set_xlabel('$x$ (m)')
ax.set_ylabel('$y$ (m)')
plt.tight_layout()
savepgf("p_contour")
