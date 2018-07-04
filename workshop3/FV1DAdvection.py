from sys import exit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'


class Grid1D:

    def __init__(self, nx, no, xmin=0.0, xmax=1.0):

        self.xmin = xmin
        self.xmax = xmax
        self.no = no
        self.nx = nx

        self.ilo = no
        self.ihi = no + nx - 1
        self.ntot = nx + 2*no

        self.dx = (xmax-xmin) / nx
        self.x = xmin + (np.arange(self.ntot) - no)*self.dx

        self.phi = pd.DataFrame(self.initialize(),
                                columns=['init', 'old', 'last', 'slope'])
        self.sf = pd.DataFrame(np.zeros([self.ntot + 1, 2]),
                               columns=['L', 'R'])

    def initialize(self):
        x_init = np.zeros([self.ntot, 4])
        x_init[int(1. / (3*self.dx)):int(2. / (3*self.dx)), :2] = 1.0
        return x_init


def fill_BCs(x):
    for n in range(g.no):
        x[g.ilo-1-n] = x[g.ihi-n]
        x[g.ihi+1+n] = x[g.ilo+n]
    return x


def minmod(a, b):
    c = a.copy()
    for i in range(a.shape[0]):
        if a[i]*b[i] < 0.0:
            c[i] = 0.0
        elif abs(a[i]) < abs(b[i]):
            c[i] = a[i]
        elif abs(b[i]) < abs(a[i]):
            c[i] = b[i]

    return c


def slope(x, method='centered'):
    if method == 'centered':
        g.phi['slope'][g.ilo-1:g.ihi+2] = 0.5 * (x[g.ilo:g.ihi+3]
                                                 - x[g.ilo-2:g.ihi+1])/g.dx
    if method == 'minmod' or method == 'MC':
        g.phi['slope'][g.ilo-1:g.ihi+2] = (x[g.ilo-1:g.ihi+2]
                                           - x[g.ilo-2:g.ihi+1])/g.dx
        a = np.array(g.phi['slope'])
        g.phi['slope'][g.ilo-1:g.ihi+2] = (x[g.ilo:g.ihi+3]
                                           - x[g.ilo-1:g.ihi+2])/g.dx
        b = np.array(g.phi['slope'])

        if method == 'MC':
            g.phi['slope'][g.ilo-1:g.ihi+2] = (x[g.ilo:g.ihi+3]
                                               - x[g.ilo-2:g.ihi+1])/g.dx
            d = np.array(g.phi['slope'])

            g.phi['slope'] = minmod(minmod(2*a, 2*b), 0.5*d)
        elif method == 'minmod':
            g.phi['slope'] = minmod(a, b)

    return np.array(g.phi['slope'])


def flux(x, c, method='centered'):
    g.sf['L'][g.ilo:g.ihi+2] = x[g.ilo-1:g.ihi+1] + \
        0.5 * g.dx * (1.0-c) * slope(x, method)[g.ilo-1:g.ihi+1]
    g.sf['R'][g.ilo:g.ihi+2] = x[g.ilo:g.ihi+2] - \
        0.5 * g.dx * (1.0+c) * slope(x, method)[g.ilo:g.ihi+2]
    side = 'L' if u > 0 else 'R'
    return np.array(g.sf[side])


def advect_explicit(x, c, scheme="predictor-corrector"):
    if scheme == "predictor-corrector":
        x[g.ilo:g.ihi+1] = x[g.ilo:g.ihi+1] -\
                           c * (flux(x, c, 'MC')[g.ilo+1:g.ihi+2]
                                - flux(x, c, 'MC')[g.ilo:g.ihi+1])
    elif scheme == "FTCS":
        x[g.ilo:g.ihi] = x[g.ilo:g.ihi] - 0.5*c*(x[g.ilo+1:g.ihi+1]
                                                 - x[g.ilo-1:g.ihi-1])
    else:
        exit("invalid method")

    return x


def advect_implicit(x, c, scheme="predictor-corrector"):
    A = np.zeros([g.ntot, g.ntot])
    if scheme == "upwind":
        np.fill_diagonal(A, 1+c)
        np.fill_diagonal(A[1:, :-1], -c)
        A[0, -1] = -c
    elif scheme == "FTCS":
        np.fill_diagonal(A, 1)
        np.fill_diagonal(A[1:, :-1], -0.5*c)
        np.fill_diagonal(A[:-1, 1:], 0.5*c)
        A[0, -1] = -0.5*c
    else:
        exit("invalid method")

    return np.linalg.solve(A, x)


def solve(x, c, scheme="predictor-corrector", solution="explicit"):
    dt = c*g.dx/u
    tmax = (g.xmax - g.xmin)/u*cycle
    t = 0.0
    while t < tmax:
        if solution == "explicit":
            x['last'] = advect_explicit(np.array(x['old']), c, scheme)
        elif solution == "implicit":
            x['last'] = advect_implicit(np.array(x['old']), c, scheme)
        x['old'] = fill_BCs(x['last'])
        t += dt
    return x


g = Grid1D(64, 2)
u = 5.0
# upwind, explicit and implicit
cycle = 5
# FTCS, explicit
# cycle = 0.1

# upwind and explicit
cos = [0.7]
# FTCS and explicit
# cos = [1.0, 0.5]

a_all = [solve(g.phi.copy(), co, 'predictor-corrector', 'explicit')
         for co in cos]
solution = pd.concat(a_all, keys=cos)

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(g.x, solution.loc[cos[0]]['init'], label='Exact', linestyle='--')
[ax.plot(g.x, solution.loc[co]['last'], label='C = ' + str(co))
 for co in cos]

plt.xlabel(r"$x$", fontsize=16)
plt.ylabel(r"$\phi$", fontsize=16)
plt.legend(frameon=False, loc="best")
plt.tight_layout()

plt.savefig("FV1D.pdf")
