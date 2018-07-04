from sys import exit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
# from numba import autojit
# import scipy.linalg as linalg



class FDGrid:

    def __init__(self, nx, no, xmin=0.0, xmax=1.0):

        self.xmin = xmin
        self.xmax = xmax
        self.no = no
        self.nx = nx

        self.ilo = no
        self.ihi = no + nx - 1
        self.ntot = nx + 2*no

        self.dx = (xmax-xmin) / nx
        self.x = xmin + (np.arange(self.ntot)-no) * self.dx

        self.phi = pd.DataFrame(self.initialize(),
                                columns=['init', 'old', 'last'])

    def initialize(self):
        x_init = np.zeros([self.ntot, 3])
        x_init[int(1. / (3*self.dx)):int(2. / (3*self.dx)), :2] = 1
        return x_init

    def fill_BCs(self):
        self.phi['last'][self.ilo-1] = self.phi['last'][self.ihi-1]
        self.phi['last'][self.ihi+1] = self.phi['last'][self.ilo+1]


# @autojit
def advect_explicit(x, c, scheme="upwind"):
    if scheme == "upwind":
        x[g.ilo:g.ihi] = x[g.ilo:g.ihi] - c*(x[g.ilo:g.ihi] -
                         x[g.ilo-1:g.ihi-1])
    elif scheme == "FTCS":
        x[g.ilo:g.ihi] = x[g.ilo:g.ihi] - 0.5*c*(x[g.ilo+1:g.ihi+1] -
                         x[g.ilo-1:g.ihi-1])
    else:
        exit("invalid method")

    g.fill_BCs()

    return x


# @autojit
def advect_implicit(x, c, scheme="upwind"):
    A = np.zeros([g.ntot, g.ntot])
    if scheme == "upwind":
        np.fill_diagonal(A, 1 + c)
        np.fill_diagonal(A[1:, :-1], -c)
        A[0, -1] = -c
    elif scheme == "FTCS":
        np.fill_diagonal(A, 1)
        np.fill_diagonal(A[1:, :-1], -0.5 * c)
        np.fill_diagonal(A[:-1, 1:], 0.5 * c)
        A[0, -1] = -0.5 * c
    else:
        exit("invalid method")

    # return linalg.lu_solve(linalg.lu_factor(A), x)
    # return linalg.solve(A, x)
    return np.linalg.solve(A, x)


def solve(x, c, scheme="upwind", solution="explicit"):
    dt = c * g.dx / u
    tmax = (g.xmax-g.xmin) / u * cycle
    t = 0.0
    while t < tmax:
        if solution == "explicit":
            x['last'] = advect_explicit(np.array(x['old']), c, scheme)
        elif solution == "implicit":
            x['last'] = advect_implicit(np.array(x['old']), c, scheme)
        x['old'] = x['last']
        t += dt
    return x


# implicit
g = FDGrid(665, 1)
# explicit
# g = FDGrid(65, 1)

u = 1.0

# upwind, explicit and implicit
cycle = 1
# FTCS, explicit
# cycle = 0.1

# upwind and explicit
# cos = [0.9, 0.5, 0.1]
# FTCS and explicit
# cos = [1.0, 0.5]
# upwind, FTCS, explicit and implicit
cos = [10.0, 1.0, 0.5]

a_all = [solve(g.phi.copy(), co, 'upwind', 'implicit') for co in cos]
solution = pd.concat(a_all, keys=cos)

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(g.x, solution.loc[cos[0]]['init'], label='Exact', linestyle='--')
[ax.plot(g.x, solution.loc[co]['last'], label='C = ' + str(co)) for co in cos]

plt.xlabel(r"$x$", fontsize=16)
plt.ylabel(r"$\phi$", fontsize=16)
plt.legend(frameon=False, loc="best")
plt.tight_layout()

plt.savefig("FD1DAdv.pdf")
