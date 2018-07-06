from sys import exit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
# from numba import autojit
# import scipy.linalg as linalg


class FDGrid1D(object):
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


class Simulation(object):
    def __init__(self, grid, u, period=5.0, scheme="FTCS", method="explicit"):
        self.grid = grid
        self.t = 0.0 # simulation time
        self.period = period # simulation time period
        self.u = u   # the constant advective velocity
        self.scheme = scheme 
        self.method = method


    # @autojit
    def advect(self, func, c):
        g = self.grid
        if self.method == "explicit":
            if self.scheme == "upwind":
                func[g.ilo:g.ihi] = func[g.ilo:g.ihi] - c*(func[g.ilo:g.ihi] -
                                func[g.ilo-1:g.ihi-1])
            elif self.scheme == "FTCS":
                func[g.ilo:g.ihi] = func[g.ilo:g.ihi] - 0.5*c*(func[g.ilo+1:g.ihi+1] -
                                func[g.ilo-1:g.ihi-1])
            else:
                exit("invalid scheme")

            g.fill_BCs()

        if self.method == "implicit":
            A = np.zeros([g.ntot, g.ntot])
            if self.scheme == "upwind":
                np.fill_diagonal(A, 1 + c)
                np.fill_diagonal(A[1:, :-1], - c)
                A[0, -1] = -c
            elif self.scheme == "FTCS":
                np.fill_diagonal(A, 1)
                np.fill_diagonal(A[1:, :-1], -0.5 * c)
                np.fill_diagonal(A[:-1, 1:], 0.5 * c)
                A[0, -1] = -0.5 * c
            else:
                exit("invalid scheme")
            # return linalg.lu_solve(linalg.lu_factor(A), x)
            # returnnp.linalg.solve(A, x) linalg.solve(A, x)
            func = np.linalg.solve(A, func)
        else:
            exit("invalid method")
        return func


    def solve(self, c=0.8):
        g = self.grid
        t = self.t 
        dt = c * g.dx / self.u
        tmax = (g.xmax-g.xmin) / self.u * self.period
        phi = g.phi.copy()
        while t < tmax:
            phi['last'] = self.advect(np.array(phi['old']), c)
            phi['old'] = phi['last']
            t += dt
        return phi


# implicit
fgrid = FDGrid1D(665, 1)
# explicit
# fgrid = FDGrid1D(65, 1)

vel = 1.0

# upwind, explicit and implicit
cycle = 1.0
# FTCS, explicit
# cycle = 0.1

s = Simulation(fgrid, vel, cycle, "upwind", "implicit")

# upwind and explicit
# cos = [0.9, 0.5, 0.1]
# FTCS and explicit
# cos = [1.0, 0.5]
# upwind, FTCS, explicit and implicit
cos = [10.0, 1.0, 0.5]

a_all = [s.solve(co) for co in cos]
solution = pd.concat(a_all, keys=cos)

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(fgrid.x, solution.loc[cos[0]]['init'], label='Exact', linestyle='--')
[ax.plot(fgrid.x, solution.loc[co]['last'], label='C = ' + str(co)) for co in cos]

plt.xlabel(r"$x$", fontsize=16)
plt.ylabel(r"$\phi$", fontsize=16)
plt.legend(frameon=False, loc="best")
plt.tight_layout()

plt.savefig("FD1DAdv.pdf")
