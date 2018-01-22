# Original code from https://github.com/ibackus/sod-shocktube which is
# based on http://cococubed.asu.edu/code_pages/exact_riemann.shtml

import sod
import numpy as np
import matplotlib.pyplot as plt

gamma = 1.4
dustFrac = 0.0
npts = 500
t = 0.007
Rs = 287.058
left_state =  (1e5, 348.432, Rs, 0)
right_state = (1e4, 278.746, Rs, 0)
geometry=(-5, 5, 0)

# left_state and right_state set pressure, density and u (velocity)
# geometry sets left boundary on 0., right boundary on 1 and initial
# position of the shock xi on 0.5
# t is the time evolution for which positions and states in tube should be 
# calculated
# gamma denotes specific heat
# note that gamma and npts are default parameters (1.4 and 500) in solve 
# function
positions, regions, values = sod.solve(left_state=left_state, \
    right_state=right_state, geometry=geometry, t=t, 
    gamma=gamma, npts=npts, dustFrac=dustFrac)
    
# Finally, let's plot solutions
xex = values['x']
pex = values['p']
rho = values['rho']
uex = values['u']
Tex = pex/(rho*Rs)

xnum, Tnum, unum, pnum= \
      np.loadtxt("postProcessing/singleGraph/0.007/data_T_mag(U)_p.xy" \
                 , dtype=np.float, usecols=(0,1,2,3), unpack=True)

'''
# Printing positions
print('Positions:')
for desc, vals in positions.items():
    print('{0:10} : {1}'.format(desc, vals))

# Printing p, rho and u for regions
print('Regions:')
for region, vals in sorted(regions.items()):
    print('{0:10} : {1}'.format(region, vals))
'''

fig1=plt.figure()
plt.plot(xex, pex, label="Exact", linewidth=1)
plt.plot(xnum, pnum, '--', label="Numerical", linewidth=2, color='r')
plt.title('Pressure fluctuations at t={0}'.format(t))
off = np.ptp(pnum) * 0.1
plt.ylim([ pnum.min()-off, pnum.max()+off ])
plt.xlim([-5,5])
plt.xlabel('x (m)')
plt.ylabel('p (Pa)')
plt.legend(loc="best")
fig1.savefig("plots/pressure.pdf")

fig2=plt.figure()
plt.plot(xex, uex, label="Exact", linewidth=1)
plt.plot(xnum, unum, '--', label="Numerical", linewidth=2, color='r')
plt.title('Velocity fluctuations at t={0}'.format(t))
off = np.ptp(unum) * 0.1
plt.ylim([ unum.min()-off, unum.max()+off ])
plt.xlim([-5,5])
plt.xlabel('x (m)')
plt.ylabel('U (m/s)')
plt.legend(loc="best")
fig2.savefig("plots/velocity.pdf")

fig3=plt.figure()
plt.plot(xex, Tex, label="Exact", linewidth=1)
plt.plot(xnum, Tnum, '--', label="Numerical", linewidth=2, color='r')
plt.title('Temperature fluctuations at t={0}'.format(t, gamma))
off = np.ptp(Tnum) * 0.1
plt.ylim([ Tnum.min()-off, Tnum.max()+off ])
plt.xlim([-5,5])
plt.xlabel('x (m)')
plt.ylabel('T (K)')
plt.legend(loc="best")
fig3.savefig("plots/temperature.pdf")

plt.close('all')
