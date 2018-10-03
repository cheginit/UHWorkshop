from fenics import *
import numpy as np
import math,sys,time,copy
nx = int(sys.argv[1])
#text = sys.argv[2]

mesh = UnitSquareMesh(nx,nx)
V = FunctionSpace(mesh, "Lagrange", 1)
#--dof--#

print "Degree-of-freedom is", V.dim()
print "nxLOG =",np.log10(nx)
f = Expression("2*pi*pi*sin(pi*x[0])*sin(pi*x[1])", degree=2)
u = TrialFunction(V)
v = TestFunction(V)

#==============;
#  Marking Bcs ;
#==============;
u0 = Constant(0.0)
bc = DirichletBC(V, u0, "on_boundary")

#==============;
#   Weak form  ;
#==============;
a = inner(grad(u), grad(v))*dx
L = f*v*dx

#==============;
#   Solution   ;
#==============;
initialtime = time.time()

A = assemble(a)
b = assemble(L)
np.savetxt("K_mat.dat", A.array())
print("K_mat: "), A.array()
np.savetxt("load_vec.dat", b.array())
print("load_vec: "), b.array()
assembletime = time.time()

sol = Function(V)
solve(a == L, sol, bc)

solvetime = time.time()
totaltime = solvetime - initialtime
print("Assembly time = %1.3e seconds" % (assembletime - initialtime))
print("Solve time = %1.3e seconds" % (solvetime - assembletime))
print("Totaltime = %1.3e seconds" % totaltime)
#=====================;
#   Dumping results   ;
#=====================;
file = File("mesh.pvd") 
file << mesh
file = File("sol.pvd")
file << sol

#=================;
#   Exact solution;
#=================;
u_ex = Function(V)
u_Exact = Expression("sin(pi*x[0])*sin(pi*x[1])",degree=5)
u_ex = interpolate(u_Exact,V)
L2_error = errornorm(u_ex,sol,norm_type='L2',degree_rise= 3)
H1_error = errornorm(u_ex,sol,norm_type='H1',degree_rise= 3)
print("L2_error is %e" % L2_error)
print("L2_errorLOG is %e" % np.log10(L2_error))
print("H1_error is %e" % H1_error)

