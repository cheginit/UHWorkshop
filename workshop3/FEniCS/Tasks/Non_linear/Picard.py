# -*- coding: utf-8 -*-
from dolfin import *
import numpy as np
import scipy.linalg as la
#=====================================;
#  Create mesh and identify boundary  ;
#=====================================;
mesh = UnitIntervalMesh(10)

#===========================;
#  Identify the boundaries  ;
#===========================;
V = FunctionSpace(mesh,"Lagrange", 1)
def left_boundary(x, on_boundary):
    return on_boundary and near(x[0],0)
def right_boundary(x, on_boundary):
    return on_boundary and near(x[0],1)
Gamma_0 = DirichletBC(V, Constant(0.0), left_boundary)
Gamma_1 = DirichletBC(V, Constant(1.0), right_boundary)
bcs = [Gamma_0, Gamma_1]

#=================================;
# Choice of nonlinear coefficient ;
#=================================;
m=2
def q(u):
    return (1+u)**m

#==================================================;
# Define variational problem for Picard iteration  ;
#==================================================;
u = TrialFunction(V)
v = TestFunction(V)
u_k = interpolate(Constant(0.0), V)  # previous (known) u
a = inner(q(u_k)*nabla_grad(u), nabla_grad(v))*dx
f = Constant(0.0)
L = f*v*dx
#====================;
# Picard iterations  ;
#====================;
u = Function(V)
eps = 1.0 # error measure ||u-u_k||
tol = 1.0E-6  # tolerance
iter = 0   # iteration counter
maxiter = 25
file = File("Solution.pvd")
while eps > tol and iter < maxiter:
	iter += 1
	solve(a == L, u, bcs)
	diff = u.vector().array() - u_k.vector().array()
	eps = np.linalg.norm(diff, ord=np.Inf)  #l_{\infty} norm
	print "iter=%d: norm=%g" % (iter, eps) 
	u_k.assign(u) # update for next iteration
	file << u
#=================;
# Find max error  ;
#=================;
u_exact = Expression("pow((pow(2, m+1)-1)*x[0] + 1, 1.0/(m+1)) - 1", m=m, degree=5)
u_e = interpolate(u_exact, V)
diff = la.norm((u_e.vector().array() - u.vector().array()), ord=np.Inf)
print "Max error:", diff
