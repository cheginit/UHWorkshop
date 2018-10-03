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

def Dq(u):
    return m*(1+u)**(m-1)
#=============================;
# Define variational problem  ;
#=============================;
v  = TestFunction(V)
du = TrialFunction(V)

u_ = Function(V)  # most recently computed solution

F  = inner(q(u_)*nabla_grad(u_), nabla_grad(v))*dx
#J = derivative(F, u_, du) # J is Gateaux derivative at u_ in the direction of du

#====================;
# Compute solution   ;
#====================;
#problem = NonlinearVariationalProblem(F, u_, bcs, J)
#solver  = NonlinearVariationalSolver(problem)
#solver.solve()
solve(F==0, u_, bcs)

file = File("Newton.pvd")
file << u_
#=================;
# Find max error  ;
#=================;
u_exact = Expression("pow((pow(2, m+1)-1)*x[0] + 1, 1.0/(m+1)) - 1", m=m, degree=5)
u_e = interpolate(u_exact, V)
diff = la.norm((u_e.vector().array() - u_.vector().array()), ord=np.Inf)
print "Max error:", diff
