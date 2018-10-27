"""
This program solves the incompressible Navier-Stokes equations
for Taylor Green Vortex problem using Artificial compressibilty method.
"""

from dolfin import *
import numpy as np
import time

nx = 10
# Load mesh from file
mesh = RectangleMesh(Point(0.0,0.0),Point(2*pi,2*pi),nx,nx)
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
plot(mesh, title='mesh')

# Define function spaces (P2-P1)
uSpace = VectorFunctionSpace(mesh, "CG", 1)
pSpace = FunctionSpace(mesh, "CG", 1)

# Define trial and test functions
u = TrialFunction(uSpace)
p = TrialFunction(pSpace)
w = TestFunction(uSpace)
q = TestFunction(pSpace)

# Set parameter values
dt = 0.0001
nu = 0.2

# Define boundary conditions
#=================================;
#  Dirichlet boundary conditions  ;
#=================================;
class Right(SubDomain):
	def inside(self,x,on_boundary):
		return on_boundary and near(x[0],1)
class Left(SubDomain):
	def inside(self,x,on_boundary):
		return on_boundary and near(x[0],0)
class Top(SubDomain):
	def inside(self,x,on_boundary):
		return on_boundary and near(x[1],1)
class Bottom(SubDomain):
	def inside(self,x,on_boundary):
		return on_boundary and near(x[1],0)

class RightTop(SubDomain):
    def inside(self,x,on_boundary):
        return x[0] < DOLFIN_EPS and x[1] < DOLFIN_EPS

right = Right()
left = Left()
top = Top()
bottom = Bottom()
rightTop = RightTop()

boundaries.set_all(0)

left.mark(boundaries,1)
right.mark(boundaries,2)
bottom.mark(boundaries,3)
top.mark(boundaries,4)
rightTop.mark(boundaries,5)

bcs_left = DirichletBC(uSpace,Constant((0.0,0.0)),boundaries,1)
bcs_right = DirichletBC(uSpace,Constant((0.0,0.0)),boundaries,2)
bcs_bottom = DirichletBC(uSpace,Constant((0.0,0.0)),boundaries,3)
bcs_top = DirichletBC(uSpace,Constant((0.0,0.0)),boundaries,4)

bcs_rt = DirichletBC(pSpace,Constant(0.0),rightTop,method= "pointwise")

ubcs = [bcs_left,bcs_right,bcs_bottom,bcs_top]
pbcs = [bcs_rt]



# Create functions
u0 = Function(uSpace)
p0 = Function(pSpace)
u1 = Function(uSpace)
p1 = Function(pSpace)

# Define coefficients
k = Constant(dt)
f = Constant((0, 0))
c = Constant(2.4083)
c1 = np.sqrt((1./nx)*(1./nx)*dt)
c2 = (np.sqrt((1./nx)*(1./nx)*dt))/2.4083
tol = 1E-9

# Initial condition
u0 = interpolate(Expression(('sin(x[0])*cos(x[1])', '-sin(x[1])*cos(x[0])')) , uSpace)
p0 = interpolate(Constant(0.0), pSpace)

# BoLM step
F1 = (1./k)*inner(u - u0, w)*dx + inner(grad(u0)*u0, w)*dx + \
     nu*inner(grad(u0), grad(w))*dx\
         -inner(div(w),p0)*dx- inner(f, w)*dx


a1 = lhs(F1)
L1 = rhs(F1)

# BoM step
F2 = + (1./k)*inner(p - p0, q)*dx + \
     c*c*inner(q,div(u1))*dx

a2 = lhs(F2)
L2 = rhs(F2)


# Create files for storing solution
ufile = File("velocityTGV.pvd")
pfile = File("PressureTGV.pvd")

Initialtime = time.time()
itr = 1
set_log_active(False)
while True:
    # Velocity step
    solve(a1==L1,u1,ubcs)
    solver_parameters={'linear_solver':'gmres','preconditioner':'ilu'}

    # Pressure step
    solve(a2==L2,p1)
    solver_parameters={'linear_solver':'gmres','preconditioner':'hypre_amg'}
 
    
    # Calculate errors
    Eu_error = c1*errornorm(u1.sub(0),u0.sub(0),norm_type='L2',degree_rise= 3)
    Ev_error = c1*errornorm(u1.sub(1),u0.sub(1),norm_type='L2',degree_rise= 3)
    Ep_error = c2*errornorm(p1,p0,norm_type='L2',degree_rise= 3)
    Edelta_error =  c1*c1*assemble( div(u1) * dx )
    ##print("Eu_error is %e" % Eu_error)
    #print("Ev_error is %e" % Ev_error)
    #print("Ep_error is %e" % Ep_error)
    #print("Edelta_error is %e" % Edelta_error)
    
    #print(max(Eu_error,Ev_error,Ep_error))
    if max(Eu_error,Ev_error,Ep_error,Edelta_error)< tol:
        break
    
    # Move to next time step

    u0.assign(u1)
    p0.assign(p1)
    itr += 1
    
    if itr%100==0:
        print("itreration= %s  error norm= %E"\
          %(itr, max(Eu_error,Ev_error,Ep_error,Edelta_error)  ) )
    
    if itr%500==0:
        ufile << u1
        pfile << p1
        Solvetime = time.time()
        print("TIME = ",Solvetime-Initialtime)




# Save to file
ufile << u1
pfile << p1

