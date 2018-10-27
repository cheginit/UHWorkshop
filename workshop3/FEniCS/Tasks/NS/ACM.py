"""
This program solves the incompressible Navier-Stokes equations
for lid-driven cavity problem using Artificial compressibilty method.
"""

from dolfin import *
import numpy as np
import time

nx = 10
# Load mesh from file
mesh = RectangleMesh(Point(0.0,0.0),Point(1.0,1.0),nx,nx)
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
dt = 0.015
nu = 0.01

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
bcs_top = DirichletBC(uSpace,Constant((1.0,0.0)),boundaries,4)

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


# Save to file
ufile << u1
pfile << p1

