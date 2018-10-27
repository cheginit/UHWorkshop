"""
This program solves the incompressible Navier-Stokes equations
for lid-driven cavity problem using Chorin's splitting method.
"""

from dolfin import *

# Load mesh from file
mesh = RectangleMesh(Point(0.0,0.0),Point(1.0,1.0),50,50)
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
plot(mesh, title='mesh')

# Define function spaces (P2-P1)
uSpace = VectorFunctionSpace(mesh, "CG", 2)
pSpace = FunctionSpace(mesh, "CG", 1)

# Define trial and test functions
u = TrialFunction(uSpace)
p = TrialFunction(pSpace)
v = TestFunction(uSpace)
q = TestFunction(pSpace)

# Set parameter values
dt = 0.0015
T = 10.0
nu = 0.001

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

bcs = [bcs_left,bcs_right,bcs_bottom,bcs_top,bcs_rt]



# Create functions
u0 = Function(uSpace)
u1 = Function(uSpace)
p1 = Function(pSpace)

# Define coefficients
k = Constant(dt)
f = Constant((0, 0))

# Tentative velocity step
F1 = (1/k)*inner(u - u0, v)*dx + inner(grad(u0)*u0, v)*dx + \
     nu*inner(grad(u), grad(v))*dx - inner(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Pressure update
a2 = inner(grad(p), grad(q))*dx
L2 = -(1/k)*div(u1)*q*dx

# Velocity update
a3 = inner(u, v)*dx
L3 = inner(u1, v)*dx - k*inner(grad(p1), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Create files for storing solution
ufile = File("velocity.pvd")
pfile = File("Pressure.pvd")

# Time-stepping
t = dt
p = Progress("Time-stepping")
while t < T + DOLFIN_EPS:

    # Compute tentative velocity step
    begin("Computing tentative velocity")
    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in bcs]
    solve(A1, u1.vector(), b1, "gmres", "ilu")
    end()

    # Pressure correction
    begin("Computing pressure correction")
    b2 = assemble(L2)
    solve(A2, p1.vector(), b2, "gmres", "hypre_amg")
    end()

    # Velocity correction
    begin("Computing velocity correction")
    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in bcs]
    solve(A3, u1.vector(), b3, "gmres", "ilu")
    end()

    # Plot solution
    #plot(p1, title="Pressure", rescale=True)
    #plot(u1, title="Velocity", rescale=True)

    # Save to file
    ufile << u1
    pfile << p1

    # Move to next time step
    u0.assign(u1)
    p.update(t / T)
    t += dt

# Hold plot
#interactive()
