"""This demo program solves Poisson's equation

    - div grad u(x, y) = f(x, y)

on the unit square with source f given by

    f(x, y) = 2*pi*pi*sin(pi*x)*sin(pi*y)

and boundary conditions given by

    u(x, y) = 0       for x = 0 or x = 1 or y=0 or y=1
    

"""


from dolfin import *

# Create mesh and define function space
mesh = UnitSquareMesh(32, 32)
boundaries = MeshFunction("size_t",mesh, mesh.topology().dim()-1)
V = FunctionSpace(mesh, "Lagrange", 1)


# Define boundary condition
class left(SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and near(x[0],0.0)

class right(SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and near (x[0],1.0)

class bottom(SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and near (x[1],0.0)

class top(SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and near (x[1],1.0)



leftBC  = left()
rightBC   = right()
bottomBC   = bottom()
topBC   = top()

boundaries.set_all(0)
leftBC.mark(boundaries,1)
rightBC.mark(boundaries,2)
bottomBC.mark(boundaries,3)
topBC.mark(boundaries,4)

bc_l = DirichletBC(V, Constant(1.0), leftBC)
bc_r = DirichletBC(V, Constant(1.0), rightBC)
bc_b = DirichletBC(V, Constant(1.0), bottomBC)
bc_t = DirichletBC(V, Constant(1.0), topBC)

bcs = [bc_l,bc_r,bc_b,bc_t]

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)

f = Expression("2*pi*pi*sin(pi*x[0])*sin(pi*x[1])", degree=2)
#g = Expression("sin(5*x[0])")
g = Constant(0.0)
a = inner(grad(u), grad(v))*dx
L = f*v*dx + g*v*ds

# Compute solution
u = Function(V)
solve(a == L, u, bcs)

# Save solution in VTK format
file = File("poisson.pvd")
file << u

