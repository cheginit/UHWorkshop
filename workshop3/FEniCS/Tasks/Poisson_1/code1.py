from fenics import *
#=====================================;
#  Create mesh and identify boundary  ;
#=====================================;
mesh = UnitSquareMesh(8, 8)

#==========================;
#  Define function spaces ;
#=========================;
V = FunctionSpace(mesh, "Lagrange", 1)

#===================================;
#  Define trial and test functions  ;
#===================================;
u = TrialFunction(V)
v = TestFunction(V)

#===========================;
#  Identify the boundaries  ;
#===========================;
u0 = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", degree=2)

#=================================;
#  Dirichlet boundary conditions  ;
#=================================;
bc = DirichletBC(V, u0, "on_boundary")

#===========================;
#  Define variational form  ;
#===========================;
f = Constant(-6.0)
a = inner(grad(u), grad(v))*dx
L = f*v*dx

#====================;
#  Compute solution  ;
#====================;
u = Function(V) 
solve(a == L, u, bc)

#=======================================;
#  Dump solution to file in VTK format  ;
#=======================================;
file = File("mesh.pvd") 
file << mesh
file = File("u.pvd")
file << u
