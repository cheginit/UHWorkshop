from fenics import *
#=====================================;
#  Create mesh and identify boundary  ;
#=====================================;
mesh = UnitSquareMesh(2,2)
#==========================;
#  Define function spaces ;
#=========================;
V = FunctionSpace(mesh, "Lagrange", 1)

#===================================;
#  Define trial and test functions  ;
#===================================;
f = Expression("2*pi*pi*sin(pi*x[0])*sin(pi*x[1])", degree=2)
u = TrialFunction(V)
v = TestFunction(V)

#==============;
#  Marking Bcs ;
#==============;
u0 = Constant(0.0)
bc = DirichletBC(V, u0, "on_boundary")

#===========================;
#  Define variational form  ;
#===========================;
a = inner(grad(u), grad(v))*dx
L = f*v*dx

#====================;
#  Compute solution  ;
#====================;
sol = Function(V) 
solve(a == L, sol, bc)
gradu = grad(sol)
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
print("H1_error is %e" % H1_error)

