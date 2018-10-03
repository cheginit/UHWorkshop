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

#=====================;
#  Define parameters  ;
#=====================;
alpha = 3; beta = 1.2
time = 0.0
g = Expression("1 + x[0]*x[0] + \
alpha*x[1]*x[1] + beta*time", alpha=alpha, beta=beta, time=time, degree=2)
#=================================;
#  Dirichlet boundary conditions  ;
#=================================;
bcs = DirichletBC(V,g,"on_boundary")

# Decide on a time step
dt = 0.3
u0 = interpolate(g, V)

#=====================================;
# Define the variational formulation  ;
#=====================================;
f = Constant(beta - 2 - 2*alpha)
a = u*v*dx + dt*inner(grad(u), grad(v))*dx
L = u0*v*dx + dt*f*v*dx  # updates at each time step


#====================;
#  March over time   ;
#====================;
file = File("Solution.pvd")
u = Function(V)
T = 2 # Set end time 
time = dt
while time <= T:
	solve(a==L,u,bcs)
	time += dt
	g.time=time
	u0.assign(u) # u0 := u
	file << u
#=======================================;
#  Dump solution to file in VTK format  ;
#=======================================;
file = File("mesh.pvd") 
file << mesh
