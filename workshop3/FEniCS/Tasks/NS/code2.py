from fenics import *
from mshr import *
#=====================================;
#  Create mesh and identify boundary  ;
#=====================================;
rec1 = Rectangle(Point(0,0), Point(1,1))
rec2 = Rectangle(Point(0.45,0.45), Point(0.55,0.55))
# Define domain and resolution
domain1 = rec1 - rec2
res = 50
# Generate mesh
mesh = generate_mesh(domain1, res)
boundaries = MeshFunction("size_t",mesh, mesh.topology().dim()-1)

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
class Outer(SubDomain):
	def inside(self,x,on_boundary):
		return on_boundary and near(x[0],0) or near(x[0],1) or  near(x[1],0) or near(x[1],1)
class Inner(SubDomain):
	def inside(self,x,on_boundary):
		return on_boundary and near(x[0],0.45) or  near(x[0],0.55) or  near(x[1],0.45) or  near(x[1],0.55)

inners = Inner()
outers = Outer()

boundaries.set_all(0)
inners.mark(boundaries,1)
outers.mark(boundaries,2)

bcs_outers = DirichletBC(V,g,inners)
bcs_inners = DirichletBC(V,Constant(0.0),outers)

bcs = [bcs_outers,bcs_inners]

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
	print '=============================='
   	print '          time =', time
    	print '=============================='
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
