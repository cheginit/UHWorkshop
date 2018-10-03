from fenics import *
#=====================================;
#  Create mesh and identify boundary  ;
#=====================================;
mesh = RectangleMesh(Point(0.0,0.0),Point(1.0,1.0),10,10)
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
plot(mesh, title='mesh')
#====================================================;
#  Define function spaces and mixed (product) space  ;
#====================================================;
uSpace = VectorFunctionSpace(mesh,"CG", 3)
pSpace = FunctionSpace(mesh,"CG", 1)
wSpace = MixedFunctionSpace([uSpace,pSpace])
#===================================;
#  Define trial and test functions  ;
#===================================;
(u,p) = TrialFunctions(wSpace)
(w,q) = TestFunctions(wSpace)

#=====================;
#  Define body force  ;
#=====================;
f = Constant((0.0,0.0))
mu = Constant(1.0)
rho = Constant(1.0)

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

right = Right()
left = Left()
top = Top()
bottom = Bottom()

boundaries.set_all(0)

left.mark(boundaries,1)
right.mark(boundaries,2)
bottom.mark(boundaries,3)
top.mark(boundaries,4)

bcs_left = DirichletBC(wSpace.sub(0),Constant((10.0,0.0)),boundaries,1)
bcs_right = DirichletBC(wSpace.sub(0),Constant((-10.0,0.0)),boundaries,2)

bcs = [bcs_left,bcs_right]
ds = Measure("ds")[boundaries]

# Traction BCs
t_up = Constant((0.0,0.0))
t_down = Constant((0.0,0.0))
#=====================;
#  Define parameters  ;
#=====================;
dt = 0.01
T = 0.3 # Set end time
#  Initial condition  
u0 = interpolate(Constant((0.0,0.0)),uSpace)
#=====================================;
# Define the variational formulation  ;
#=====================================;
a = rho*dot(u,w)*dx\
+ dt*rho*inner(grad(u0)*u0, w)*dx\
+dt*2*mu*inner(grad(u)+grad(u).T, grad(w))*dx\
-dt *p*div(w)*dx\
-dt*div(u)*q*dx
L = dot(u0,w)*dx\
+dt*inner(f,w)*dx\
+dt*inner(t_down,w) * ds(3)+ dt*inner(t_up,w) * ds(4)

F= a-L 

#====================;
#  March over time   ;
#====================;
#  Save solution in VTK format

ufile_pvd = File("velocity.pvd")
pfile_pvd = File("pressure.pvd")

sol = Function(wSpace)
time = 0.0
while time <= T:
    print '=============================='
    print '          time =', time
    print '=============================='
    solve(F==0,sol,bcs, solver_parameters={"newton_solver":{"relative_tolerance": 1e-6}})
    vel , pres = sol.split(True)
    u0.assign(vel) # u0 := u
    time += dt
    ufile_pvd << vel
    pfile_pvd << pres
