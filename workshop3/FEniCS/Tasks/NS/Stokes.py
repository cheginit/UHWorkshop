from dolfin import *
#===========================;
#  Load computational mesh  ;
#===========================;
mesh = RectangleMesh(Point(0.0,0.0),Point(1.0,1.0),10,10)
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)

#==========================;
#  Define function spaces  ;
#==========================;
velSpace = VectorFunctionSpace(mesh,"CG",3)
pSpace   = FunctionSpace(mesh,"CG",1)
wSpace   = MixedFunctionSpace([velSpace,pSpace])

#==============================;
#  Define variational problem  ;
#==============================;
(v,p) = TrialFunctions(wSpace)
(tau,q) = TestFunctions(wSpace)

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


# Define new measures associated with the interior domains and
# exterior boundaries
ds = Measure("ds")[boundaries]

#========================;
#  Prescribed IBVP data  ;
#========================;

#  Define forcing function  
rhob = Expression(("0.0","0.0"))

#  Initial condition  
u0 = interpolate(Expression(("0.0","0.0")),velSpace)

#  Pressure BCs 
pout_t = Constant((0.0,0.0))

pout_b = Constant((0.0,0.0))

#==================================;
#  Time discretization parameters  ;
#==================================;
deltat, TotalTime = 0.01, 0.3
dt = Constant(deltat)



a = dot(v,tau) * dx + dt * ( inner(grad(v) + (grad(v)).T,grad(tau))*dx - div(tau)*p*dx - q*div(v)*dx )

L = dt * ( inner(rhob, tau)*dx + inner(pout_b,tau) * ds(3) +\
    inner(pout_t,tau) * ds(4) ) + dot(tau,u0) * dx

#===============================;
#  Save solution in VTK format  ;
#===============================;
ufile_pvd = File("Pul_Vel_velocity.pvd")
pfile_pvd = File("Pul_Vel_pressure.pvd")

#===================;
#  March over time  ;
#===================;
time =0.0
U = Function(wSpace)
while time < TotalTime:
    solve(a == L, U, bcs)
    vel, pres = U.split(True)
    u0.assign(vel)
    ufile_pvd << vel
    pfile_pvd << pres
    time += deltat

#=============================;
#  Plot interactive solution  ;
#=============================;
#plot(vel,title="Velocity")
#plot(pres,title="Pressure")
#interactive()

