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
sol = Function(wSpace)
#(u,p) = TrialFunctions(wSpace)
(u,p) = (as_vector((sol[0], sol[1])), sol[2])

(w,q) = TestFunctions(wSpace)

#=====================;
#  Define body force  ;
#=====================;
f = Constant((0.0,0.0))
nu = Constant(0.01)
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

bcs_left = DirichletBC(wSpace.sub(0),Constant((0.0,0.0)),boundaries,1)
bcs_right = DirichletBC(wSpace.sub(0),Constant((0.0,0.0)),boundaries,2)
bcs_bottom = DirichletBC(wSpace.sub(0),Constant((0.0,0.0)),boundaries,3)
bcs_top = DirichletBC(wSpace.sub(0),Constant((1.0,0.0)),boundaries,4)

bcs_rt = DirichletBC(wSpace.sub(1),Constant(0.0),rightTop,method= "pointwise")

bcs = [bcs_left,bcs_right,bcs_bottom,bcs_top,bcs_rt]


#=====================================;
# Define the variational formulation  ;
#=====================================;
F=+ rho*inner(grad(u)*u, w)*dx\
+nu*inner(grad(u)+grad(u).T, grad(w))*dx\
-p*div(w)*dx\
-div(u)*q*dx-inner(f,w)*dx



#====================;
#  March over time   ;
#====================;
#  Save solution in VTK format

ufile_pvd = File("velocityMixed.pvd")
pfile_pvd = File("pressureMixed.pvd")


dw = TrialFunction(wSpace)
dF = derivative(F, sol, dw)
nsproblem = NonlinearVariationalProblem(F, sol, bcs, dF)
solver = NonlinearVariationalSolver(nsproblem)
solver.solve()
vel , pres = sol.split(True)
#(u,p) = sol.split()

ufile_pvd << vel
pfile_pvd << pres
