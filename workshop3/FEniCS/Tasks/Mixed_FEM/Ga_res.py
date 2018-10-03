"""
Governing equations: Steady-state Darcy equations
Reservoir Problem  
Formulation: Galerkin formulation 

Written By: Sarraf
"""
from dolfin import *

#=====================================;
#  Create mesh and identify boundary  ;
#=====================================;
mesh = RectangleMesh(Point(0.0,0.0),Point(2.0,1.0),50,25)
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)

#====================================================;
#  Define function spaces and mixed (product) space  ;
#====================================================;
vSpace = VectorElement("Lagrange", mesh.ufl_cell(), 2)
pSpace = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
TH = vSpace * pSpace
wSpace = FunctionSpace(mesh, TH)
#===================================;
#  Define trial and test functions  ;
#===================================;
(v,p) = TrialFunctions(wSpace)
(w,q) = TestFunctions(wSpace)

#=====================;
#  Define body force  ;
#=====================;
rhob = Constant((0.0,0.0))

#============================;
#  Define medium properties  ;
#============================; 
alpha = Constant(1.0)

#=======================;
#  Boundary conditions  ;
#=======================;
class right(SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and near (x[0],2)

class left(SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and near(x[0],0.0)

class bottom(SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and near(x[1],0.0)

class topRigid(SubDomain):              
    def inside(self,x,on_boundary):
        return on_boundary and near(x[1],1.0) and \
            ( x[0] <= 0.95 or x[0] >= 1.05 )

class topProduction(SubDomain):              
    def inside(self,x,on_boundary):
        return on_boundary and near(x[1],1.0) and \
            (x[0] > 0.95) and (x[0] < 1.05)

right = right()
left = left()
bottom = bottom()
topRigid = topRigid()
topProduction = topProduction()

boundaries.set_all(0)
left.mark(boundaries,1)
right.mark(boundaries,2)
bottom.mark(boundaries,3)
topRigid.mark(boundaries,4)
topProduction.mark(boundaries,5)

#-------------------------------------------------;
#  Prescribe normal velocity boundary conditions  ;
#-------------------------------------------------;
bc_Bottom = DirichletBC(wSpace.sub(0).sub(1),Constant(0.0),bottom)
bc_Top_Rigid = DirichletBC(wSpace.sub(0).sub(1),Constant(0.0),topRigid)

#-------------------------------;
#  Combine boundary conditions  ;
#-------------------------------;
bcs = [bc_Bottom,bc_Top_Rigid]

#===========================;
#  Define variational form  ;
#===========================;                                                                            
ds = Measure("ds")[boundaries]
n = FacetNormal(mesh)

p_Production = Constant(101325.0)
p_Injection = Constant(101325000.0)

a = dot(alpha * v,w)*dx - div(w)*p*dx - div(v)*q*dx

L = dot(rhob,w)*dx - dot(w,n) * p_Production * ds(5) - dot(w,n) * p_Injection * ds(2) - dot(w,n) * p_Injection * ds(1)

#====================;
#  Compute solution  ;
#====================;
sol = Function(wSpace)
solve(a == L, sol, bcs)
(v,p) = sol.split()

#============================================;
#  Calculate flux along the production well  ;
#============================================;             
flux = assemble( dot(v,n)*ds(5) ) 
print "production flux = ", flux  

#=======================================;
#  Dump solution to file in VTK format  ;
#=======================================;
file = File('Velocity.pvd')
file << v

file = File('Pressure.pvd')
file << p

