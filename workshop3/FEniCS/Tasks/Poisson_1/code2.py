from fenics import *
import numpy as np
from mpl_toolkits.mplot3d import Axes3D                                          
import matplotlib.pyplot as plt       
#=====================================;
#  Create mesh and identify boundary  ;
#=====================================;
mesh = UnitSquareMesh(8, 8)
boundaries = MeshFunction("size_t",mesh, mesh.topology().dim()-1)

#==========================;
#  Define function spaces ;
#=========================;
V = FunctionSpace(mesh, "Lagrange", 1)
V_ext = VectorFunctionSpace(mesh, "CG", 1)


n = V.dim()                                                                      
d = mesh.geometry().dim()                                                        

#===================================;
#  Define trial and test functions  ;
#===================================;
u = TrialFunction(V)
v = TestFunction(V)

#===========================;
#  Identify the boundaries  ;
#===========================;
class left(SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and near(x[0],0.0)

class right(SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and near (x[0],1.0)

leftBC  = left()
rightBC   = right()

boundaries.set_all(0)
leftBC.mark(boundaries,1)
rightBC.mark(boundaries,2)

u_left = Constant(1.0)
u_right = Expression("sin(x[1])", degree=2)
#=================================;
#  Dirichlet boundary conditions  ;
#=================================;
bc_left = DirichletBC(V,u_left,leftBC)
bc_right = DirichletBC(V, u_right,rightBC)
bc = [bc_right,bc_left]

f = Constant(-6.0)

#===========================;
#  Define variational form  ;
#===========================;
ds = Measure("ds")[boundaries]
n = FacetNormal(mesh)

a = inner(grad(u), grad(v))*dx
L = f*v*dx

#====================;
#  Compute solution  ;
#====================;
u = Function(V) 
solve(a == L, u, bc)

gradu = grad(u)
gradu_proj = project(gradu, V_ext)
#=====================================;
#  Calculate flux at production well  ;
#=====================================;
flux_form = dot(grad(u),n)*ds(2)
output_flux = assemble(flux_form)
print "output_flux = ", output_flux

#=======================================;
#  Dump solution to file in VTK format  ;
#=======================================;
file = File("mesh.pvd") 
file << mesh
file = File("u.pvd")
file << u
file = File("gu.pvd")
file << gradu_proj

VEC = np.zeros((V.dim(), 4))
counter = 0
for v in vertices(mesh):
    x = v.point().x()
    y = v.point().y()
    VEC[counter] = [x, y, gradu_proj(x,y)[0],gradu_proj(x,y)[1]]
    counter += 1
np.savetxt("Val.txt", VEC)

