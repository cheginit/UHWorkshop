"""
Galerkin formulation for hyperelasticity (quasi-static)

Machine Problem #3

Written By: Professor Kalyana B. Nakshatrala 

Last code modification: December 3, 2014
Code last tested on January 15, 2016

"""
from dolfin import *
import numpy as np 

#=========================;
#  Select material model  ;
#=========================;
material_model = "Venant-Kirchhoff"

#========================================;
#  Create computational domain and mesh  ;
#========================================;
nx, ny = 20, 5
# mesh = RectangleMesh(0.0,0.0,1.0,0.1,nx,ny)
mesh = RectangleMesh(Point(0.0,0.0),Point(1.0,0.1),nx,ny)

#=========================;
#  Define function space  ;
#=========================;
V = VectorFunctionSpace(mesh,"Lagrange",1)

#=====================;
#  Define boundaries  ;
#=====================;
def left(x,on_boundary): return near(x[0],0.0)

# Define Dirichlet boundary (x = 0)
c = Expression(("0.0","0.0"))

bcl = DirichletBC(V, c, left)
bcs = [bcl]

# Define functions
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
u  = Function(V)                 # Displacement from previous iteration
B  = Constant((0.0,-10000.0))    # Body force per unit volume (rho = 1000)
T  = Constant((0.0,0.0))         # Traction force on the boundary

# Kinematics
I = Identity(V.cell().d)    # Identity tensor
F = I + grad(u)             # Deformation gradient
C = F.T*F                   # Right Cauchy-Green tensor
E = (C - I) / 2             # St. Venant-Kirchhoff strain 

# Invariants of deformation tensors
Ic = tr(C)
J  = det(F)

# Elasticity parameters
lmbda, mu = 1e7, 1e7 

#================================;
#  Stored strain energy density  ;
#================================;
if material_model.lower() == "neo-hookean": 
    psi = (mu/2)*(Ic - 3) - mu*ln(J) + (lmbda/2)*(ln(J))**2
elif material_model.lower() == "venant-kirchhoff": 
    psi = mu * tr(E * E) + (lmbda/2)*(tr(E))**2
else: 
    print "Error in material_model"

# Total potential energy
Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds

# Compute first variation of Pi (directional derivative about u in the direction of v)
F = derivative(Pi, u, v)

# Compute Jacobian of F
J = derivative(F, u, du)

#=============================;
#  Solve variational problem  ;
#=============================;
solve(F == 0, u, bcs, J=J)

#===============================;
#  Save solution in VTK format  ;
#===============================;
file = File("displacement.pvd");
file << u;

#========================================;
#  Obtain top tip vertical displacement  ;
#========================================; 
nodal_displacement = interpolate(u,V)
displacement_array = np.reshape(nodal_displacement.vector().array(),(mesh.num_vertices(),2))
print "top tip vertical displacement = ", displacement_array[mesh.num_vertices()-1,1]

#==========================;
#  Plot and hold solution  ;
#==========================;
plot(u, mode = "displacement", interactive = True)
