from mshr import *
from dolfin import *
from numpy import cos , pi

mesh_1 = UnitSquareMesh(10, 10, "right/left")
file = File("mesh_1.pvd")
file << mesh_1
#
def mesh2(lx,ly, Nx,Ny):
    m = UnitSquareMesh(Nx, Ny)
    x = m.coordinates()

    #Refine on top and bottom sides
    x[:,1] = (x[:,1] - 0.5) * 2.
    x[:,1] = 0.5 * (cos(pi * (x[:,1] - 1.) / 2.) + 1.)   

    #Scale
    x[:,0] = x[:,0]*lx
    x[:,1] = x[:,1]*ly

    return m

mesh_2 = mesh2(4.0,1.0, 15,15)
file = File("mesh_2.pvd")
file << mesh_2
#
# Generating geometry
ret1 = Rectangle(Point(0,0), Point(3,3))
circ = Circle(Point(1.5,1.5),0.5)
point = Point(0.5,1)
# Define domain and resolution
domain1 = ret1 - circ
res = 50
# Generate mesh
mesh_3 = generate_mesh(domain1, res)
file = File("mesh_3.pvd")
file << mesh_3



