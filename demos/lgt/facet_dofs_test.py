#import netgen.gui

#from lgt import *

from ngsolve import *

from ngsolve.meshes import *

from ngsolve.webgui import Draw

#import matplotlib.pyplot as plt

#from su2_yangmills import *


ne=2
mesh = MakeStructured2DMesh(quads=True, nx=ne, ny=ne, periodic_x=False, periodic_y=False)

fesT = H1(mesh, order=1)
fesL = FacetFESpace(mesh, order=0)

gfT = GridFunction(fesT)
gfT.vec[:] = 0
gfT.vec[0] = 1

Draw(gfT, mesh)

gfT = GridFunction(fesT)


gfU = GridFunction(fesU)
e_id = 0
edge = mesh.edges[e_id]
print(edge)
gfU.components[0].vec[e_id] = 1

v_of_e = edge.vertices

for v in v_of_e:
    gfT.vec[v.nr] = 1

Draw(gfT, mesh)

a = BilinearForm(fesL)
u,v = fesL.TnT()
a += u*v*dx(element_boundary=True)
#a += u*v*dx(skeleton=True)
f = LinearForm(fesL)
f += gfT*v*dx(element_boundary=True)
a.Assemble()
f.Assemble()
ainv = a.mat.Inverse(fesL.FreeDofs())

gfsol = GridFunction(fesL)
gfsol.vec.data = ainv*f.vec

#Draw(gfsol, mesh)
print(gfsol.vec)
#for e in mesh.edges:
#    #print(e.nr)
#    print(gfU.components[0].vec[e.nr])
