#import netgen.gui

#from lgt import *

from ngsolve import *

from ngsolve.meshes import *

from ngsolve.webgui import Draw

#import matplotlib.pyplot as plt

#from su2_yangmills import *


ne=2
mesh = MakeStructured2DMesh(quads=True, nx=ne, ny=ne, periodic_x=False, periodic_y=False)

fesW = H1(mesh, order=1)

#Draw(mesh)

#for v in mesh.vertices:
#    print (v, v.point, v.edges)
#    vedges = []
#    for edge_id in v.edges:
#        #print(edge_id)
#        edge = mesh[edge_id]
#        if v in edge.vertices:
#            vedges.append(edge)
#    print (v, v.point, vedges)

    #print (v, v.point, v.faces)


gfL = GridFunction(fesW)

gfL.vec[:] = 0

gfL.vec[0] = 1

Draw(gfL, mesh)

print("\n")

#for e in mesh.edges:
#    print (e, e.vertices)

#print("\n")

#for f in mesh.faces:
#    print (f, f.edges)

#print("\n")

input()
