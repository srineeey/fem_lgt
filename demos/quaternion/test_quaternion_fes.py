from ngsolve.bla import Vec3D
from liegroupfem import *
from math import *

q1 = Quaternion(0, (1,0,0))
q2 = Quaternion(1, (0,1,0))



print ("q1 =", q1)
print ("q2 =", q2)
print ("q1+q2 =", q1+q2)
print ("q1*q2 =", q1*q2)
print ("q2*q1 =", q2*q1)
print ("q2.conjugate =", q2.conjugate())

"""
my_vec = Vec3D()
print(my_vec)
rot_vec = g*my_vec
print(rot_vec)
"""


uq1 = GetRotQuaternion(0.5*pi, (1.,0.,0.))
uq2 = GetRotQuaternion(0.25*pi, (0,0.,1.))
print("uq1=", uq1)
print("uq2=", uq2)
grot1 = SO3(uq1)
grot2 = SO3(uq2)

print("grot1 =\n")
print(grot1)
print("grot1.Inverse() =\n")
print(grot1.Inverse())

lievec1 =  Log(grot1)
print("lievec1 =\n")
print(lievec1)

grotexp1 = Exp(lievec1)

print("grotexp1 =\n")
print(grotexp1)


print("grot2 =\n")
print(grot2)

lievec2 =  Log(grot2)
print("lievec2 =\n")
print(lievec2)

grotexp2 = Exp(lievec2)

print("grotexp2 =\n")
print(grotexp2)


uq3= uq1*uq2*(uq1.conjugate())
grot3 = SO3(uq3)
print("grot3 = grot1*grot2*grot1.Inverse()\n")
print(grot3)

lievec3 =  Log(grot3)
print("lievec3 =\n")
print(lievec3)

grotexp3 = Exp(lievec3)

print("grotexp3 =\n")
print(grotexp3)

print("SLERP Interpolation:\n")

uq1 = GetRotQuaternion(0.5*pi, (1.,0.,0.))
uq2 = GetRotQuaternion(0.5*pi, (0.,0.,1.))
grot1 = SO3(uq1)
grot2 = SO3(uq2)
ginter = Interpolate(0.5, grot1, grot2)

print("grot1 =\n")
print(grot1)
print("grot2 =\n")
print(grot2)

print("ginter = \n")
print(ginter)

print(Log(grot1))
print(Log(ginter))
print(Log(grot2))

print(Log(grot1).getvec())
print(Log(ginter).getvec())
print(Log(grot2).getvec())
