#from ngsolve import CF, Trace, Conj, cos, sin, InnerProduct, Cross, Norm
from ngsolve import *

#from ngsolve.meshes import *
#from ngsolve.solvers import *
#from ngsolve.webgui import Draw
#from netgen.geom2d import unit_square


import numpy as np







def paulimat(i):
    if (i == 0):
        return CF( (1,0,0,1) ).Reshape( (2,2) )
    elif (i == 1):
        return CF( (0,1,1,0) ).Reshape( (2,2) )
    elif (i == 2):
        return CF( (0,-1j,1j,0) ).Reshape( (2,2) )
    elif (i == 3):
        return CF( (1,0,0,-1) ).Reshape( (2,2) )
    else:
        raise ValueError(f"there is no {i}-th pauli matrix")

# quaternion-su(2) isomorphism based on
# 1/i * sigma_k = -i sigma_k = k
def su2mat_from_q(q):
    su2mat = (q[0]*CF( (1, 0, 0, 1) )).Reshape( (2,2) )
    for i in range(1,4):
        su2mat = su2mat - 1j*q[i]*paulimat(i)
        
    return su2mat.Reshape( (2,2) )
    

    
def matmul(A,B):
    matlist = []
    for pos in [(0,0), (0,1), (1,0), (1,1)]:
        #print(i,j)
        i = pos[0]
        j = pos[1]
        #matlist.append(InnerProduct(A[i,:], B[:,j]))
        matlist.append(A[i,0]*B[0,j]+A[i,1]*B[1,j])
        
        
    return CF( tuple(matlist) ).Reshape( (2,2) )
   
def tr_ah(U):
    return 0.5*CF(
        tuple(
            [
                (Trace(paulimat(k).Reshape((2,2))*U.Reshape((2,2)))).imag for k in range(1,4)
            ]
        )
    )





"""Quaternion functions"""

#def su2_from_q(q):
#    return 2*q

#def q_from_SU2(q):
#    return 0.5*q

def q_from_vecCF(v):
    return CF(
            (
                0,
                CF(tuple(v))
            )
    )

def vec_from_q(q):
    return q[1:]

def axis(theta, phi):
    return np.array(
        [
            sin(theta)*sin(phi),
            sin(theta)*cos(phi),
            cos(theta)
        ]
    )

def get_rotqCF(phi, n):
    return CF( 
        tuple([cos(0.5*phi), sin(0.5*phi)*CF(tuple(n))])
    )


def get_rotq(phi, n):
    return np.array([cos(0.5*phi), *[sin(0.5*phi)*n_i for n_i in n]])


#def qexpCF(q):
#    
#    vqnorm = Norm(q[1:])
#    expq0 = Exp(q[0])
#    
#    if vqnorm < 10**(-10):
#        return CF( (expq0,0,0,0) )
#    else:
#        return expq0*CF( 
#            tuple([Cos(vqnorm), (Sin(vqnorm)/vqnorm)*q[1:]])
#        )

def qexp(q):
    
    vqnorm = np.linalg.norm(q[1:])
    #expq0 = Exp(q[0])
    expq0 = np.exp(q[0])
    return expq0*np.array(
            [ cos(vqnorm), *((sin(vqnorm)/vqnorm)*q[1:]) ], dtype=type(vqnorm)
        )


def qexpCF(q):
    
    vqnorm = Norm(q[1:])
    #expq0 = Exp(q[0])
    expq0 = exp(q[0])
    return expq0*CF( 
            tuple([cos(vqnorm), (sin(vqnorm)/vqnorm)*q[1:]])
        )

def qlogCF(q):
    
    qnorm = Norm(q)
    vqnorm = Norm(q[1:])
    
    #if (qnorm < 10**(-10)) or (vqnorm < 10**(-10)):
    #    raise ValueError(f"unable to take ln of quaternion with vanishing norm")
    #else:
    #    return CF( 
    #        tuple([Log(qnorm), Arccos(q[0]/qnorm)/vqnorm*q[1:]])
    #    )
    
    #return CF( 
    #        tuple([log(qnorm), acos(q[0]/qnorm)/vqnorm*q[1:]])
    #    )

    return CF( 
            tuple([log(qnorm), acos(vqnorm/qnorm)/vqnorm*q[1:]])
        )

def qmul(q1,q2):
    
    return np.array(
            [
                q1[0] * q2[0] - np.inner(q1[1:], q2[1:]),
                *(
                q1[0] * np.array(q2[1:]) + q2[0] * np.array(q1[1:]) + np.cross(q1[1:], q2[1:]) 
                )
            ]
    )

def qmulCF(q1,q2):
    
    return CF(
        (
            q1[0] * q2[0] - InnerProduct(q1[1:], q2[1:]),
            CF( q1[0] * q2[1:] + q2[0] * q1[1:] + Cross(q1[1:], q2[1:]) )
        )
    )




def qconj(q):
    return np.array([q[0], *(-q[1:])])


def qconjCF(q):
    return CF( 
        (
            q[0],
            -CF(q[1:])
        )
    )