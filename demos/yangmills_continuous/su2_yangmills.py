from ngsolve import *

from ngsolve.meshes import *

from ngsolve.solvers import *
from ngsolve.webgui import Draw

import matplotlib.pyplot as plt

sys.path.insert(1, '../quaternion/')

from quaternion_cf import *

#from lgt import *

liedim = 3
g = 1.
#g = 0.


"""
general functions
"""

def generate_scenes(scenes_dict, mesh):
    scenes = []
    for name in scenes_dict.keys():
        print(name)
        scenes.append(Draw (scenes_dict[name][0], mesh, **scenes_dict[name][1]))
    print("\n")
    return scenes


def gaussxy(mu, sigma2):
    #return 1/sqrt(sigma2) * exp(-1/sigma2 *( (x-mu[0])*(x-mu[0]) + (y-mu[1])*(y-mu[1])))
    return gaussx(mu[0], sigma2) * gaussy(mu[1], sigma2)

def gaussx(mu, sigma2):
    return 1./sqrt(2*pi*sigma2) * exp(-0.5/sigma2 *( (x-mu)*(x-mu)))

def gaussy(mu, sigma2):
    return 1./sqrt(2*pi*sigma2) * exp(-0.5/sigma2 *( (y-mu)*(y-mu)))

def gaussxyperx(mu, sigma):
    return gaussxy(mu, sigma)+gaussxy([mu[0]-1,mu[1]], sigma)+gaussxy([mu[0]+1,mu[1]], sigma)

def gaussxypery(mu, sigma):
    return gaussxy(mu, sigma)+gaussxy([mu[0],mu[1]-1], sigma)+gaussxy([mu[0],mu[1]+1], sigma)

def gaussxyper(mu, sigma):
    return gaussxyperx(mu, sigma)+gaussxyperx([mu[0],mu[1]+1], sigma)+gaussxyperx([mu[0],mu[1]-1], sigma)

def square_per_func(func, x0=[0,0], zones=([1,0],[0,1],[1,1],[1,-1])):
    #TODO: truncate function
    
#     per_func_CF = func(x0)
#     for d in range(len(wavel)):
#         per_func_CF = square_perd_func(per_func, d=d, x0=x0 ,wavel=wavel)

    per_func = func(x0)
    for direction in zones:
        x0shifted = [x0[d] + direction[d] for d in range(len(x0))]
        per_func += func(x0shifted)
        x0shifted = [x0[d] - direction[d] for d in range(len(x0))]
        per_func += func(x0shifted)
        
    return per_func

# def square_perd_func(func, d=0, x0=[0,0], wavel=[1,1]):
#     #TODO: truncate function
    
#     perd_func = func(x0)
#     x0shifted = x0
#     x0shifted[d] = x0[d] + wavel[d]
#     perd_func += func(x0shifted)

#     x0shifted = x0
#     x0shifted[d] = x0[d] - wavel[d]
#     perd_func += func(x0shifted)

#     return perd_func



# Using Form_space**3 for su(2) valued forms
# the first index is algebra valued, the second is form valued

# in the following the basis for the Lie Algebra is chosen via {1/(2i) sigma_k}
# where sigma_k are the Pauli matrices (traceless, hermitian)
# this is chosen such that the Lie Bracket becomes a cross product

# (Lie algebra valued) g-forms are represented by CFs of shape [g_dim,space_dim] 
# order of indices for Lie algebra valued CFs therefore is
# FIRST Lie Algebra THEN space
# during operations the order of indices needs to be transposed 

#Extended form Lie Bracket for g-0-form with g-0-form (regular Lie Bracket)
def Bracket_Wedge0f0f(u_gf, v_gf):
    return CF(Cross (u_gf, v_gf))


#Extended form Lie Bracket for g-1-form with g-0-form
def Bracket_Wedge1f0f(u_gf, v_gf):
    #return CF ( (Cross (u_gf[:,0], v_gf), Cross (u_gf[:,1], v_gf)) ).Reshape((2,3)).trans
    return CF ( tuple([ Cross (u_gf[:,d], v_gf) for d in range(u_gf[0,:].dim) ]) ).Reshape((u_gf[0,:].dim,3)).trans
    # TODO: how does reshape work?
    #return CF ( tuple([ Cross (u_gf[:,d], v_gf) for d in range(u_gf[0,:].dim) ]) ).Reshape((3, u_gf[0,:].dim))



#Extended form Lie Bracket for g-1-forms
#wedge product does not divide by overcounting factor
def Bracket_Wedge1f1f(u_gf, v_gf):
    # spatial dim = 2 case
    if ( (u_gf[0,:].dim == 2) and (v_gf[0,:].dim == 2) ):
        return CF ( (Cross(u_gf[:,0], v_gf[:,1]) - Cross(u_gf[:,1], v_gf[:,0])) ).Reshape((3,1))
    
    # spatial dim = 3 case
    elif ( (u_gf[0,:].dim == 3) and (v_gf[0,:].dim == 3) ):
        # TODO: how does reshape work?
        return CF(
            (
            Cross(u_gf[:,1], v_gf[:,2]) - Cross(u_gf[:,2], v_gf[:,1]),
            Cross(u_gf[:,2], v_gf[:,0]) - Cross(u_gf[:,0], v_gf[:,2]),
            Cross(u_gf[:,0], v_gf[:,1]) - Cross(u_gf[:,1], v_gf[:,0])
            ).Reshape((u_gf[0,:].dim,3)).trans
        )
#         return CF(
#             (
#             2*Cross(u_gf[:,1], v_gf[:,2]),
#             2*Cross(u_gf[:,2], v_gf[:,0]),
#             2*Cross(u_gf[:,0], v_gf[:,1])
#             ).Reshape((u_gf[0,:].dim,3)).trans
#         )
    else:
        raise ValueError(f"no implementation for spatial dims {u_gf[0,:].dim} and {v_gf[0,:].dim}")

#TODO: what inner product?
# for Frobenius su2 inner product is the same as R3 for coefficients 
# in the Basis {1/(2i) sigma_k} this is not entirely correct: there should be a (global!) prefactor of 1/8 
def LieInnerProduct(f,g):
    return InnerProduct(f,g)
    #if (f.dims = g.dims):
        #fg = CF(0)
        #for d in range(f[:,0].dim):
        #    fg += InnerProduct(f[d,:], g[d,:])
        #return fg
    #else:
        #RaiseError("dims dont match")
        




"""
Connection specific functions
"""

def grad_A(gauge_A_gf, g_form_gf):
    return grad(g_form_gf) + g*CF(Bracket_Wedge1f0f(gauge_A_gf, g_form_gf))
        
def curl_A(gauge_A_gf, g_form_gf):
    #return curl(g_form_gf) + 0.5*g*CF(Bracket_Wedge1f1f(gauge_A_gf, g_form_gf)) 
    return curl(g_form_gf) + g*CF(Bracket_Wedge1f1f(gauge_A_gf, g_form_gf)) 
    
    
    
    # take not of the 1/2 factor
    # there are two slightly different formulations
    # using forms:
    # d A + 0.5* [A \wedge A]
    # using indices:
    # \partial_\mu A_\nu - \partial_\nu A_\mu + 0.5* ([A_\mu, A_\nu] - [A_\nu, A_\mu])
    # \partial_\mu A_\nu - \partial_\nu A_\mu + [A_\mu, A_\nu]    
def F_A(gauge_A_gf):
    return curl(gauge_A_gf) + 0.5*g*CF(Bracket_Wedge1f1f(gauge_A_gf, gauge_A_gf)) 
    #return curl(gauge_A_gf) + g*CF(Bracket_Wedge1f1f(gauge_A_gf, gauge_A_gf)) 


"""
2D+z effective functions
"""

# Pay attention to the order of index types
def gVector_3DCF(xy_CF, z_CF):
    return CF(
        (xy_CF.trans, z_CF)
    ).Reshape( (3,3) ).trans
   


# TODO: 3D -> 2D effective model, no d_z
# 2Dz functions rely on above grad_A curl_A functions in 2D!


#'overloads' that take 3D CFs
#def grad_A3D_2Dz(gauge_A_gf, g_form_gf): 
#    return grad_A_2Dz(gauge_A_gf[:,0:2], gauge_A_gf[:,2], g_form_gf)
#    
#def curl_A3D_2Dz(gauge_A_gf, g_form_gf):
#    return curl_A_2Dz(gauge_A_gf[:,0:2], gauge_A_gf[:,2], g_form_gf[:,0:2], g_form_gf[:,2])
#
#def F_A3D_2Dz(gauge_A_gf):
#    return F_A_2Dz(gauge_A_gf[:,0:2], gauge_A_gf[:,2])


# TODO: what return type?
# g_form_gf is a function, no z component
def grad_A_2Dz(gauge_A_gf, gauge_Az_gf, g_form_gf): 
    grad_A_g_form_list = []
    
    # gradient components are added after present components in g_form_gf
    grad_A_2D = grad_A(gauge_A_gf, g_form_gf)
    
    # take note of index order of CFs:
    # FIRST Lie Algebra THEN space
    # whereas with this list:
    # FIRST space THEN Lie Algebra
    grad_A_g_form_list.append( grad_A_2D[:,0] )
    grad_A_g_form_list.append( grad_A_2D[:,1] )
    
    # this would be if g_form_gf had a z dependance
    #grad_A_g_form_gf.append(dz(g_form_gf) + g*CF(Bracket_Wedge0f0f(gauge_Az_gf, g_form_gf)) ) 
    # but everything is assumed to be z-symmetric
    # therefore dz(...) = 0:
    grad_A_g_form_list.append( g*Bracket_Wedge0f0f(gauge_Az_gf, g_form_gf) )
    
    #return CF( tuple(grad_A_g_form_list) )
    return CF( tuple(grad_A_g_form_list) ).Reshape((3,3)).trans
        
    
def curl_A_2Dz(gauge_A_gf, gauge_Az_gf, g_form_gf, gz_form_gf):
    
    # take note of index order of CFs:
    # FIRST Lie Algebra THEN space
    # whereas with this list:
    # FIRST space THEN Lie Algebra
    curl_A_g_form_list = []

    # this would be if g_form_gf had z dependance
    #curl_A_g_form_list.append( - dz(g_form_gf[1]) + dy(gz_form_gf)  )
    #curl_A_g_form_list.append( dz(g_form_gf[0]) - dx(gz_form_gf)  )
    # but everything is assumed to be z-symmetric
    # therefore:
    curl_A_g_form_list.append(grad(gz_form_gf)[:,1])
    curl_A_g_form_list.append(-grad(gz_form_gf)[:,0])

    #no 0.5 factor (present in curvature F_A!)
    # careful with the signs: they come from the epsilon symbol in the wedge product 
    curl_A_g_form_list[0] += g*(
        (-1.)*Bracket_Wedge0f0f(gauge_Az_gf, g_form_gf[:,1]) + Bracket_Wedge0f0f(gauge_A_gf[:,1], gz_form_gf) 
    )
    curl_A_g_form_list[1] += g*(
        Bracket_Wedge0f0f(gauge_Az_gf, g_form_gf[:,0]) - Bracket_Wedge0f0f(gauge_A_gf[:,0], gz_form_gf)
    )

    #z component corresponds to 2D curl_A
    curl_A_g_form_list.append( curl_A(gauge_A_gf, g_form_gf))
    return CF( tuple(curl_A_g_form_list) ).Reshape((3,3)).trans


# this function is virtually identical to curl_A_2Dz
# except for the 1/2 factor present in curvature
#TODO: 1/2 curl factor?
def F_A_2Dz(gauge_A_gf, gauge_Az_gf):
 
    # take note of index order of CFs:
    # FIRST Lie Algebra THEN space
    # whereas with this list:
    # FIRST space THEN Lie Algebra
    curl_A_g_form_list = []

    # this would be if g_form_gf had z dependance
    #curl_A_g_form_list.append( - dz(g_form_gf[1]) + dy(gz_form_gf)  )
    #curl_A_g_form_list.append( dz(g_form_gf[0]) - dx(gz_form_gf)  )
    # but everything is assumed to be z-symmetric
    # therefore:
    curl_A_g_form_list.append(grad(gauge_Az_gf)[:,1])
    curl_A_g_form_list.append(-grad(gauge_Az_gf)[:,0])

    #0.5 factor present in curvature F_A
    # as both summands in one term are identical, it could be dropped in principle...
    # careful with the signs: they come from the epsilon symbol in the wedge product 
    curl_A_g_form_list[0] += 0.5*g*(
     (-1.)*Bracket_Wedge0f0f(gauge_Az_gf, gauge_A_gf[:,1]) + Bracket_Wedge0f0f(gauge_A_gf[:,1], gauge_Az_gf) 
    )
    curl_A_g_form_list[1] += 0.5*g*(
     Bracket_Wedge0f0f(gauge_Az_gf, gauge_A_gf[:,0]) - Bracket_Wedge0f0f(gauge_A_gf[:,0], gauge_Az_gf)
    )


    #z component ALMOST corresponds to 2D curl_A
    # there is no 0.5 present in curl_A!
    # this is wrong:
    #curl_A_g_form_list.append( curl_A(gauge_A_gf, gauge_A_gf)) 
    # this is correct:
    curl_A_g_form_list.append( F_A(gauge_A_gf) )
    #curl_A_g_form_list.append( curl(gauge_A_gf) + 0.5*Bracket_Wedge1f1f(gauge_A_gf, gauge_A_gf))
    
    
    
    return CF( tuple(curl_A_g_form_list) ).Reshape((3,3)).trans

    #return curl_A_2Dz(gauge_A_gf, gauge_Az_gf, gauge_A_gf, gauge_Az_gf)



   

"""
Functions for gauging CFs
"""

# gauge a 0form using a quaternion CF
def gauge_0f(g_0f_CF, phi_CF, n_CF):
    # 1: choose Lie group element to gauge with
    # 2: convert Lie group element into quaternion space using an isomorphism
    # (we directly parametrize gauge group element using 'angle' and 'axis')
    qgauge_CF = get_rotqCF(phi_CF,n_CF)
    #print(qgauge_CF[0, ...].dim)
    # 2: convert Lie algebra element E into quaternion space using a isomorphism
    # this ismorphism (1/i) sigma_k --> k is chosen
    # look above to see what basis for the Lie Algebra is chosen ...
    # it should be (1/(2i)) sigma_k
    # therefore the transformation of components from basis (1/(2i)) sigma_k to (1/(i)) sigma_k, is:
    qg_0f_CF = q_from_vecCF(0.5*g_0f_CF)
    # 3: apply gauge by q E q.conj with quaternion multiplication
    qg_0f_CF = qmulCF( qmulCF(qgauge_CF, qg_0f_CF), qconjCF(qgauge_CF) )
    # 4: convert E back to Lie algebra using opposite component transformation
    g_gauged_0f_CF = vec_from_q(2*qg_0f_CF)
    
    return g_gauged_0f_CF
    
# gauge a 1form using a quaternion CF
def gauge_1f(g_1f_CF, phi_CF, n_CF):
    gauged_g_1f_list = []
    for l in range(g_1f_CF.dims[1]):
        gauged_g_1f_list.append( gauge_0f(g_1f_CF[:,l], phi_CF, n_CF) )
    
    return CF ( tuple(gauged_g_1f_list) ).Reshape( (g_1f_CF.dims[1],3) ).trans

# gauge gauge field A using a quaternion CF
def gauge_A(A_CF, phi_CF, n_CF):
    # 1: choose Lie group element to gauge with
    
    if A_CF.dim > 3:
        gAginv_CF = gauge_1f(A_CF, phi_CF, n_CF)
    elif A_CF.dim == 3:
        gAginv_CF = gauge_0f(A_CF, phi_CF, n_CF)
    
    # dont forget to gauge A with -dg*g_inv!
    # special formula for SU(2):
    #dgginv = grad(phi_CF)*n_CF + 0.5 * sin(2*phi_CF)*grad(n_CF)
    
    gradphi_CF = CF( (phi_CF.Diff(x), phi_CF.Diff(y)) )
    gradn_CF = CF( (n_CF.Diff(x), n_CF.Diff(y)) )
    
    # index associated with n = Lie algebra index
    # index associated with grad = spatial index
    # order of indices for Lie algebra valued CFs is
    # FIRST Lie Algebra THEN space
    dgginv_CF = OuterProduct(n_CF, gradphi_CF) + sin(phi_CF)*gradn_CF
    
    # take care of the sign convention for gauges
    #return gAginv_CF + dgginv_CF
    return gAginv_CF - dgginv_CF

"""
Evaluating observables
"""

# energy is just square value of fields
# charge is solved for in a secondary problem
