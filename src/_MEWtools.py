#from __future__ import division
import numpy as np
import sympy as sym
#import time


r, θ, φ, t = sym.symbols(  'r θ φ t'    ,   real=True  )




############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
#                                 Love Number / y-Function Calculator                                      #
def lovefunc(L_, rin_, pin_, uin_, nin_, w_):                                                              #
    
    G = 6.67e-11
    
    
    # Number of layers in the body
    layno = len(rin_)
    
    
    # Find Purely fluid or elastic layers and create the index lists
    fluidlayerlist   = np.where(np.array(uin_) == 0)[0]
    elasticlayerlist = np.where(np.array(nin_) == 0)[0]
    
    
    # If there are no fluid layers, then create this list instead
    if len(fluidlayerlist) == 0:
        fluidlayerlist = [None]
    
    
    # Turn rigidity into viscoelastic parameter for those layers
    uV=[]
    for i in range(layno):
        if i in elasticlayerlist:
            uV.append(uin_[i])
        else:
            uV.append(uin_[i] / (1-1j*(uin_[i]/(nin_[i]*w_)) ))
    
    
    # Make a new lists of inputs
    rlist=rin_
    plist=pin_
    ulist=uV
    
    
    # Gravitational Acceleration Function
    def gfunc(r_, lay_):
        
        rsw0 = np.insert(rlist,0,0)
        
        gcalc = plist[lay_] * (r_**3 - rsw0[lay_]**3)
        
        for i in range(0,lay_):
            gcalc += plist[i] * (rsw0[i+1]**3 - rsw0[i]**3)
        
        return ((4*np.pi*G)/(3*r_**2)) * gcalc
    
    
    # Fundamental Matrix for Elastic Solutions
    def Emat(r_, g_, rho_, mu_):
        return sym.Matrix([
                           [                                             r_**(1+L_) ,                                r_**(-1+L_) , 0                      ,                                               r_**(-L_)   ,                                r_**(-2-L_) , 0                  ],
                           [ ((3+L_)/(L_*(1+L_)))                      * r_**(1+L_) , (1/L_)                       * r_**(-1+L_) , 0                      , -((-2+L_)/(L_*(1+L_)))                      * r_**(-L_)   , -(1/(1+L_))                  * r_**(-2-L_) , 0                  ],
                           [ (2*mu_*((-3+L_*(-1+L_))/L_) + r_*rho_*g_) * r_**(L_)   , (2*mu_*(-1+L_) + r_*rho_*g_) * r_**(-2+L_) , rho_     * r_**(L_)    , (2*mu_*((1-L_*(3+L_))/(1+L_)) + r_*rho_*g_) * r_**(-1-L_) , (2*mu_*(-2-L_) + r_*rho_*g_) * r_**(-3-L_) , rho_ * r_**(-1-L_) ],
                           [ 2*mu_*((2+L_)/(1+L_))                     * r_**(L_)   , 2*mu_*((-1+L_)/(L_))         * r_**(-2+L_) , 0                      , 2*mu_*((-1+L_)/(L_))                        * r_**(-1-L_) , 2*mu_*((2+L_)/(1+L_))        * r_**(-3-L_) , 0                  ],
                           [ 0                                                      , 0                                          ,            r_**(L_)    , 0                                                         , 0                                          ,        r_**(-1-L_) ],
                           [ 4*np.pi*G*rho_                           * r_**(1+L_) , 4*np.pi*G*rho_              * r_**(-1+L_) , (1+2*L_) * r_**(-1+L_) , 4*np.pi*G*rho_                             * r_**(-L_)   , 4*np.pi*G*rho_              * r_**(-2-L_) , 0                  ]
                           ])
    
    
    # Fundamental Matrix for Elastic Solutions Invserse
    def EmatINV(r_, g_, rho_, mu_):
        return sym.Matrix([
                           [ -(((2*mu_*(2+L_) - r_*rho_*g_)*L_*(1+L_))/(2*mu_*(3+4*L_*(2+L_))))         * r_**(-1-L_) , ((L_**2 *(1+L_)*(2+L_)) / (3+4*L_*(2+L_)))  * r_**(-1-L_) , -((L_*(1+L_))/(2*mu_*(3+4*L_*(2+L_)))) * r_**(-L_)  , ((L_**2 *(1+L_))/(2*mu_*(3+4*L_*(2+L_))))   * r_**(-L_)  , ((L_*(1+L_)*rho_)/(2*mu_*(3+4*L_*(2+L_))))  * r_**(-L_)  , 0                          ],
                           [  ((L_*(2*mu_*(L_*(3+L_)-1)-(1+L_)*r_*rho_*g_))/(2*mu_*(-1+4*L_**2)))       * r_**(1-L_)  , -((L_*(-1+L_)*(1+L_)**2)/(-1+4*L_**2))       * r_**(1-L_)  ,  ((L_*(1+L_))/(2*mu_*(-1+4*L_**2)))    * r_**(2-L_) , -(((-2+L_)*L_*(1+L_))/(2*mu_*(-1+4*L_**2))) * r_**(2-L_) , -((L_*(1+L_)*rho_)/(2*mu_*(-1+4*L_**2)))    * r_**(2-L_) , 0                          ],
                           [ -((4*np.pi*G*rho_)/(1+2*L_))                                              * r_**(1-L_)  , 0                                                         , 0                                                   , 0                                                        , 0                                                        , (1/(1+2*L_)) * r_**(1-L_)  ],
                           [  ((L_*(1+L_)*(2*mu_*(-1+L_) + r_*rho_*g_))/(2*mu_*(-1+4*L_**2)))           * r_**(L_)    , ((L_*(-1+L_)*(1+L_)**2)/(-1+4*L_**2))        * r_**(L_)    , -((L_*(1+L_))/(2*mu_*(-1+4*L_**2)))    * r_**(1+L_) , -((L_*(1+L_)**2)/(2*mu_*(-1+4*L_**2)))      * r_**(1+L_) , ((L_*(1+L_)*rho_)/(2*mu_*(-1+4*L_**2)))     * r_**(1+L_) , 0                          ],
                           [ -(((1+L_)*(2*mu_*(-3+L_*(-1+L_))+ L_*r_*rho_*g_))/(2*mu_*(3+4*L_*(2+L_)))) * r_**(2+L_)  , -((L_**2 *(1+L_)*(2+L_))/((3+4*L_*(2+L_)))) * r_**(2+L_)  , ((L_*(1+L_))/(2*mu_*(3+4*L_*(2+L_))))  * r_**(3+L_) , ((L_*(1+L_)*(3+L_))/(2*mu_*(3+4*L_*(2+L_)))) * r_**(3+L_) , -((L_*(1+L_)*rho_)/(2*mu_*(3+4*L_*(2+L_)))) * r_**(3+L_) , 0                          ],
                           [  ((4*np.pi*G*rho_)/(1+2*L_))                                              * r_**(2+L_)  , 0                                                         , 0                                                   , 0                                                        ,                                               r_**(1+L_) , -(1/(1+2*L_)) * r_**(2+L_) ]
                           ])
    
    
    # Fundamental Matrix for Fluid Solutions
    def Fmat(r_, g_, rho_):
        return sym.Matrix([
                           [ -(1/g_)                                              * r_**L_      , -(1/g_)                                               * r_**(-1-L_) ],
                           [ ((4*np.pi*G*rho_*r_ - (4+L_)*g_)/(L_*(1+L_)*g_**2)) * r_**L_      , ((4*np.pi*G*rho_*r_ + (-3+L_)*g_)/(L_*(1+L_)*g_**2)) * r_**(-1-L_) ],
                           [ 0                                                                  , 0                                                                   ],
                           [ 0                                                                  , 0                                                                   ],
                           [                                                        r_**L_      ,                                                         r_**(-1-L_) ],
                           [ ((-4*np.pi*G*rho_*r_ + (1+2*L_)*g_)/(g_))           * r_**(-1+L_) , -((4*np.pi*G*rho_)/(g_))                             * r_**(-1-L_) ]
                           ])
    
    
    # Fundamental Matrix for Fluid Solutions Inverse
    def FmatINV(r_, g_, rho_):
        return sym.Matrix([
                           [ ((4*np.pi*G*rho_)/((1+2*L_)*g_))                  * r_**(1-L_) ,  (1/(1+2*L_)) * r_**(1-L_) ],
                           [ (((1+2*L_)*g_ - 4*np.pi*G*rho_*r_)/((1+2*L_)*g_)) * r_**(1+L_) , -(1/(1+2*L_)) * r_**(2+L_) ]
                           ])
    
    
    # Make empty lists for boundary conditions resulting from elastic->fluid interface
    fluidBC     = []
    fluidBCvals = sym.Matrix(0, 1, [])
    
    
    # Make empty lists for each layer's fundimental matrix
    mats=[None]*(layno)
    
    
    # Determine if the core is fluid or not & produce core matrix
    if fluidlayerlist[0]==0:
        mats[0] = Fmat(r, gfunc(r, 0), plist[0])[:,0].applyfunc(sym.expand)
    else:
        mats[0] = Emat(r, gfunc(r, 0), plist[0], ulist[0])[:,0:3].applyfunc(sym.expand)
    
    
    # Produce the matrix for all the other layers
    for i in range(1,layno):
        if i-1 in fluidlayerlist and i in fluidlayerlist:
            
            mats[i] = (Fmat(r, gfunc(r, i), plist[i]) * FmatINV(rlist[i-1], gfunc(rlist[i-1], i-1), plist[i]) * mats[i-1][4:,:].subs(r,rlist[i-1])).applyfunc(sym.expand)
        #            print "fluid fluid"
        
        
        elif i-1 in fluidlayerlist and i not in fluidlayerlist:
            
            mattemp = mats[i-1][:,:].subs(r,rlist[i-1])
            mattemp[1:2,:]=sym.zeros(1,len(mattemp[1:2,:]))
            mattemp2 = sym.Matrix([[-1,0],[0,1],[-plist[i-1]*gfunc(rlist[i-1],i-1),0],[0,0],[0,0],[-4 *np.pi *G *plist[i-1],0]]).col_insert(0,mattemp)
            
            
            mats[i] = (Emat(r, gfunc(r, i), plist[i], ulist[i]) * EmatINV(rlist[i-1], gfunc(rlist[i-1], i-1), plist[i], ulist[i]) * mattemp2).applyfunc(sym.expand)
        #            print "fluid elastic"
        
        
        elif i-1 not in fluidlayerlist and i in fluidlayerlist:
            
            fluidBC.append((mats[i-1][2:3,:] - plist[i] * gfunc(rlist[i-1], i-1) * ((mats[i-1][4:5,:] / gfunc(rlist[i-1], i-1)) + mats[i-1][0:1,:])).subs(r,rlist[i-1]))
            fluidBC.append(mats[i-1][3:4,:].subs(r,rlist[i-1]))
            
            fluidBCvals = sym.Matrix([[0]]).row_insert(0,fluidBCvals)
            fluidBCvals = sym.Matrix([[0]]).row_insert(0,fluidBCvals)
            
            mats[i] = (Fmat(r, gfunc(r, i), plist[i]) * FmatINV(rlist[i-1], gfunc(rlist[i-1], i-1), plist[i]) * (mats[i-1][5:6,:] - 4*np.pi*G*plist[i]*((mats[i-1][4:5,:] / gfunc(rlist[i-1], i-1)) + mats[i-1][0:1,:])).subs(r,rlist[i-1]).row_insert(0,mats[i-1][4:5,:].subs(r,rlist[i-1]))).applyfunc(sym.expand)
        #            print "elastic fluid"
        
        
        elif i-1 not in fluidlayerlist and i not in fluidlayerlist:
            
            mats[i] = (Emat(r, gfunc(r, i), plist[i], ulist[i]) * EmatINV(rlist[i-1], gfunc(rlist[i-1], i-1), plist[i], ulist[i]) * mats[i-1].subs(r,rlist[i-1])).applyfunc(sym.expand)
    #            print "elastic elastic"
    
    
    
    # Give that matrix right at the body surface
    surfmat = mats[-1].subs(r,rlist[-1]).applyfunc(sym.expand)
    #    print "surface matrix"
    
    
    # Pad the list of elastic->fluid interface boundary condition with zeros to the right (so all rows have the same ammt of columns for BC constants)
    for i in range(len(fluidBC)):
        fluidBC[i] = sym.Matrix(np.lib.pad(fluidBC[i],((0,0),(0,len(surfmat[2:3,:])-len(fluidBC[i]))),'constant'))
    
    # Recast the e->f BC list as a matrix
    fluidBC = sym.Matrix(fluidBC)
    
    
    # Join the surface BCs to the e->f BCs
    if fluidlayerlist[-1]==layno-1:
        BCmat  = fluidBC.col_join(surfmat[5:6,:])
        BCvals = fluidBCvals.col_join(sym.Matrix([[-(2*L_+1)/rlist[-1]]]))
    elif fluidBC == sym.Matrix(0, 0, []):
        BCmat  = sym.Matrix([surfmat[2:3,:],surfmat[3:4,:],surfmat[5:6,:]])
        BCvals = sym.Matrix([[0],[0],[-(2*L_+1)/rlist[-1]]])
    else:
        BCmat  = fluidBC.col_join( sym.Matrix([surfmat[2:3,:],surfmat[3:4,:],surfmat[5:6,:]]))
        BCvals = fluidBCvals.col_join(sym.Matrix([[0],[0],[-(2*L_+1)/rlist[-1]]]))
    #    print "BC stuff"
    


    # Solve for the constants of integration
    consts = (BCmat.inv() * BCvals).applyfunc(sym.simplify)
    #    print "constants"
    
    

    
    # Calculate the Love numbers
    lovenums = [sym.expand(gfunc(rlist[-1], layno-1)*(surfmat[0:1,:]*consts)[0]) ,
                sym.expand(gfunc(rlist[-1], layno-1)*(surfmat[1:2,:]*consts)[0]) ,
                sym.expand(-((surfmat[4:5,:]*consts)[0] + 1)) ]
#    print "love numbers"
                
                
                
    # Calculate the y-functions in every layer
    yfuncs = []
    for i in range(len(mats)):
        yfuncs.append(sum((mats[i]*sym.Matrix(consts[:len(mats[i][0:1,:])])).applyfunc(sym.expand).tolist(),[]))
#    print "y functions"
    
    

    
    


    return lovenums


############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
#                                 Tide Calculator    

#n, S, N, Ln, opn, t = sym.symbols('n S N Ln opn t',real=True)

def freqsinorder(n_, S_, N_, Ln_, opn_):
    return np.abs([ 2*n_ - 2*S_ - 2*N_ + 0*Ln_ + 0*opn_,
             1*n_ - 2*S_ - 2*N_ + 0*Ln_ + 0*opn_,
             3*n_ - 2*S_ - 2*N_ + 0*Ln_ + 0*opn_,
             0*n_ + 2*S_ + 2*N_ + 0*Ln_ - 2*opn_,
             1*n_ - 2*S_ - 2*N_ + 0*Ln_ + 2*opn_,
             1*n_ + 2*S_ + 2*N_ + 0*Ln_ - 2*opn_,
             2*n_ + 2*S_ + 2*N_ + 0*Ln_ - 4*opn_,
             3*n_ + 2*S_ + 2*N_ + 0*Ln_ - 4*opn_,
             1*n_ + 2*S_ + 2*N_ + 0*Ln_ - 4*opn_,
            -2*n_ + 2*S_ + 2*N_ + 1*Ln_ - 0*opn_,
             2*n_ - 2*S_ - 2*N_ + 1*Ln_ - 0*opn_,
            -2*n_ - 2*S_ - 2*N_ + 1*Ln_ + 4*opn_,
             0*n_ + 2*S_ + 2*N_ + 1*Ln_ - 2*opn_,
             0*n_ - 2*S_ - 2*N_ + 1*Ln_ + 2*opn_,
             2*n_ + 2*S_ + 2*N_ + 1*Ln_ - 4*opn_,
             2*n_ - 1*S_ - 1*N_ + 0*Ln_ - 1*opn_,
             3*n_ - 1*S_ - 1*N_ + 0*Ln_ - 1*opn_,
             1*n_ - 1*S_ - 1*N_ + 0*Ln_ - 1*opn_,
             0*n_ + 1*S_ + 1*N_ + 0*Ln_ - 1*opn_,
             1*n_ + 1*S_ + 1*N_ + 0*Ln_ - 1*opn_,
             1*n_ - 1*S_ - 1*N_ + 0*Ln_ + 1*opn_,
             2*n_ + 1*S_ + 1*N_ + 0*Ln_ - 3*opn_,
             3*n_ + 1*S_ + 1*N_ + 0*Ln_ - 3*opn_,
             1*n_ + 1*S_ + 1*N_ + 0*Ln_ - 3*opn_,
             2*n_ + 1*S_ + 1*N_ + 1*Ln_ - 3*opn_,
             0*n_ + 1*S_ + 1*N_ + 1*Ln_ - 1*opn_,
             2*n_ - 1*S_ - 1*N_ + 1*Ln_ - 1*opn_,
            -2*n_ + 1*S_ + 1*N_ + 1*Ln_ + 1*opn_,
             0*n_ - 1*S_ - 1*N_ + 1*Ln_ + 1*opn_,
            -2*n_ - 1*S_ - 1*N_ + 1*Ln_ + 3*opn_,
             1*n_ + 0*S_ + 0*N_ + 0*Ln_ + 0*opn_,
             2*n_ + 0*S_ + 0*N_ + 0*Ln_ - 2*opn_,
             1*n_ + 0*S_ + 0*N_ + 0*Ln_ - 2*opn_,
             3*n_ + 0*S_ + 0*N_ + 0*Ln_ - 2*opn_,
             0*n_ + 0*S_ + 0*N_ + 0*Ln_ + 0*opn_])


def potsinorder(G_, M_, R_, a_, e_, o_, op_, LA_, Lp_, n_, S_, N_, Ln_, opn_, t_ ,theta_, phi_):
    return [  ((3*G_*M_*R_**2)/(4*a_**3))   *              sym.cos(o_/2)**4   *                sym.sin(theta_)**2                     *    sym.cos(freqsinorder(n_, S_, N_, Ln_, opn_)[0]*t_ - 2*phi_),
             -((3*G_*M_*R_**2)/(8*a_**3))   * e_ *         sym.cos(o_/2)**4   *                sym.sin(theta_)**2                     *    sym.cos(freqsinorder(n_, S_, N_, Ln_, opn_)[1]*t_ - 2*phi_),
             ((21*G_*M_*R_**2)/(8*a_**3))   * e_ *         sym.cos(o_/2)**4   *                sym.sin(theta_)**2                     *    sym.cos(freqsinorder(n_, S_, N_, Ln_, opn_)[2]*t_ - 2*phi_),
              ((3*G_*M_*R_**2)/(8*a_**3))   *              sym.sin(o_)**2     *                sym.sin(theta_)**2                     *    sym.cos(freqsinorder(n_, S_, N_, Ln_, opn_)[3]*t_ + 2*phi_ - 2*op_),
              ((9*G_*M_*R_**2)/(16*a_**3))  * e_ *         sym.sin(o_)**2     *                sym.sin(theta_)**2                     *    sym.cos(freqsinorder(n_, S_, N_, Ln_, opn_)[4]*t_ - 2*phi_ + 2*op_),
              ((9*G_*M_*R_**2)/(16*a_**3))  * e_ *         sym.sin(o_)**2     *                sym.sin(theta_)**2                     *    sym.cos(freqsinorder(n_, S_, N_, Ln_, opn_)[5]*t_ + 2*phi_ - 2*op_),
              ((3*G_*M_*R_**2)/(4*a_**3))   *              sym.sin(o_/2)**4   *                sym.sin(theta_)**2                     *    sym.cos(freqsinorder(n_, S_, N_, Ln_, opn_)[6]*t_ + 2*phi_ - 4*op_),
             ((21*G_*M_*R_**2)/(8*a_**3))   * e_ *         sym.sin(o_/2)**4   *                sym.sin(theta_)**2                     *    sym.cos(freqsinorder(n_, S_, N_, Ln_, opn_)[7]*t_ + 2*phi_ - 4*op_),
             -((3*G_*M_*R_**2)/(8*a_**3))   * e_ *         sym.sin(o_/2)**4   *                sym.sin(theta_)**2                     *    sym.cos(freqsinorder(n_, S_, N_, Ln_, opn_)[8]*t_ + 2*phi_ - 4*op_),
              ((3*G_*M_*R_**2)/(4*a_**3))   * e_ * LA_ *   sym.cos(o_/2)**4   *                sym.sin(theta_)**2                     *    sym.cos(freqsinorder(n_, S_, N_, Ln_, opn_)[9]*t_ + 2*phi_ + 0*op_ + Lp_),
             -((3*G_*M_*R_**2)/(4*a_**3))   * e_ * LA_ *   sym.cos(o_/2)**4                *   sym.sin(theta_)**2                     *    sym.cos(freqsinorder(n_, S_, N_, Ln_, opn_)[10]*t_ - 2*phi_ + 0*op_ + Lp_),
             -((3*G_*M_*R_**2)/(4*a_**3))   * e_ * LA_ *   sym.sin(o_/2)**4                *   sym.sin(theta_)**2                     *    sym.cos(freqsinorder(n_, S_, N_, Ln_, opn_)[11]*t_ - 2*phi_ + 4*op_ + Lp_),
              ((3*G_*M_*R_**2)/(8*a_**3))   * e_ * LA_ *   sym.sin(o_)**2                  *   sym.sin(theta_)**2                     *    sym.cos(freqsinorder(n_, S_, N_, Ln_, opn_)[12]*t_ + 2*phi_ - 2*op_ + Lp_),
             -((3*G_*M_*R_**2)/(8*a_**3))   * e_ * LA_ *   sym.sin(o_)**2                  *   sym.sin(theta_)**2                     *    sym.cos(freqsinorder(n_, S_, N_, Ln_, opn_)[13]*t_ - 2*phi_ + 2*op_ + Lp_),
              ((3*G_*M_*R_**2)/(4*a_**3))   * e_ * LA_ *   sym.sin(o_/2)**4                *   sym.sin(theta_)**2                     *    sym.cos(freqsinorder(n_, S_, N_, Ln_, opn_)[14]*t_ + 2*phi_ - 4*op_ + Lp_),
             -((3*G_*M_*R_**2)/(16*a_**3))  *              (2*sym.sin(o_) + sym.sin(2*o_)) *   sym.sin(2*theta_)                      *    sym.sin(freqsinorder(n_, S_, N_, Ln_, opn_)[15]*t_ - 1*phi_ - 1*op_),
             -((21*G_*M_*R_**2)/(32*a_**3)) * e_ *         (2*sym.sin(o_) + sym.sin(2*o_)) *   sym.sin(2*theta_)                      *    sym.sin(freqsinorder(n_, S_, N_, Ln_, opn_)[16]*t_ - 1*phi_ - 1*op_),
              ((3*G_*M_*R_**2)/(32*a_**3))  * e_ *         (2*sym.sin(o_) + sym.sin(2*o_)) *   sym.sin(2*theta_)                      *    sym.sin(freqsinorder(n_, S_, N_, Ln_, opn_)[17]*t_ - 1*phi_ - 1*op_),
             -((3*G_*M_*R_**2)/(8*a_**3))   *              sym.sin(2*o_)                   *   sym.sin(2*theta_)                      *    sym.sin(freqsinorder(n_, S_, N_, Ln_, opn_)[18]*t_ + 1*phi_ - 1*op_),
             -((9*G_*M_*R_**2)/(16*a_**3))  * e_ *         sym.sin(2*o_)                   *   sym.sin(2*theta_)                      *    sym.sin(freqsinorder(n_, S_, N_, Ln_, opn_)[19]*t_ + 1*phi_ - 1*op_),
              ((9*G_*M_*R_**2)/(16*a_**3))  * e_ *         sym.sin(2*o_)                   *   sym.sin(2*theta_)                      *    sym.sin(freqsinorder(n_, S_, N_, Ln_, opn_)[20]*t_ - 1*phi_ + 1*op_),
             -((3*G_*M_*R_**2)/(2*a_**3))   *              sym.sin(o_) * sym.sin(o_/2)**2  *   sym.sin(theta_) * sym.cos(theta_)      *    sym.sin(freqsinorder(n_, S_, N_, Ln_, opn_)[21]*t_ + 1*phi_ - 3*op_),
             -((21*G_*M_*R_**2)/(4*a_**3))  * e_ *         sym.sin(o_) * sym.sin(o_/2)**2  *   sym.sin(theta_) * sym.cos(theta_)      *    sym.sin(freqsinorder(n_, S_, N_, Ln_, opn_)[22]*t_ + 1*phi_ - 3*op_),
             ((3*G_*M_*R_**2)/(4*a_**3))    * e_ *         sym.sin(o_) * sym.sin(o_/2)**2  *   sym.sin(theta_) * sym.cos(theta_)      *    sym.sin(freqsinorder(n_, S_, N_, Ln_, opn_)[23]*t_ + 1*phi_ - 3*op_),
             -((3*G_*M_*R_**2)/(4*a_**3))   * e_ * LA_ *   sym.sin(o_) * sym.sin(o_/2)**2  *   sym.sin(theta_) * sym.cos(theta_)      *    sym.sin(freqsinorder(n_, S_, N_, Ln_, opn_)[24]*t_ + 1*phi_ - 3*op_ + Lp_),
             -((3*G_*M_*R_**2)/(16*a_**3))   * e_ * LA_ *  sym.sin(2*o_)                  *   sym.sin(2*theta_)                       *    sym.sin(freqsinorder(n_, S_, N_, Ln_, opn_)[25]*t_ + 1*phi_ - 1*op_ + Lp_),
              ((3*G_*M_*R_**2)/(32*a_**3))   * e_ * LA_ *  (2*sym.sin(o_) + sym.sin(2*o_)) *   sym.sin(2*theta_)                      *    sym.sin(freqsinorder(n_, S_, N_, Ln_, opn_)[26]*t_ - 1*phi_ - 1*op_ + Lp_),
              ((3*G_*M_*R_**2)/(32*a_**3))   * e_ * LA_ *  (2*sym.sin(o_) + sym.sin(2*o_)) *   sym.sin(2*theta_)                      *    sym.sin(freqsinorder(n_, S_, N_, Ln_, opn_)[27]*t_ + 1*phi_ + 1*op_ + Lp_),
             -((3*G_*M_*R_**2)/(16*a_**3))   * e_ * LA_ *  sym.sin(2*o_)                   *   sym.sin(2*theta_)                      *    sym.sin(freqsinorder(n_, S_, N_, Ln_, opn_)[28]*t_ - 1*phi_ + 1*op_ + Lp_),
             -((3*G_*M_*R_**2)/(8*a_**3))   * e_ * LA_ *  sym.sin(o_) * sym.sin(o_/2)**2   *   sym.sin(2*theta_)                      *    sym.sin(freqsinorder(n_, S_, N_, Ln_, opn_)[29]*t_ - 1*phi_ + 3*op_ + Lp_),
             -((3*G_*M_*R_**2)/(32*a_**3))   * e_ *       (1 + 3*sym.cos(2*o_))            *   (1 + 3*sym.cos(2*theta_))              *    sym.cos(freqsinorder(n_, S_, N_, Ln_, opn_)[30]*t_ ),
              ((3*G_*M_*R_**2)/(8*a_**3))   *             sym.sin(o_)**2                   *   (1 - 3*sym.cos(theta_)**2)             *    sym.cos(freqsinorder(n_, S_, N_, Ln_, opn_)[31]*t_ - 2*op_),
              ((3*G_*M_*R_**2)/(32*a_**3))  * e_ *        sym.sin(o_)**2                   *   (1 + 3*sym.cos(2*theta_))              *    sym.cos(freqsinorder(n_, S_, N_, Ln_, opn_)[32]*t_ - 2*op_),
             ((21*G_*M_*R_**2)/(16*a_**3))  * e_ *        sym.sin(o_)**2                   *   (1 - 3*sym.cos(theta_)**2)             *    sym.cos(freqsinorder(n_, S_, N_, Ln_, opn_)[33]*t_ - 2*op_),
             -((G_*M_*R_**2)/(32*a_**3))    *             (1 + 3*sym.cos(2*o_))            *   (1 + 3*sym.cos(2*theta_))              *    sym.cos(freqsinorder(n_, S_, N_, Ln_, opn_)[34]*t_ )]


############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
#                                 Create the Satellite Class                                               #
class satellite:                                                                                           #

    def __init__(self,rin=[],pin=[],uin=[],nin=[],G=0,M=0,a=0,e=0,o=0,op=0,on=0,Ω=0,NSRn=0,libA=0,libp=0,libn=0):
        
        
        #All The Inputs
        self.rin = rin
        self.pin = pin
        self.uin = uin
        self.nin = nin

        self.G = G
        self.M = M
        self.a = a

        self.e = e



        self.o    = o
        self.op   = op
        self.on   = on
        

        self.NSRn    = NSRn


        self.libA   = libA
        self.libp   = libp
        self.libn   = libn


        # Body Radius
        self.R = self.rin[-1]

        # Body Mass
        self.mass = (4*np.pi/3) * pin[0]*rin[0]**3
        for i in range(1,len(rin)):
            self.mass += (4*np.pi/3) * pin[i]*(rin[i]**3 - rin[i-1]**3)
         
        # Mean Motion
        self.n = np.sqrt(self.G *(self.mass + self.M) / self.a**3)

        # grav accel
        self.g = self.G*self.mass / self.R**2


        if Ω == 0:
            self.Ω = self.n
        else:
            self.Ω = Ω
        


        self.freqlist = freqsinorder(self.n, self.Ω, self.NSRn, self.libn, self.on)
        self.potslist = potsinorder(self.G, self.M, self.R, self.a, self.e, self.o, self.op, self.libA, self.libp, self.n, self.Ω, self.NSRn, self.libn, self.on, t ,θ, φ)


        self.freqlist2 = []
        self.potslist2 = []
        self.static = []
        for i in range(len(self.freqlist)):
            if self.potslist[i] == 0:
                continue
            else:
                if self.freqlist[i] == 0:
                    self.static.append(self.potslist[i])
                else:
                    self.potslist2.append(self.potslist[i])
                    self.freqlist2.append(self.freqlist[i])

        self.static = [sum(self.static)]


        self.lovelist = []
        for i in self.freqlist2:
            self.lovelist.append(lovefunc(2, self.rin, self.pin, self.uin, self.nin, i))




        self.tt = 0
        self.pp = 0
        self.tp = 0
        for i in range(len(self.freqlist2)):

            #Surface Rigidity
            if self.nin[-1] == 0:
                self.uR = self.uin[-1]
            else:
                self.uR =  self.uin[-1] / (1-1j*(self.uin[-1]/(self.nin[-1]*self.freqlist2[i])))
            
            POT = self.potslist2[i]
            h = self.lovelist[i][0]
            l = self.lovelist[i][1]

            self.tt += ((2.0*self.uR)/(self.R*self.g)) * ((3.0*h -  6.0*l)*POT + l*sym.diff(POT,θ,θ))
            self.pp += ((2.0*self.uR)/(self.R*self.g)) * ((3.0*h - 12.0*l)*POT - l*sym.diff(POT,θ,θ))
            self.tp += ((2.0*self.uR)/(self.R*self.g)) *l * sym.csc(θ)*(sym.diff(POT,θ,φ) - sym.cot(θ)*sym.diff(POT,φ))

        self.ttR = sym.re(self.tt)
        self.ppR = sym.re(self.pp)
        self.tpR = sym.re(self.tp)       
        
        self.PC1  = (1/2) * (self.ttR + self.ppR + sym.sqrt(4*self.tpR**2 + (self.ttR-self.ppR)**2))
        self.PC2  = (1/2) * (self.ttR + self.ppR - sym.sqrt(4*self.tpR**2 + (self.ttR-self.ppR)**2))
        self.PCΨ  = (1/2) * sym.atan( (2*self.tpR)/(self.ttR-self.ppR))
        self.PCΨ2 = (1/2) * sym.atan2((2*self.tpR),(self.ttR-self.ppR))


