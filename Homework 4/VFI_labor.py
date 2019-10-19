import numpy as np
import sympy as sp
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import quantecon as qe
from scipy.optimize import brentq
from scipy.optimize import minimize

### Based on code by ALbert Rodriguez Sala ###

y = 1   # Normalization
k = 1   # From production function
theta = 0.679
beta = 0.988
delta = 0.013
kappa = 5.24
nu = 2.0

## Function definitions
# Production function
def outp(k,h):
    prod = pow(k,1-theta)*pow(h,theta)
    return(prod)

# Utility Function
def u(k1,h):
    c = outp(k1,h)+(1-delta)*k1-K
    ut = np.log(c)-kappa*pow(h,1+1/nu)/(1+1/nu)
    return(ut)

## Steady state variables
kss = ((1/beta-1+delta)/(1-theta))**(-1/theta)
i_ss = delta*kss

def steadys(x): #x[0] is c_ss and x[1] is h_ss
    f1 = x[0] - outp(kss,x[1])+i_ss
    f2 = pow(theta/(kappa*x[0])*pow(x[1],1-theta),1/(1/nu-theta+1))
    return[f1,f2]

steady = fsolve(steadys, [1,1])
css = steady[0]
hss = steady[1]

## Discretization of the state space defined by k
kmin = 0.5*kss # Not too close to zero
kmax = 1.5*kss # close to steady state
hmin = 0.5*hss # Not too close to zero
hmax = 1.5*hss # close to steady state
nval = 100 # Number of values in the grid(well, then squared)
K = np.linspace(kmin, kmax, nval)
H = np.linspace(hmin, hmax, nval)

##
V0=np.zeros(nval)

## Definition of the Chi matrix
g_pol = np.zeros(nval) # Store vector for optimal decisions rule for capital
h_pol = np.zeros(nval) # Store vector for optimal decisions rule for labor
c_pol = np.zeros(nval)

time1 = time.time()

def bellman_operator(V,return_policies=False):
    V_upd = np.zeros((nval))
    M_chi = np.empty([nval,nval,nval])
    for i, k in enumerate(K):
        for j,h in enumerate(H):
            M_chi[i,:,j] =  u((outp(k,h) +(1-delta)*k - K),K,h) +beta*V

        V_upd[i] = np.nanmax(M_chi[i,:,:])
        g_pol[i] = K[np.unravel_index(np.argmax(M_chi[i,:,:], axis=None), M_chi[i,:,:].shape)[0]]
        h_pol[i] = H[np.unravel_index(np.nanargmax(M_chi[i,:,:], axis=None), M_chi[i,:,:].shape)[1]]
        c_pol[i] =  outp(k,h_pol[i]) +(1-delta)*k - g_pol[i]

    if return_policies==True:
        return V_upd, g_pol, h_pol, c_pol
    else:
        return V_upd

V = qe.compute_fixed_point(bellman_operator, V0, max_iter=2000, error_tol=0.001, print_skip=20)
V, g_k, g_h, g_c = bellman_operator(V, return_policies=True)

elapsed = round(time.time() - time1,2)
print('Time elapsed is', elapsed,'seconds')

## Checking the bounds
print(V) # If 1 or nval-1 appear in the g_pol vector, bounds are too tight
plt.plot(K,V)
plt.show()
