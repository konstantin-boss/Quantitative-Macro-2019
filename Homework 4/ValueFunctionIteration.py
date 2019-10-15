import numpy as np
import sympy as sp
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import seaborn as sb
import math
import quantecon as qec
from interpolation import interp
from numba import njit, prange
from quantecon.optimize.scalar_maximization import brent_max
import time

## Setting the parameter and steady state values
y = 1   # Normalization
h = 1
k = 1   # From production function
theta = 0.679
beta = 0.988
delta = 0.013
kappa = 0 # I block the labor part in utility for now
nu = 2.0

## Function definitions
# Production function
def outp(k,h):
    prod = pow(k,1-theta)*pow(h,theta)
    return(prod)

# Utility Function
def u(k1,k2,h):
    c = outp(k1,h)+(1-delta)*k1-k2
    ut = np.log(c)-kappa*pow(h,1+1/nu)/(1+1/nu)
    return(ut)

## Steady state variables
kss = ((1/beta-1+delta)/(1-theta))**(-1/theta)
i_ss = delta*kss
c_ss = outp(kss,h)-i_ss

## Discretization of the state space defined by k
kmin = 0.5*kss # Not too close to zero
kmax = 1.5*kss # close to steady state
nval = 100 # Number of values in the grid(well, then squared)
K = np.linspace(kmin, kmax, nval)

## Initial guess
V0 = np.zeros(nval)

## Definition of the return matrix
omega = -10000 # A very negative value to penalize negative consumption values
M_ret = np.empty([nval,nval])
for i,k1 in enumerate(K):
    for j,k2 in enumerate(K):
        if k2>outp(k1,h)+(1-delta)*k1:
            M_ret[i,j] = omega
        else:
            M_ret[i,j] = u(k1,k2,h)

'''## Computing the matrix Chi
M_chi = np.empty([nval,nval])
for i in range(nval):
    for j in range(nval):
        M_chi[i,j] = M_ret[i,j]+beta*V0[j]

## Compute the updated value function
V_upd = np.empty(nval)
for i in range(nval):
    V_upd[i] = max(M_chi[i])'''

## Putting it all together
epsilon = 0.005 # An arbitrary stopping parameter
V_s = V0
V_upd =np.ones(nval) # I just do not want the while loop to stop immediately
counter = 0
time1 = time.time()
g_pol = np.empty(nval) # Store vector for optimal decisions rule
while max(V_upd - V_s)>epsilon:
    counter +=1
    if counter == 1:
        V_s = V0
    else:
        V_s = V_upd

    M_chi = np.empty([nval,nval])
    for i in range(nval):
        for j in range(nval):
            M_chi[i,j] = M_ret[i,j]+beta*V_s[j]

    V_upd = np.empty(nval)
    for i in range(nval):
        V_upd[i] = max(M_chi[i])
        g_pol[i] = np.argmax(M_chi[i])

elapsed = round(time.time() - time1,2)
print('The number of iterations is',counter, 'and the time elapsed is', elapsed,'seconds')

## Checking the bounds
print(g_pol) # If 1 or nval-1 appear in the g_pol vector, bounds are too tight

## Exercise 1b) - 1f) - Iterations taking into account monotonicity of optimal decisions rule
epsilon = 0.005 # An arbitrary stopping parameter
V_s = V0
V_upd =np.ones(nval) # I just do not want the while loop to stop immediately
counter2 = 0 # To count the number of iterations
time2 = time.time() # To clock the time it takes for the iteration to converge
g_pol = np.empty(nval) # Store vector for optimal decisions rule

while max(V_upd - V_s)>epsilon:
    counter2 +=1
    if counter2 == 1: # In first iteration use initial guess
        V_s = V0
    else:
        V_s = V_upd # In following iterations update the guess

    M_chi = np.empty([nval,nval])

    for i in range(nval):
        for j in range(nval):
                if j>=1:
                    M_chi[i,j] = M_ret[i,j]+beta*V_s[j]
                    if M_chi[i,j] < M_chi[i,j-1]: # Here I exploit the concavity property, but I do not know how to break out of the inner loop so I just write zeros
                        M_chi[i,j] = 0
                else: # When j==0, we cannot exloit concavity yet as there is no comparison value
                    M_chi[i,j] = M_ret[i,j]+beta*V_s[j]


    V_upd = np.empty(nval)
    for i in range(nval):
        V_upd[i] = max(M_chi[i])
        g_pol[i] = K[np.unravel_index(np.argmax(M_chi[i,:], axis=None), M_chi[i,:].shape)[0]]


elapsed2 = round(time.time() - time2,2)
print('The number of iterations is',counter2, 'and the time elapsed is', elapsed2,'seconds')
plt.plot(K, V_upd)
plt.xlabel('Capital')
plt.ylabel('Value')
plt.show()
