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


## Function definitions
# Production function
def outp(k,h):
    prod = pow(k,1-theta)*pow(h,theta)
    return(prod)

# Utility Function
def u(k1,k2,h):
    c = outp(k1,h)+(1-delta)*k1-k2
    ut = np.log(c)
    return(ut)


## Steady state variables
kss = ((1/beta-1+delta)/(1-theta))**(-1/theta)
i_ss = delta*kss
c_ss = outp(kss,h)-i_ss

## Discretization of the state space defined by k
kmin = 0.5*kss # Not too close to zero
kmax = 1.5*kss # close to steady state
nval = 200 # Number of values in the grid(well, then squared)
K = np.linspace(kmin, kmax, nval)



## Initial guess
V0 = np.zeros(nval)

epsilon = 0.005 # An arbitrary stopping parameter
V_s = V0
V_upd =np.ones(nval) # I just do not want the while loop to stop immediately
g_pol = np.zeros(nval) # Store vector for optimal decisions rule
omega = -10000 # A very negative value to penalize negative consumption values
counter = 0
time1 = time.time()

omega = -10000 # A very negative value to penalize negative consumption values
M_ret = np.empty([nval,nval])
for i,k1 in enumerate(K):
    for j,k2 in enumerate(K):
        if k2>outp(k1,h)+(1-delta)*k1:
            M_ret[i,j] = omega
        else:
            M_ret[i,j] = u(k1,k2,h)

## I run this code 50 times to get a stabel policy function and a first approximation of V
while counter<50:
    counter +=1
    if counter == 1:
        V_s = V0
    else:
        V_s = V_upd

    M_chi = np.empty([nval,nval])
    for i in range(nval):
        for j in range(nval):
            t = int(g_pol[i])  # If the polcy of the iteration before says argmax is at t, only consider positions in Chi that are larger than t
            if t+j<nval:
                M_chi[i,t+j] = M_ret[i,j+t]+beta*V_s[t+j]
            else:
                continue

        g_pol[i] = np.argmax(M_chi[i])

    V_upd = np.empty(nval)
    for i in range(nval):
        V_upd[i] = max(M_chi[i])


#### Now for the policy function iteration
V1 = V_upd
V = V_upd
new_g_pol = np.empty(len(K))

iterations = [50,20,10,5]

for T in iterations:
    for it in range(T):
        for i,k in enumerate(K):
            kpol = int(g_pol[i])
            if K[kpol]<outp(k,h)+(1-delta*k):
                k2 = K[kpol]
                V_upd[i] = u(k,k2,h)+beta*V[kpol]
            else:
                V_upd[i] = omega
        if max(V_upd-V)<epsilon:
            V_final = V_upd
            print('The value function is approximated as V_final.')
            break
        else:
            V = V_upd

    for i,k1 in enumerate(K):
        for j, k2 in enumerate(K):
            M_chi[i,j] = M_ret[i,j]+beta*V[j]
        new_g_pol[i] = np.argmax(M_chi[i])
        V_upd[i] = np.nanmax(M_chi[i])
    if max(V_upd-V)<epsilon:
        V_final = V_upd
        print('The value function is approximated as V_final.')
        break
    else:
        V = V_upd
    if (new_g_pol==g_pol).all():
        break
    else:
        new_g_pol=g_pol

plt.plot(V_final)
plt.show()








































## The return matrix
'''omega = -10000 # A very negative value to penalize negative consumption values
M_ret = np.empty([nval,nval])
for i,k1 in enumerate(K):
    for j,k2 in enumerate(K):
        if k2>outp(k1,h)+(1-delta)*k1:
            M_ret[i,j] = omega
        else:
            M_ret[i,j] = u(k1,k2,h)

##
epsilon = 0.005 # An arbitrary stopping parameter
V_s = V0
V_upd =np.ones(nval) # I just do not want the while loop to stop immediately
counter2 = 0 # To count the number of iterations
time2 = time.time() # To clock the time it takes for the iteration to converge
g_pol = np.empty(nval) # Store vector for optimal decisions rule
gplo = np.empty(nval)
while max(V_upd - V_s)>epsilon:
    counter2 +=1
    if counter2 == 1: # In first iteration use initial guess
        V_s = V0
    else:
        V_s = V_upd # In following iterations update the guess

    M_chi = np.empty([nval,nval])

    for i in range(nval):
        for j in range(nval):
                M_chi[i,j] = M_ret[i,j]+beta*V_s[j]


    V_upd = np.empty(nval)
    for i in range(nval):
        V_upd[i] = max(M_chi[i])
        gplo[i] = np.argmax(M_chi[i,:])
        g_pol[i] = K[np.unravel_index(np.argmax(M_chi[i,:], axis=None), M_chi[i,:].shape)[0]]


elapsed2 = round(time.time() - time2,2)

plt.plot(K, V_upd)
plt.xlabel('Capital')
plt.ylabel('Value')
plt.show()

print(gplo)'''
