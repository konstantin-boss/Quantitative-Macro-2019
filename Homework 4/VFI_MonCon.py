import numpy as np
import sympy as sp
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import time

import numpy as np
import sympy as sp
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
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
nval = 200 # Number of values in the grid(well, then squared)
K = np.linspace(kmin, kmax, nval)

## Initial guess
V0 = np.zeros(nval)

## Definition of the Chi matrix
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


while max(V_upd - V_s)>epsilon:
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
                M_chi[i,t+j] = M_ret[i,j+t]+beta*V_s[t+j] # Write Chi matrix only for positions larger than t
                if (j>=1 and M_chi[i,j+t]<M_chi[i,t+j-1]): # Check if among those positions we get that there is a preceding one that is larger
                    continue                                # If yes, skip the for loop for this j and go on, if no, do nothing

        g_pol[i] = np.argmax(M_chi[i])

    V_upd = np.empty(nval)
    for i in range(nval):
        V_upd[i] = max(M_chi[i])


elapsed = round(time.time() - time1,2)
print('The number of iterations is',counter, 'and the time elapsed is', elapsed,'seconds')

## Checking the bounds
print(g_pol) # If 1 or nval-1 appear in the g_pol vector, bounds are too tight


plt.plot(K, V_upd)
plt.xlabel('Capital')
plt.ylabel('Value')
plt.show()
