import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.misc import derivative
from mpl_toolkits.mplot3d import Axes3D as plt3d
from matplotlib import cm


theta = 0.679
beta = 0.988
delta = 0.013
h = 1
kappa = 0 # I block the labor part in utility for now
nu = 2.0

# Production function
def outp(k,h):
    prod = pow(k,1-theta)*pow(h,theta)
    return(prod)

# Utility Function
def u(k1,k2,h):
    c = outp(k1,h)+(1-delta)*k1-k2
    ut = np.log(c)-kappa*pow(h,1+1/nu)/(1+1/nu)
    return(ut)

def cheby_nodes(n):
    nodes_cheb = np.empty(n)
    for j in range(1,n+1):
        nodes_cheb[j-1] = np.cos((2*j-1)/(2*n)*np.pi)
    return(nodes_cheb)

# Steady states
kss = ((1/beta-1+delta)/(1-theta))**(-1/theta)
i_ss = delta*kss
c_ss = outp(kss,h)-i_ss
kmin = 0.5*kss # Not too close to zero
kmax = 1.5*kss # close to steady state
nval = 20 # Number of values in the grid(well, then squared)



## Step 1: Creating the nodes 
n = nval                          # number of desired nodes
a = kmin                          # lower bound for h and k
b = kmax                          # upper bound for h and k
nodes = cheby_nodes(n)

## Step 2: Adjusting the nodes for the intervals
nodes1 = (nodes + 1)*0.5*(b-a)+a
nodes2 = (nodes + 1)*0.5*(b-a)+a

## Initial guess V0
V0 = np.zeros(nval)

## Finding the Chebyshev polynomials
def cheby_poly(d,x): # d=degree, x=nodes
    psi = []
    psi.append(np.ones(len(x)))
    psi.append(x)
    for i in range(1,d):
        p = 2*x*psi[i]-psi[i-1]
        psi.append(p)
    pol_d = np.matrix(psi[d])
    return pol_d


def coeffs(V, d, x, theta): # V0 = initial guess, d = degree, x = nodes
    F = np.zeros(nval)
    for i,k in enumerate(x):
        for j in range(d+1):
            F[i] += theta[j]*cheby_poly(j,)
        F[i] = F[i]-V[i]

test = coeffs(V0,3,nodes1,np.ones(nval))
