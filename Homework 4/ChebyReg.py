import numpy as np
import matplotlib.pyplot as plt
import math
import time

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

# Steady state
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



def cheby_nodes(n):
    nodes_cheb = np.empty(n)
    for j in range(1,n+1):
        nodes_cheb[j-1] = np.cos((2*j-1)/(2*n)*np.pi)
    return(nodes_cheb)

def cheby_basis(j,x):
    if j == 0:
        cheby_polyms = np.empty(1)
        cheby_polyms[0] = 1
    elif j == 1:
        cheby_polyms = np.empty(2)
        cheby_polyms[0] = 1
        cheby_polyms[1] = x
    else:
        t = j+1
        cheby_polyms = np.empty(t)
        cheby_polyms[0] = 1
        cheby_polyms[1] = x
        for i in range(2,t):
            cheby_polyms[i] = 2*x*cheby_polyms[i-1]-cheby_polyms[i-2]
    return(cheby_polyms)

def psi(i,m): # i = degree of Cheby polynomial, m = number of nodes, fun = function
    z = cheby_nodes(m)
    nodes1 = (z + 1)*0.5*(b-a)+a
    psi_i = np.empty(m)
    l = -1
    for k in nodes1:
        l += 1
        psi_i[l] = cheby_basis(i,k)[i]
    return(psi_i)

def numerator(i,fval,m):
    fv = fval
    psii = psi(i,m)
    numer = np.dot(fv,psii)
    return(numer)

def denominator(i,m):
    denom = np.dot(psi(i,m), psi(i,m))
    return(denom)

def theta_cheb(i,m): # finds the thetas for a Cheby polynomial of order i fiven a function
    t = i+1
    thet = np.empty(t)
    for p in range(t):
        thet[p] = numerator(p,fval,m)/denominator(p,m)
    return(thet)

def cheby_approx(i,x_grid,m):
    l = len(x_grid)
    fx = np.empty(l)
    j = -1
    thetas = theta_cheb(i,m) # dimension i (degree, e.g.2)+1 to include 0
    for k in x_grid:
        j += 1
        psy = cheby_basis(i,(2*(k-a)/(b-a)-1))
        t = np.dot(psy, thetas)
        fx[j] = t
    return(fx)


d = 3 # degree of Cheby polynomial
m = 50 # Number of nodes
nval = 50

epsilon = 0.5 # tolerance level
p = 0 # counting the iterations
V_upd = np.ones(nval) # starting value for updatable function
V_s = np.zeros(nval) # initial guess
K = np.linspace(kmin, kmax, nval)

g_pol =np.empty(nval)
omega = -10000 # A very negative value to penalize negative consumption values
M_ret = np.empty([nval,nval])
for i,k1 in enumerate(K):
    for j,k2 in enumerate(K):
        if k2>outp(k1,h)+(1-delta)*k1:
            M_ret[i,j] = omega
        else:
            M_ret[i,j] = u(k1,k2,h)

time1 = time.time()
counter= 0




while max(V_upd - V_s)>epsilon:
    counter +=1
    if counter == 1:
        fval = np.zeros(m)
    else:
        fval = V_upd
        V_s = V_upd

    M_chi = np.empty([nval,nval])
    for i in range(nval):
        for j in range(nval):
            M_chi[i,j] = M_ret[i,j]+beta*cheby_approx(d,K,m)[j]


        g_pol[i] = np.argmax(M_chi[i])

    V_upd = np.empty(nval)
    for i in range(nval):
        V_upd[i] = max(M_chi[i])


elapsed = round(time.time() - time1,2)
print('The number of iterations is',counter, 'and the time elapsed is', elapsed,'seconds')

## Checking the bounds
print(g_pol) # If 1 or nval-1 appear in the g_pol vector, bounds are too tight
