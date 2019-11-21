import numpy as np
import sympy as sp
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import time

## Setting the parameter and steady state values
theta = 0.679
beta = 0.988
yuck = 0.0000001
grdmax = 30
grdmin = np.sqrt(yuck)
varepsi = 0.01
g = 0.01
muepsi = -varepsi/2
el = 7
ngr = 20
curv = 3.0
r = 0.02
omega = -10000

## Function definitions

# Utility Function
def u(c):
    if abs(theta-1.0)<np.sqrt(yuck):
        ut = np.log(c)
    else:
        ut = 1/(1-theta)*pow(c,1-theta)
    return(ut)

# Making the grids
def makegrid(x1,x2,ngr,curv):
    grd = np.empty(ngr)
    scale=x2-x1
    grd[0] = x1
    grd[ngr-1] = x2
    for igr in range(1, ngr-1):
        grd[igr] = x1+scale*pow((igr-1.0)/(ngr-1.0),curv)
    return(grd)

np.random.seed(13)
x1 = makegrid(grdmin,grdmax,ngr,curv)
x2 = makegrid(grdmin,grdmax,ngr,curv)
err = np.random.normal(muepsi,np.sqrt(varepsi),el)
err = np.exp(err)



W = np.empty([el,ngr])
V0 = np.zeros(ngr*el)
pi = np.array([0.0006, 0.03075, 0.2401, 0.4571, 0.2401, 0.03075,0.0006])
epsilon = 0.005 # An arbitrary stopping parameter
V_s = V0
Vupd = np.ones(ngr*el)
xpol = np.empty(ngr*el)
counter = 0
time1 = time.time()

while max(Vupd - V_s)>epsilon:
    counter +=1
    if counter == 1:
        V_s = V0
    else:
        V_s = Vupd

    M = np.empty([ngr*el,ngr])
    for i,z1 in enumerate(err):
        for j,xc1 in enumerate(x1):
            for k,xc2 in enumerate(x2):
                if xc1-(xc2-z1)/(1-r)<np.sqrt(yuck): # rule out negative consumption
                    M[j+(len(x1)*i),k] = omega
                elif xc1-(xc2-z1)/(1-r) > xc1: # occasionally binding borrowing constraint
                    M[j+(len(x1)*i),k] = u(xc1)
                else:
                    M[j+(len(x1)*i),k] = u(xc1-(xc2-z1)/(1-r)*(1+g))

    for j in range(el):
        for i in range(ngr):
            W[j,i] = np.dot(pi,V_s[i::ngr])

    ## Compute Matrix Chi
    Chi = np.empty([ngr*el,ngr])
    for i in range(ngr*el):
        for j in range(ngr):
            k = 0
            Chi[i,j] = M[i,j]+beta*pow(1+g,1-theta)*W[k,j]
        if k < el:
            k += 1
        else:
            k = 0
        xpol[i] = np.argmax(Chi[i])

    for i in range(ngr*el):
        Vupd[i] = np.max(Chi[i])





elapsed = round(time.time() - time1,2)
print('The number of iterations is',counter, 'and the time elapsed is', elapsed,'seconds')

## Checking the bounds
v1 = Vupd[0:ngr]


ci = np.empty(ngr)
xpol = xpol.astype(int)
for d,c in enumerate(xpol[0:ngr]):
    ci[d] = (1+g)/(1+r)*(err[1]-x2[c])+(1+r)/(1+g)*x1[c]

 # If 1 or nval-1 appear in the g_pol vector, bounds are too tight
plt.plot(x1,v1)
plt.show()

plt.plot(x1,ci)
plt.show()
