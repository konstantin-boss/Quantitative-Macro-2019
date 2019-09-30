import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.misc import derivative
from mpl_toolkits.mplot3d import Axes3D as plt3d
from matplotlib import cm

alpha = 0.5
sig = 0.25
k = np.linspace(0,10,50)
h = np.linspace(0,10,50)

def f_kh(k,h):
    return ((1-alpha)*k**((sig-1)/sig)+alpha*h**((sig-1)/sig))**(sig/(sig-1))

def cheby_nodes(n):
    nodes_cheb = np.empty(n)
    for j in range(1,n+1):
        nodes_cheb[j-1] = np.cos((2*j-1)/(2*n)*np.pi)
    return(nodes_cheb)

## Step 1: Creating the nodes on [0,10] x [0,10]
n = 20                          # number of desired nodes
a = 0                           # lower bound for h and k
b = 10                          # upper bound for h and k
nodes = cheby_nodes(n)

## Step 2: Adjusting the nodes for the intervals
nodes1 = (nodes + 1)*0.5*(b-a)+a
nodes2 = (nodes + 1)*0.5*(b-a)+a

## Step 3: Evaluating f_kh at the nodes
X, Y = np.meshgrid(nodes1, nodes2)
w = f_kh(X, Y)



## Step 4: Computing the Chebyshev coefficients with degrees i = 3-15
def cheby_poly(d,x): # d=degree, x=nodes
    psi = []
    psi.append(np.ones(len(x)))
    psi.append(x)
    for i in range(1,d):
        p = 2*x*psi[i]-psi[i-1]
        psi.append(p)
    pol_d = np.matrix(psi[d])
    return pol_d

def coeff(x,d,w):
    t = len(x)
    thetas = np.empty((d+1)*(d+1))
    thetas.shape = (d+1,d+1)
    for i in range(d+1):
        for j in range(d+1):
            thetas[i,j] = (np.sum(np.array(w)*np.array((np.dot(cheby_poly(i,x).T,cheby_poly(j,x)))))
                          /np.array((cheby_poly(i,x)*cheby_poly(i,x).T)*(cheby_poly(j,x)*cheby_poly(j,x).T)))
    return thetas


def cheby_approx(x,y,thetas,d):
    f = []
    in1 = (2*(x-a)/(b-a)-1)
    in2 = (2*(y-a)/(b-a)-1)
    for u in range(d):
        for v in range(d):
                f.append(np.array(thetas[u,v])*np.array((np.dot(cheby_poly(u,in1).T,cheby_poly(v,in2)))))
    f_sum = sum(f)
    return f_sum

# Approximation for degree = 3

d = 3
x = nodes1
y = nodes2

thetas = coeff(x,d,w)
approx3 = cheby_approx(x,y,thetas,d)

## Plotting
X, Y = np.meshgrid(nodes1, nodes2)

# Real figure
real = f_kh(X, Y)
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, real, cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax.set_xlabel('k')
ax.set_ylabel('h')
ax.set_zlabel('f_kh');
plt.show()

# Approximation degree 3 plot
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, approx3, cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax.set_xlabel('k')
ax.set_ylabel('h')
ax.set_zlabel('f_kh');
plt.show()

# Approximation of degree 10
d = 10

thetas = coeff(x,d,w)
approx10 = cheby_approx(x,y,thetas,d)
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, approx10, cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax.set_xlabel('k')
ax.set_ylabel('h')
ax.set_zlabel('f_kh')
plt.show()

# Approximation of degree 15
d = 15

thetas = coeff(x,d,w)
approx15 = cheby_approx(x,y,thetas,d)
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, approx15, cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax.set_xlabel('k')
ax.set_ylabel('h')
ax.set_zlabel('f_kh');
plt.show()

fig = plt.figure()
error3 = abs(w - approx3)
ax.set_xlabel('k')
ax.set_ylabel('h')
ax.set_zlabel('f_kh')
ax.plot_surface(X, Y, error3, cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax.set(title='Error 3')
plt.show()
