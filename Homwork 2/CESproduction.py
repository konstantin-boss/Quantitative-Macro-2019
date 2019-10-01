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
w = f_kh(nodes1[:,None],nodes1[None,:])
w = np.matrix(w)



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

thetas = coeff(nodes,d,w)
approx3 = cheby_approx(x,y,thetas,d)


## Plotting the true function and approximations
X, Y = np.meshgrid(nodes1, nodes2)

# Real figure
real = f_kh(X, Y)
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, real, cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.title('The real CES production function')
ax.set_xlabel('k')
ax.set_ylabel('h')
ax.set_zlabel('f_kh');
plt.show()

# Approximation degree 3 plot
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, approx3, cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.title('Approximation of degree 3')
ax.set_xlabel('k')
ax.set_ylabel('h')
ax.set_zlabel('f_kh');
plt.show()

# Approximation of degree 10
d = 10

thetas = coeff(nodes,d,w)
approx10 = cheby_approx(x,y,thetas,d)
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, approx10, cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.title('Approximation of degree 10')
ax.set_xlabel('k')
ax.set_ylabel('h')
ax.set_zlabel('f_kh')
plt.show()

# Approximation of degree 15
d = 15

thetas = coeff(nodes,d,w)
approx15 = cheby_approx(x,y,thetas,d)
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, approx15, cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.title('Approximation of degree 15')
ax.set_xlabel('k')
ax.set_ylabel('h')
ax.set_zlabel('f_kh');
plt.show()

# Plotting the exact isoquants of the real function
percentiles = np.array([5,10,25,50,75,90,95])
j = -1
levels = np.empty(len(percentiles))
for p in percentiles:
    j += 1
    levels[j] = np.percentile(real,p)

plt.contour(X,Y,real,levels)
plt.title('Isoquants for selected percentiles of real CES production')
plt.show()

## Plotting the isoquants of the approximations and the errors
# Approximation 3
j = -1
levels3 = np.empty(len(percentiles))
for p in percentiles:
    j += 1
    levels3[j] = np.percentile(approx3,p)

plt.contour(X,Y,approx3,levels3)
plt.title('Isoquants for selected percentiles of order 3 approximation')
plt.show()

errorlevels3 = abs(levels - levels3)
plt.plot(errorlevels3)
plt.title('Errors between percentiles of real and order 3')
plt.show()

# Approximation 10
j = -1
levels10 = np.empty(len(percentiles))
for p in percentiles:
    j += 1
    levels10[j] = np.percentile(approx10,p)

plt.contour(X,Y,approx10,levels10)
plt.title('Isoquants for selected percentiles of order 10 approximation')
plt.show()

errorlevels10 = abs(levels - levels10)
plt.plot(errorlevels10)
plt.title('Errors between percentiles of real and order 10')
plt.show()

# Approximation 15
j = -1
levels15 = np.empty(len(percentiles))
for p in percentiles:
    j += 1
    levels15[j] = np.percentile(approx15,p)

plt.contour(X,Y,approx15,levels10)
plt.title('Isoquants for selected percentiles of order 15 approximation')
plt.show()

errorlevels15 = abs(levels - levels15)
plt.title('Errors between percentiles of real and order 15')
plt.plot(errorlevels15)
plt.show()


## Plotting the approximation errors for all three approximations
errors3 = abs(real - approx3)
errors10 = abs(real - approx10)
errors15 = abs(real - approx15)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, errors3, cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.title('Errors between real and order 3 approximation in fval')
ax.set_xlabel('k')
ax.set_ylabel('h')
ax.set_zlabel('f_kh');
plt.show()

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, errors10, cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.title('Errors between real and order 10 approximation in fval')
ax.set_xlabel('k')
ax.set_ylabel('h')
ax.set_zlabel('f_kh');
plt.show()


fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, errors15, cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.title('Errors between real and order 10 approximation in fval')
ax.set_xlabel('k')
ax.set_ylabel('h')
ax.set_zlabel('f_kh');
plt.show()
