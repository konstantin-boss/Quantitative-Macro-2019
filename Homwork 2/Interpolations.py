import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.misc import derivative

x_lims=[-1,1]

def f(x):
    return np.exp(1/x)

def g(x):
    return 1/(1+25*x**2)

def h(x):
    return (x+abs(x))/2

def even_space(a,b,n,fun):
    nodes_ev = np.empty(n)
    for i in range(1,n+1):
        nodes_ev[i-1] = a + (i-1)/(n-1)*(b-a)
    fvals = fun(nodes_ev)
    return(nodes_ev,fvals)

def cheby_nodes(n,fun):
    nodes_cheb = np.empty(n)
    for j in range(1,n+1):
        nodes_cheb[j-1] = np.cos((2*j-1)/(2*n)*np.pi)
    fvals = fun(nodes_cheb)
    return(nodes_cheb, fvals)

def xmatrix(nodes):
    order = len(nodes)
    xmat = np.empty([order,order])
    xmat[:,0] = 1
    l =0
    for k in nodes:
        l+=1
        for i in range(1,order):
            xmat[l-1,i] = k**i
    return(xmat)


def thetas(xmatrix,fvals):
    inv = np.linalg.inv(xmatrix)
    thetas = np.dot(inv,fvals)
    return(thetas)

def monomial_values(x,theta):
    order = len(theta)
    mon = np.empty(order)
    tot = []
    for k in x:
        for i in range(order):
            mon[i] = theta[i]*k**i
            sum = np.sum(mon)
        tot.append(sum)
    return(tot)

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

def psi(i,fun,m): # i = degree of Cheby polynomial, m = number of nodes, fun = function
    z = cheby_nodes(m,fun)[0]
    psi_i = np.empty(m)
    l = -1
    for k in z:
        l += 1
        psi_i[l] = cheby_basis(i,k)[i]
    return(psi_i)

def numerator(i,fun,m):
    fv = cheby_nodes(m,fun)[1]
    psii = psi(i,fun,m)
    numer = np.dot(fv,psii)
    return(numer)

def denominator(i,fun,m):
    denom = np.dot(psi(i,fun,m), psi(i,fun,m))
    return(denom)

def theta_cheb(i,fun,m): # finds the thetas for a Cheby polynomial of order i fiven a function
    t = i+1
    thet = np.empty(t)
    for p in range(t):
        thet[p] = numerator(p,fun,m)/denominator(p,fun,m)
    return(thet)

def cheby_approx(i,fun,x_grid,a,b,m):
    l = len(x_grid)
    fx = np.empty(l)
    j = -1
    thetas = theta_cheb(i,fun,m) # dimension i (degree, e.g.2)+1 to include 0
    for k in x_grid:
        j += 1
        psy = cheby_basis(i,(2*(k-a)/(b-a)-1))
        t = np.dot(psy, thetas)
        fx[j] = t
    return(fx)

### Exercise 3.1 ###
## Ramp function ##
orders = [4,6,11]
polys =[]
x = np.linspace(-1,1,50)

for z in orders:
    even_nodes = even_space(-1,1,z,h)
    nodes = even_nodes[0]
    fvals = even_nodes[1]
    xmatrics = xmatrix(nodes)
    theta = thetas(xmatrics,fvals)
    cub = monomial_values(x,theta)
    polys.append(cub)

fig = plt.figure()
fig.suptitle('Ramp function, evenly-spaced nodes, cubic, and monomials')
ax1 = plt.subplot(211)
real_h =h(x)
cubic, = plt.plot(polys[0], label='Cubic')
five, = plt.plot(polys[1], label='Degree 5')
ten, = plt.plot(polys[2], label='Degree 10')
ramp_real, = plt.plot(real_h, label ='ramp(x)')
plt.title('The ramp function with evenly spaced nodes')
plt.legend(handles = [cubic, five, ten, ramp_real], loc='upper left')


errors3 = abs(real_h - polys[0])
errors5 = abs(real_h - polys[1])
errors10 = abs(real_h - polys[2])

plt.subplot(212, sharex=ax1)
er3, = plt.plot(errors3, label = 'Cubic')
er5, = plt.plot(errors5, label ='Monomial 5')
er10, = plt.plot(errors10, label ='Monomial 10')
plt.title('Errors Ramp, Evenly spaced, Cubic and Monomials')
plt.legend(handles = [er3, er5, er10], loc='upper center')
plt.show()

## Runge function
orders = [4,6,11]
polys =[]
x = np.linspace(-1,1,50)
real_g =g(x)
for z in orders:
    even_nodes = even_space(-1,1,z,g)
    nodes = even_nodes[0]
    fvals = even_nodes[1]
    xmatrics = xmatrix(nodes)
    theta = thetas(xmatrics,fvals)
    cub = monomial_values(x,theta)
    polys.append(cub)

fig = plt.figure()
fig.suptitle('Runge function, evenly-spaced nodes, cubic, and monomials')
ax1 = plt.subplot(211)
cubic, = plt.plot(polys[0], label='Cubic')
five, = plt.plot(polys[1], label='Degree 5')
ten, = plt.plot(polys[2], label='Degree 10')
ramp_real, = plt.plot(real_g, label ='Runge(x)')
plt.title('The Runge function with evenly spaced nodes')
plt.legend(handles = [cubic, five, ten, ramp_real], loc='upper center')

errors3 = abs(real_g - polys[0])
errors5 = abs(real_g - polys[1])
errors10 = abs(real_g - polys[2])

plt.subplot(212, sharex=ax1)
er3, = plt.plot(errors3, label = 'Cubic')
er5, = plt.plot(errors5, label ='Monomial 5')
er10, = plt.plot(errors10, label ='Monomial 10')
plt.title('Errors Runge, Evenly spaced, Cubic and Monomials')
plt.legend(handles = [er3, er5, er10], loc='upper center')
plt.show()

### Exercise 3.2
## Ramp function
orders = [4,6,11]
polys =[]
x = np.linspace(-1,1,50)
for z in orders:
    cheb_nodes = cheby_nodes(z,h)
    nodes = cheb_nodes[0]
    fvals = cheb_nodes[1]
    xmatrics = xmatrix(nodes)
    theta = thetas(xmatrics,fvals)
    cub = monomial_values(x,theta)
    polys.append(cub)

real_h =h(x)
fig = plt.figure()
fig.suptitle('Ramp function, Cheby nodes, cubic, and monomials')
ax1 = plt.subplot(211)
real_h =h(x)
cubic, = plt.plot(polys[0], label='Cubic')
five, = plt.plot(polys[1], label='Degree 5')
ten, = plt.plot(polys[2], label='Degree 10')
ramp_real, = plt.plot(real_h, label ='ramp(x)')
plt.title('The ramp function with Cheby nodes')
plt.legend(handles = [cubic, five, ten, ramp_real], loc='upper left')

errors3 = abs(real_h - polys[0])
errors5 = abs(real_h - polys[1])
errors10 = abs(real_h - polys[2])

plt.subplot(212, sharex=ax1)
er3, = plt.plot(errors3, label = 'Cubic')
er5, = plt.plot(errors5, label ='Monomial 5')
er10, = plt.plot(errors10, label ='Monomial 10')
plt.title('Errors Ramp, Cheby nodes, Cubic and Monomials')
plt.legend(handles = [er3, er5, er10], loc='upper left')
plt.show()

## The Runge function
orders = [4,6,11]
polys =[]
x = np.linspace(-1,1,50)
for z in orders:
    cheb_nodes = cheby_nodes(z,g)
    nodes = cheb_nodes[0]
    fvals = cheb_nodes[1]
    xmatrics = xmatrix(nodes)
    theta = thetas(xmatrics,fvals)
    cub = monomial_values(x,theta)
    polys.append(cub)

real_g =g(x)
fig = plt.figure()
fig.suptitle('Runge function, Cheby nodes, cubic, and monomials')
ax1 = plt.subplot(211)
cubic, = plt.plot(polys[0], label='Cubic')
five, = plt.plot(polys[1], label='Degree 5')
ten, = plt.plot(polys[2], label='Degree 10')
ramp_real, = plt.plot(real_g, label ='Runge(x)')
plt.title('The Runge function with Cheby nodes')
plt.legend(handles = [cubic, five, ten, ramp_real], loc='upper left')

errors3 = abs(real_g - polys[0])
errors5 = abs(real_g - polys[1])
errors10 = abs(real_g - polys[2])

plt.subplot(212, sharex=ax1)
er3, = plt.plot(errors3, label = 'Cubic')
er5, = plt.plot(errors5, label ='Monomial 5')
er10, = plt.plot(errors10, label ='Monomial 10')
plt.title('Errors Runge, Cheby nodes, Cubic and Monomials')
plt.legend(handles = [er3, er5, er10], loc='upper left')
plt.show()

### Exercise 3.3
## Chebyshev interpolation nodes and Chebyshev polynomial of orders 3,5 and 10
a = -1 # lower interval bound
b = 1 # upper interval bound

# The ramp function
x_grid = np.linspace(-1,1,50)
real_ramp = h(x_grid)
m = 4
i = 3
ch_app_3 = cheby_approx(i,h,x_grid,a,b,m)
i = 5
m = 6
ch_app_5 = cheby_approx(i,h,x_grid,a,b,m)
i = 10
m = 11
ch_app_10 = cheby_approx(i,h,x_grid,a,b,m)

fig = plt.figure()
fig.suptitle('Ramp function, Cheby nodes, Cheby polynomials')
ax1 = plt.subplot(211)
app3, = plt.plot(ch_app_3, label = 'Cubic Chevy')
app5, = plt.plot(ch_app_5, label = '5 Chevy')
app10, = plt.plot(ch_app_10, label = '10 Chevy')
realr, = plt.plot(real_ramp,label = 'Ramp')
plt.title('Ramp function, Cheby nodes, Cheby polys')
plt.legend(handles = [app3, app5, app10, realr], loc='upper left')
plt.legend(handles = [cubic, five, ten, ramp_real], loc='upper left')

errors3 = abs(real_h - ch_app_3)
errors5 = abs(real_h - ch_app_5)
errors10 = abs(real_h - ch_app_10)

plt.subplot(212, sharex=ax1)
er3, = plt.plot(errors3, label = 'Chevy Cubic')
er5, = plt.plot(errors5, label ='Chevy 5')
er10, = plt.plot(errors10, label ='Chevy 10')
plt.title('Errors Ramp, Cheby Nodes, Chevy Polys')
plt.legend(handles = [er3, er5, er10], loc='upper left')
plt.show()

# The Runge function
real_runge = g(x_grid)
i = 3
m = 4
ch_app_3 = cheby_approx(i,g,x_grid,a,b,m)
i = 5
m = 6
ch_app_5 = cheby_approx(i,g,x_grid,a,b,m)
i = 10
m = 11
ch_app_10 = cheby_approx(i,g,x_grid,a,b,m)

fig = plt.figure()
fig.suptitle('Runge function, Cheby nodes, Cheby polynomials')
ax1 = plt.subplot(211)
app3, = plt.plot(ch_app_3, label = 'Cubic Chevy')
app5, = plt.plot(ch_app_5, label = '5 Chevy')
app10, = plt.plot(ch_app_10, label = '10 Chevy')
realr, = plt.plot(real_runge,label = 'Runge')
plt.title('Runge function, Cheby nodes, Cheby polys')
plt.legend(handles = [app3, app5, app10, realr], loc='upper left')

plt.subplot(212, sharex=ax1)
errors3 = abs(real_g - ch_app_3)
errors5 = abs(real_g - ch_app_5)
errors10 = abs(real_g - ch_app_10)
er3, = plt.plot(errors3, label = 'Chevy Cubic')
er5, = plt.plot(errors5, label ='Chevy 5')
er10, = plt.plot(errors10, label ='Chevy 10')
plt.title('Errors Runge, Cheby Nodes, Chevy Polys')
plt.legend(handles = [er3, er5, er10], loc='upper left')
plt.show()
