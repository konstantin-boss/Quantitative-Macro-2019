import sympy as sy
import numpy as np
import matplotlib.pyplot as plt
import math

# Taylor approximation at x0 of the function 'function'
def taylor(function, x0, n):
    i = 0
    p = 0
    while i <= n:
        p = p + (function.diff(x, i).subs(x, x0))/(math.factorial(i))*(x - x0)**i
        i += 1
    return p

# Define the parameters and the function to approximate
x = sy.Symbol('x')
f = x**0.321
def fx(x):
    return x**0.321

x0 = 1
a = 0
b = 4
N = 50
x_grid = np.linspace(a,b,N)
degrees = [1,2,5,10]

f_taylor = []
for i in degrees:
    f_approx = taylor(f,x0,i)
    for k in x_grid:
        pt = f_approx.subs(x, k)
        f_taylor.append(pt)
taylor1 = f_taylor[0:49]
taylor2 = f_taylor[50:99]
taylor5 = f_taylor[100:149]
taylor20 = f_taylor[150:199]
real = fx(x_grid)

tay1, = plt.plot(taylor1, label='Degree 1')
tay2, = plt.plot(taylor2, label='Degree 2')
tay5, = plt.plot(taylor5, label='Degree 5')
tay20, = plt.plot(taylor20, label='Degree 20')
pow_real, = plt.plot(real, label ='Power function')
plt.title('A power function with Taylor approximations')
plt.legend(handles = [tay1, tay2, tay5, tay20, pow_real], loc='upper left')
plt.ylim(0, 2)
plt.show()


## Exercise 2
# Function definition
x = sy.Symbol('x', real=True)
def g(x):
	return 0.5*(x + abs(x))

# Setting parameters
x0 = 2
a = -2
b = 6
N = 50
x_grid = np.linspace(a,b,N)
degrees = [1,2,5,10]

# Approximations with different degrees
f_taylor = []
for i in degrees:
    f_approx = taylor(g(x),x0,i)
    for k in x_grid:
        pt = f_approx.subs(x, k)
        f_taylor.append(pt)
taylor1 = f_taylor[0:49]
taylor2 = f_taylor[50:99]
taylor5 = f_taylor[100:149]
taylor20 = f_taylor[150:199]

real = g(x_grid) # real functon
tay1, = plt.plot(taylor1, label='Degree 1')
tay2, = plt.plot(taylor2, label='Degree 2')
tay5, = plt.plot(taylor5, label='Degree 5')
tay20, = plt.plot(taylor20, label='Degree 20')
ramp_real, = plt.plot(real, label ='ramp(x)')
plt.title('The ramp function with Taylor approximations')
plt.legend(handles = [tay1, tay2, tay5, tay20, ramp_real], loc='upper left')
plt.ylim(-1, 2)
plt.show()
