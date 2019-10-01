import sympy as sy
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.misc import derivative
import cmath

def f(x):
	return x**0.321

x_grid = np.linspace(0,4,100)
fun = f
x0 = 1


# Taylor approximation at x0
def taylor(fun, x0, x, n):
	h = x - x0
	fx = np.empty(n+1)
	if n == 0:
		fx = fun(x0)
	else:
		for i in range(1, n+1):
			t = fun(x0)
			fx[0] = t
			fx[i] = derivative(fun, x0, n=i, order=25).real/math.factorial(i)*h**i
			fxsum = sum(fx)
	return(fxsum)

n = 1
j = -1
taylors1 = np.empty(len(x_grid))
for x in x_grid:
	j += 1
	taylors1[j] = taylor(fun,x0,x,n)

n = 2
j = -1
taylors2 = np.empty(len(x_grid))
for x in x_grid:
	j += 1
	taylors2[j] = taylor(fun,x0,x,n)

n = 5
j = -1
taylors5 = np.empty(len(x_grid))
for x in x_grid:
	j += 1
	taylors5[j] = taylor(fun,x0,x,n)

n = 20
j = -1
taylors20 = np.empty(len(x_grid))
for x in x_grid:
	j += 1
	taylors20[j] = taylor(fun,x0,x,n)

real = f(x_grid)

tay1, = plt.plot(taylors1, label = 'Order 1')
tay2, = plt.plot(taylors2, label = 'Order 2')
tay5, = plt.plot(taylors5, label ='Order 5')
tay20, = plt.plot(taylors20, label ='Order 20')
real, = plt.plot(real, label ='Actual function')
plt.title('Taylor approximations of power function')
plt.legend(handles = [tay1, tay2, tay5, tay20, real], loc='upper left')
plt.show()

## Exercise 2
def h(x):
    return (x+abs(x))/2

fun = h
x_grid = np.linspace(-2,6,100)
x0 = 2

n = 1
j = -1
taylors1 = np.empty(len(x_grid))
for x in x_grid:
	j += 1
	taylors1[j] = taylor(fun,x0,x,n)

n = 2
j = -1
taylors2 = np.empty(len(x_grid))
for x in x_grid:
	j += 1
	taylors2[j] = taylor(fun,x0,x,n)

n = 5
j = -1
taylors5 = np.empty(len(x_grid))
for x in x_grid:
	j += 1
	taylors5[j] = taylor(fun,x0,x,n)

n = 20
j = -1
taylors20 = np.empty(len(x_grid))
for x in x_grid:
	j += 1
	taylors20[j] = taylor(fun,x0,x,n)

real = h(x_grid)

tay1, = plt.plot(taylors1, label = 'Order 1')
tay2, = plt.plot(taylors2, label = 'Order 2')
tay5, = plt.plot(taylors5, label ='Order 5')
tay20, = plt.plot(taylors20, label ='Order 20')
real, = plt.plot(real, label ='Actual function')
plt.title('Taylor approximations of Ramp function')
plt.legend(handles = [tay1, tay2, tay5, tay20, real], loc='lower center')
plt.show()
