import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import math

theta = 0.67    # Given
y = 1           # Normalization of output
kss1 = 4           # Since we have that k/y=4
delta = 1/16    # Since i/y=0.25 and i=delta*k in steady state
i = delta*kss1     # In steay state
h = 0.31        # Given
c = y-i         # Goods market clearing

## Part a)
# From the production function we get the following
z1 = ((1/kss1**(1-theta))**(1/theta))/h
print('The intitial steady state level of productivitiy is',z1)
print('The intitial steady state level of output is',y)
print('The intitial steady state level of consumption is',c)
print('The intitial steady state level of investment is',i)

## Determining beta
beta = 1/((1-theta)*kss1**-theta*(z1*h)**theta+1-delta)
print('The discount factor beta is',beta)

## Part b)
z2 = 2*z1       # Doubling z permanently
k2 = (z2*h)/(1/(beta*(1-theta))-(1-delta)/(1-theta))**(1/theta)
i2 = delta*k2
y2 = k2**(1-theta)*(h*z2)**theta
c2 = y2-i2
print(' ')
print('New steady state for double productivity:')
print('The new steady state level of capital',k2)
print('The new steady state level of productivitiy is',z2)
print('The new steady state level of output is',y2)
print('The new steady state level of consumption is',c2)
print('The new steady state level of investment is',i2)

## Part c)
# From the first order condition we get the following Euler equation which governs the dynamics
def euler(k0,k1,k2):
    return pow(k1,1-theta)*pow(z2*h,theta)+(1-delta)*k1-k2               \
    -beta*((1-theta)*pow(k1,-theta)*pow(z2*h,theta)+(1-delta))      \
    *(pow(k0,1-theta)*pow(z2*h,theta)-k1+(1-delta)*k0)

size = 100
def ktrans(z):
    F = np.zeros(size)
    z = z
    F[0] = euler(4,z[1],z[2])
    z[size-1] = k2
    F[size-2] = euler(z[size-3], z[size-2], z[size-1])
    for i in range(1,size-2):
        F[i] = euler(z[i],z[i+1],z[i+2])
    return F
z = np.ones(size)*4 # The initial steady state
k = fsolve(ktrans, z)
k[0] = 4

plt.plot(k)
plt.title('The transition path of capital between steady states')
plt.xlabel('Time')
plt.ylabel('Capital')
plt.show()


### Plotting the time paths for savings, consumption, and output. I skip labor as it is constant.
savingspath = np.empty(size)
# We lose one degree of freedom due to the t+1 structure of investment (=savings)
for i in range(size-1):
        savingspath[i] = k[i+1]-(1-delta)*k[i]
        savingspath[size-1] = savingspath[size-2]

outputpath = pow(k,1-theta)*pow(z2*h, theta)
consumptionpath = outputpath - savingspath

cons, = plt.plot(consumptionpath, label='Consumption')
sav, = plt.plot(savingspath, label='Savings')
outp,= plt.plot(outputpath,label='Output')
plt.title('Different transition paths without shock')
plt.legend(handles = [cons, sav, outp], loc='center right')
plt.xlabel('Time')
plt.ylabel('Levels')
plt.show()

### Unexpected shocks
ks = np.empty(size) # new ksurprise
ts = 10 # period of shock
ks[0:ts] = k[0:ts]


def euler_new(k0,k1,k2):
    return pow(k1,1-theta)*pow(z1*h,theta)+(1-delta)*k1-k2               \
    -beta*((1-theta)*pow(k1,-theta)*pow(z1*h,theta)+(1-delta))      \
    *(pow(k0,1-theta)*pow(z1*h,theta)-k1+(1-delta)*k0)

def ktrans_new(z):
    F = np.zeros(size-ts)
    z = z
    F[0] = euler_new(ks[ts-1],z[1],z[2])
    z[size-1-ts] = kss1
    F[size-2-ts] = euler_new(z[size-3-ts], z[size-2-ts], z[size-1-ts])
    for i in range(1,size-2-ts):
        F[i] = euler_new(z[i],z[i+1],z[i+2])
    return F
z = np.ones(size-ts)*8 # The initially assumed steady state
k = fsolve(ktrans_new, z)
k[0] = ks[ts-1]
ks[ts:size] = k

plt.plot(ks)
plt.title('Transition path of capital with technology shock at t=10')
plt.xlabel('Time')
plt.ylabel('Capital')
plt.show()

savingspath_new = np.empty(size)
# We lose one degree of freedom due to the t+1 structure of investment (=savings)
for i in range(size-1):
        savingspath_new[i] = ks[i+1]-(1-delta)*ks[i]
        savingspath_new[size-1] = savingspath_new[size-2]

z_new1 = np.ones(ts)*z2
z_new2 = np.ones(size-ts)*z1
znew = np.concatenate((z_new1,z_new2))

outputpath_new = pow(ks,1-theta)*pow(znew*h, theta)
consumptionpath_new = outputpath_new - savingspath_new

cons, = plt.plot(consumptionpath_new, label='Consumption')
sav, = plt.plot(savingspath_new, label='Savings')
outp,= plt.plot(outputpath_new,label='Output')
plt.title('Different transition paths with shock at t=10')
plt.legend(handles = [cons, sav, outp], loc='upper right')
plt.xlabel('Time')
plt.ylabel('Levels')
plt.show()
