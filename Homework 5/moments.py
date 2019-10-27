import numpy as np
import random
from scipy.optimize import fsolve
import math
import scipy.integrate as integrate
import scipy.special as special

gamma = 0.6
np.random.seed(13)

### Checking the mean of ln(k) which I calculate as -0.5. Does this give E[k] =1?
lnk = np.random.normal(-0.5,1,10000000)
ek = np.exp(lnk)
q = np.mean(ek)
print(q) # If this is equal or close to 1 we are good



### Finding the mean of ln(z). Does this give E[s] = 1?
def integrand(s,mu):
    return (1-gamma)/(np.sqrt(2*math.pi))*np.exp(-0.5*((1-gamma)*np.log(s)-mu)**2)

def expints(mu):
    return integrate.quad(integrand, 0, np.inf, args=(mu))[0]

def int(mu):
    f1 = 1-expints(mu)
    return(f1)

eq =fsolve(int,0.5)
print(eq) # This should be the mean that gives us E[s] = 1

test = expints(eq)
print(test) # Evaluating the integral for this mu should give us 1 as we have log-normal density


lnz = np.random.normal(eq,1,1000000)
z = np.exp(lnz)
s = pow(z, 1/(1-gamma))
m = np.mean(s)
print(m) # This should be close to 1
