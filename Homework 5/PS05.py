import numpy as np
import sympy as sp
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import seaborn as sns; sns.set(color_codes=True)
import time
import random
import numba
from numba import jit

### Exercise 1.1: Simulating data and plotting PDFs

### Since the logged variables are jointly normally distributed, we need the covariance matrix (identitiy) and expectations
### to fully describe the distribution. Derivation is done in a separate script called moments.py

np.random.seed(13)
#cov = np.identity(2)
#cov = np.array([[1,0.5],[0.5,1]])
cov = np.array([[1,-0.5],[-0.5,1]])
mu = np.array([-0.5, -1.25])
size = 1000000 # at 10 million my RAM is overloaded

### If a vector X is normally distributed, then exp(X) is lognormally distributed with the same mean and variance

log_data = np.random.multivariate_normal(mu,cov, size=size)
level_data = np.exp(log_data)
k = level_data[:,1]
z = level_data[:,0]
lnk = log_data[:,1]
lnz = log_data[:,0]


### Plotting the joint density functions for levels and for logs
## First levels
sns.jointplot(k,z,kind="hex").set_axis_labels("Capital", "Productivity")
plt.show()

sns.jointplot(lnk,lnz,kind="hex").set_axis_labels("Log Capital", "Log Productivity")
plt.show()
'''### Plotting the raw joint density of lognormal variables does not make much sense as in 10,000,000 observations there will be massive outliers
### I atempt to get rid of these outliers for plotting purposes

meank = np.mean(k)
sdk = np.std(k)
final_k = [x for x in k if (x > meank - 2 * sdk)]
final_k = [x for x in final_k if (x < meank + 2 * sdk)]


meanz = np.mean(z)
sdz = np.std(z)
final_z = [x for x in z if (x > meanz - 2 * sdz)]
final_z = [x for x in final_z if (x < meanz + 2 * sdz)]

leng = min(len(final_k), len(final_z))

### Of course now the two vectors will have different lengths so I will just arbitrarily pick the first observations that match from both for plotting
knew = np.asarray(final_k[0:leng])
znew = np.asarray(final_z[0:leng])
### Now try to plot the joint density for levels again
sns.jointplot(k,z,kind="hex").set_axis_labels("Capital", "Productivity")
plt.show()
'''
### Exercise 1.2: Compute firm output
#gamma = 0.6
gamma = 0.8 ### For Exercise 2
s = pow(z,1/(1-gamma))
y = s*pow(k,gamma)


### Exercise 1.3: Solve the maximization problem
K = np.sum(k)


def k0(k0):
    f1 = 0
    for i in range(len(s)):
        f1 += k0*s[i]/s[0]
    f2 = K-f1
    return(f2)

kstart = fsolve(k0,1)
k_e = np.empty(size)
for i in range(size):
    k_e[i] = kstart*s[i]/s[0]

### Exercise 1.4: Compare the efficient allocation to the 'data'
# First the efficient one
sns.regplot(x=s, y=k_e, color="g")
plt.show()

# Now the supposedly inefficient one
sns.regplot(x=s, y=k, color="g")
plt.show()


### Exercise 1.5: Compute the gains from reallocation
def Y(k):
    return sum(pow(s,1-gamma)*pow(k,gamma))

gains1 = 100*(Y(k_e)/Y(k)-1)
print('The gains from reallocation could give an output which is higher by a factor of', round(gains1,2),'% for gamma =',gamma,'.')


### Exercise 3 ###

### Exercise 3.1: Random sampling

nsample = 1000 # Number of sample points for our random sample
#nsample = 100
#nsample = 1000
#nsample = 1000000

lnk_pop = list(lnk)
lnz_pop = list(lnz)
s_pop = list(s)
lk_sam = random.sample(lnk_pop,nsample)
lz_sam = random.sample(lnz_pop,nsample)


### Calculating moments from the random sample for comparison
var_sam_k = np.var(lk_sam)
var_pop_k = np.var(lnk)

print('The variance of log capital in the random sample is', round(var_sam_k,2), 'compared to', round(var_pop_k,2), 'in the data.')

var_sam_z = np.var(lz_sam)
var_pop_z = np.var(lnz)

print('The variance of log capital in the random sample is', round(var_sam_z,2), 'compared to', round(var_pop_z,2), 'in the data.')

corr_zk_sam = np.corrcoef(lk_sam, lz_sam)[1,0]
corr_zk_pop = np.corrcoef(lnk,lnz)[1,0]

print('The Pearson correlation coefficient of log capital and log productivity in the random sample is', round(corr_zk_sam,2), 'compared to', round(corr_zk_pop,2), 'in the data.')

### Exercise 1.2: Comparing the reallocation gains between sample and population
## Making level data
ksam = np.exp(lk_sam)
zsam = np.exp(lz_sam)
ssam = pow(zsam,1/(1-gamma))
ysam = ssam*pow(ksam,gamma)

## Calculating the optimal k_e in the sample
Ksam = np.sum(ksam)


def k0(k0):
    f1 = 0
    for i in range(len(ssam)):
        f1 += k0*ssam[i]/ssam[0]
    f2 = Ksam-f1
    return(f2)

kstartsam = fsolve(k0,1)
k_e_sam = np.empty(nsample)
for i in range(nsample):
    k_e_sam[i] = kstartsam*ssam[i]/ssam[0]

# Plotting the sampled observations
# First the efficient one
sns.regplot(x=ssam, y=k_e_sam, color="g")
plt.show()

# Now the supposedly inefficient one
sns.regplot(x=ssam, y=ksam, color="g")
plt.show()


def Y(k):
    return sum(pow(ssam,1-gamma)*pow(k,gamma))

gains_sam = 100*(Y(k_e_sam)/Y(ksam)-1)
print('For the sample, the gains from reallocation could give an output which is higher by a factor of', round(gains_sam,2),'% for gamma =',gamma,'which compares to', round(gains1,2), 'in the population.')


### Repeating 3.1 and 3.2 1000 times:
reps = 1000
gain_sim = np.empty(reps)
for i in range(reps):
    lk_sam = random.sample(lnk_pop,nsample)
    lz_sam = random.sample(lnz_pop,nsample)

    ksam = np.exp(lk_sam)
    zsam = np.exp(lz_sam)
    ssam = pow(zsam,1/(1-gamma))

    Ksam = np.sum(ksam)
    def k0(k0):
        f1 = 0
        for i in range(len(ssam)):
            f1 += k0*ssam[i]/ssam[0]
        f2 = Ksam-f1
        return(f2)

    kstartsam = fsolve(k0,1)
    k_e_sam = np.empty(nsample)
    for i in range(nsample):
        k_e_sam[i] = kstartsam*ssam[i]/ssam[0]
    gains_sam = 100*(Y(k_e_sam)/Y(ksam)-1)

    gain_sim[i] = gains_sam

plt.hist(gain_sim, bins = 100)
plt.show()

### Calculate probability that sample gain is within 10% of actual misallocation
actual = gains1
interval = 0.1
ub = actual*(1+interval)
lb = actual*(1-interval)

hits = []
for i in range(reps):
    if (gain_sim[i]<ub) and (gain_sim[i]>lb):
        hit = gain_sim[i]
        hits.append(hit)

tot = len(hits)
prob = round(tot/reps*100,3)

print('The probability of finding a sample gain within a', int(interval*100), '% interval around the population gain is given by', prob, '%.')
