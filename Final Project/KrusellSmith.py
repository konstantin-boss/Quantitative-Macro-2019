import numpy as np
from scipy import optimize
from scipy.integrate import fixed_quad
import matplotlib.pyplot as plt
import quantecon
import statsmodels.api as sm

np.random.seed(13)

## Exercise 1: Simple Variant of Krusell-Smith Algorithm ##
## Part 1.2 ##

# Initializing parameters
alpha = 0.3
t = 40
beta = 0.99**t
BETA = beta/(1+beta)
tau = 0
lam = 0.5
theta = 0.5
g = 0
T = 50000
burn = 500
ups = np.empty(T)
ups[0] = 1
for i in range(1,T):
    ups[i] = (1+g)*ups[i-1]

sd_zeta = np.sqrt(40*0.02**2)
sd_rho = np.sqrt(40*0.08**2)
sd_eta = np.sqrt(40*0.15**2)
mu_zeta = 1
mu_rho = 1
mu_eta = 1
zeta = np.exp(np.random.normal(mu_zeta, sd_zeta,2))
rho = np.exp(np.random.normal(mu_rho, sd_rho,2))
eta = quantecon.quad.qnwnorm(11, mu=mu_eta, sig2=sd_eta**2)
zeta_sam = np.random.choice(zeta,T)
rho_sam = np.random.choice(rho,T)
zlow = np.array([zeta[0], rho[0]])
zhigh = np.array([zeta[1], rho[1]])
zett = np.array([zlow,zhigh])


rho_prob = np.array([0.5, 0.5])
eta_prob = eta[1]
phi = 0
for ir, r in enumerate(rho):
    for ie, e in enumerate(np.exp(eta[0])):
        phi += rho_prob[ir]*eta_prob[ie]*(1+(1-alpha)/(alpha*(1+lam)*r)*(lam*e+tau*(1+lam*(1-e))))**(-1)
s = beta * phi/(1+beta*phi)


k = np.empty(T)
k[0] = np.exp(1/(1-alpha)*(np.log(s)+np.log(1-tau)+np.log(1-alpha)))

for i in range(1,T):
    k[i] = np.exp(np.log(s)+np.log(1-tau)+np.log(1-alpha)+np.log(zeta_sam[i-1])+alpha*np.log(k[i-1]))

# Exercise 1.3
# Part a)
# Calculating Psi
PSI0 = np.empty([len(zett),len(zett)])
for i in range(len(zett)):
    PSI0[i,0] = np.log(s)+np.log(1-tau)+np.log(1-alpha)+np.log(zett[i,0])
    PSI0[i,1] = alpha
print(PSI0)

# Calculating kstar
kstar = np.exp((np.log(s)+np.log(1-tau)+np.log(1-alpha))/(1-alpha))

# Part b)
def solveHH(PSI):
    kgrid = np.linspace(0.5*kstar,1.5*kstar,5)
    k1 = np.empty([len(PSI),len(kgrid)])
    for i,k in enumerate(kgrid):
        for j in range(len(PSI)):
            k1[j,i] = np.exp(PSI[j,0]+PSI[j,1]*np.log(k))
    savgrid = np.empty([len(zett),len(kgrid)])
    for i in range(len(zett)):
        for j,k in enumerate(kgrid):
            def EE(a):
                ee = 1-beta*((1-tau)*(1-alpha)*k**alpha*zett[i,0]-a)*(1+alpha*k1[i,j]**(alpha-1)*zett[i,0]*zett[i,1])/(a*(1+alpha*k1[i,j]**(alpha-1)*zett[i,0]*zett[i,1])+lam*(1-alpha)*k1[i,j]**alpha*zett[i,0]*(1-tau)+(1-lam)*tau*(1-alpha)*k1[i,j]**alpha*zett[i,0]*(1+lam)/(1-lam))
                return(ee)
            savgrid[i,j] = optimize.brentq(EE,0.0001,1)/((1-tau)*(1-alpha)*k**alpha*zett[i,0])
    return(savgrid)


def simModel(savgrid):
    s1 = np.random.choice(savgrid[0],T)
    s2 = np.random.choice(savgrid[1],T)
    sav = np.array([s1,s2])
    ksim = np.empty([len(zett),T])
    ksim[:,0] = kstar
    c1 = np.empty([len(zett),T])
    c2 = np.empty([len(zett),T])
    u = np.empty(len(zett))
    for i in range(len(zett)):
        for j in range(1,T):
            ksim[i,j] = np.exp(np.log(sav[i,j])+np.log(1-tau)+np.log(1-alpha)+np.log(zeta_sam[j])+alpha*np.log(ksim[i,j-1]))
            c1[i,j] = (1-sav[i,j])*(1-tau)*(1-alpha)*zeta_sam[j]*ksim[i,j]**alpha
            c2[i,j] = beta*c1[i,j]*(1+alpha*ksim[i,j]**(alpha-1)*zeta_sam[j]*rho_sam[j])
    return(ksim[0], ksim[1], c1,c2)

def regCoef(kl,kh): # Has to return 2x2 matrix
    lnk_l = np.log(kl)[0:T-1]
    lnk_h = np.log(kh)[0:T-1]
    lnk1_h = np.log(kh)[1:T]
    lnk1_l = np.log(kl)[1:T]

    Xh = lnk_h[burn:]
    Xh = sm.add_constant(Xh)
    Xl = lnk_l[burn:]
    Xl = sm.add_constant(Xh)

    yh = lnk1_h[burn:]
    yl = lnk1_l[burn:]

    model_l = sm.OLS(yl, Xl)
    results_l = model_l.fit()
    model_h = sm.OLS(yh, Xh)
    results_h = model_h.fit()
    psi_l = np.array(results_l.params)
    psi_h = np.array(results_h.params)
    PSI = np.array([psi_l, psi_h])
    return(PSI)



count = 0
omega = 0.9
PSIold = PSI0
PSInew = np.array([[3,2],[1,2]])
epsilon = 0.0001

while np.max(np.abs(PSInew-PSIold))>epsilon:
    count += 1
    if count == 1:
        savgrid = solveHH(PSI0)
        sim = simModel(savgrid)
        c1 = sim[2]
        c1l = c1[0]
        c1l = c1l[burn:]
        c1h = c1[1]
        c1h = c1h[burn:]
        c2 = sim[3]
        c2l = c2[0]
        c2l = c2l[burn:]
        c2h = c2[1]
        c2h = c2h[burn:]
        utill = 1/(T-burn)*np.sum((1-BETA)*np.log(c1l)+BETA*np.log(c2l))
        utilh = 1/(T-burn)*np.sum((1-BETA)*np.log(c1h)+BETA*np.log(c2h))
        PSIint = regCoef(sim[0],sim[1])
        PSInew = omega*PSIint+(1-omega)*PSI0
    else:
        savgrid = solveHH(PSInew)
        sim = simModel(savgrid)
        c1 = sim[2]
        c1l = c1[0]
        c1l = c1l[burn:]
        c1h = c1[1]
        c1h = c1h[burn:]
        c2 = sim[3]
        c2l = c2[0]
        c2l = c2l[burn:]
        c2h = c2[1]
        c2h = c2h[burn:]
        utill = 1/(T-burn)*np.sum((1-BETA)*np.log(c1l)+BETA*np.log(c2l))
        utilh = 1/(T-burn)*np.sum((1-BETA)*np.log(c1h)+BETA*np.log(c2h))
        PSIint = regCoef(sim[0],sim[1])
        PSIold = PSInew
        PSInew = omega*PSIint+(1-omega)*PSIold
    print(count)

print(PSInew)
print(savgrid)
print(utill)
print(utilh)
