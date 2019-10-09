import numpy as np
import sympy as sp
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import math
import intervals as I
## Throughout this exercise I will assume that a countr has a union mass of
## households/agents with shares of high and low productivity which can be set
## I will also assume that k^bar is the given capital stock in a country which
## consists of the (given) capital holdings of the low and high types

## Setting parameters
pi_a        = 0.5 # The proportion of workers with a high productivity in country A
pi_b        = 0.5
nu          = 1
kappa       = 5
sigma       = 0.8
eta_h_a     = 5.5
eta_l_a     = 0.5
eta_h_b     = 5.5
eta_l_b     = 2.5
Z           = 1
theta       = 0.6
k_bar       = 2
lambda_a    = 0.95
lambda_b    = 0.84
phi         = 0.2
k_l_a       = 1 # I assign an initial level of capital to each type in a country
k_h_a       = 1
k_l_b       = 1
k_h_b       = 1

## Defining the model's fucntion for checking market clearance and utility maximization later
def production(K,H):
    return Z*pow(K,1-theta)*pow(H,theta)

def utility(c,h):
    return pow(c,1-sigma)/(1-sigma)-kappa*pow(h,1+1/nu)/(1+1/nu)

## Euler equations, budget constraints and firm FOCs for both countries and types
## The h subscript refers to a guy with high prodcutivity, l to low productivity
# h_l = x[0]
# h_h = x[1]
# c_l = x[2]
# c_h = x[3]
# w = x[4]
# r = x[5]

def f(x):
    f1 = -kappa*pow(x[0],1/nu)+pow(x[2],-sigma)*x[4]*eta_l_a*lambda_a*(1-phi)*pow(x[4]*x[0]*eta_l_a,-phi)       # Euler low types
    f2 = -kappa*pow(x[1],1/nu)+pow(x[3],-sigma)*x[4]*eta_h_a*lambda_a*(1-phi)*pow(x[4]*x[1]*eta_h_a,-phi)       # Euler high types
    f3 = -x[5]+Z*(1-theta)*pow(k_bar,-theta)*pow(x[0]*eta_l_a + x[1]*eta_h_a, theta)  # Firm FOC rate of return
    f4 = -x[4]+Z*(theta)*pow(k_bar,1-theta)*pow(x[0]*eta_l_a + x[1]*eta_h_a, theta-1) # Firm FOC wage
    f5 = -x[2] + lambda_a*pow(x[4]*x[0]*eta_l_a, 1-phi)+x[5]*pow(k_l_a,eta_l_a)                         # BC low types
    f6 = -x[3] + lambda_a*pow(x[4]*x[1]*eta_h_a, 1-phi)+x[5]*pow(k_h_a,eta_h_a)                         # BC high types
    return[f1,f2,f3,f4,f5,f6]

equilibrium_a = fsolve(f, [1,1,1,1,1,1])

def g(x):
    g1 = -kappa*pow(x[0],1/nu)+pow(x[2],-sigma)*x[4]*eta_l_b*lambda_b*(1-phi)*pow(x[4]*x[0]*eta_l_b,-phi)       # Euler low types
    g2 = -kappa*pow(x[1],1/nu)+pow(x[3],-sigma)*x[4]*eta_h_b*lambda_b*(1-phi)*pow(x[4]*x[1]*eta_h_b,-phi)       # Euler high types
    g3 = -x[5]+Z*(1-theta)*pow(k_bar,-theta)*pow(x[0]*eta_l_b + x[1]*eta_h_b, theta)      # Firm FOC rate of return
    g4 = -x[4]+Z*(theta)*pow(k_bar,1-theta)*pow(x[0]*eta_l_b + x[1]*eta_h_b, theta-1)     # Firm FOC wage
    g5 = -x[2] + lambda_b*pow(x[4]*x[0]*eta_l_b, 1-phi)+x[5]*pow(k_l_b,eta_l_b)                         # BC low types
    g6 = -x[3] + lambda_b*pow(x[4]*x[1]*eta_h_b, 1-phi)+x[5]*pow(k_h_b,eta_h_b)                         # BC high types
    return[g1,g2,g3,g4,g5,g6]

equilibrium_b = fsolve(g, [1,1,1,1,1,1])

print('A: The equilibrium consumption of low types is', round(equilibrium_a[2],2))
print('A: The equilibrium consumption of high types is', round(equilibrium_a[3],2))
print('A: The equilibrium labor supply of low types is', round(equilibrium_a[0],2))
print('A: The equilibrium labor supply of high types is', round(equilibrium_a[1],2))
print('A: The equilibrium rate of return is', round(equilibrium_a[5],2))
print('A: The equilibrium wage is', round(equilibrium_a[4],2))
print('B: The equilibrium consumption of low types is', round(equilibrium_b[2],2))
print('B: The equilibrium consumption of high types is', round(equilibrium_b[3],2))
print('B: The equilibrium labor supply of low types is', round(equilibrium_b[0],2))
print('B: The equilibrium labor supply of high types is', round(equilibrium_b[1],2))
print('B: The equilibrium rate of return is', round(equilibrium_b[5],2))
print('B: The equilibrium wage is', round(equilibrium_b[4],2))

## Part 2: A capital union - need to solve 16 unknowns
# h_l_a = x[0]      h_l_b = x[8]
# h_h_a = x[1]      h_h_b = x[9]
# c_l_a = x[2]      c_l_b = x[10]
# c_h_a = x[3]      c_h_b = x[11]
# k_ls_a = x[4]     k_ls_b = x[12]
# k_hs_a = x[5]     k_hs_b = x[13]
# w_a = x[6]        w_b = x[14]
# r_a = x[7]        r_b = x[15]

def h(x):
    h1 = -kappa*pow(x[0],1/nu)+pow(x[2],-sigma)*x[6]*eta_l_a*lambda_a*(1-phi)*pow(x[6]*x[0]*eta_l_a,-phi)       # Euler low types country A
    h2 = -kappa*pow(x[1],1/nu)+pow(x[3],-sigma)*x[6]*eta_h_a*lambda_a*(1-phi)*pow(x[6]*x[1]*eta_h_a,-phi)       # Euler high types country A
    h3 = -x[7]+Z*(1-theta)*pow(x[4]+x[5]+(k_l_b-x[12])+(k_h_b-x[13]),-theta)*pow(x[0]*eta_l_a + x[1]*eta_h_a, theta)     # Firm FOC rate of return A
    h4 = -x[6]+Z*(theta)*pow(x[4]+x[5]+(k_l_b-x[12])+(k_h_b-x[13]),1-theta)*pow(x[0]*eta_l_a + x[1]*eta_h_a, theta-1)    # Firm FOC wage A
    h5 = -x[15]+x[7]*eta_l_a*pow(x[4],eta_l_a-1)                                                                           # Euler for low type capital A
    h6 = -x[15]+x[7]*eta_h_a*pow(x[5],eta_h_a-1)                                                                           # Euler for high type capital A
    h7 = -x[2] + lambda_a*pow(x[6]*x[0]*eta_l_a, 1-phi)+x[7]*pow(x[4],eta_l_a)+x[15]*(k_l_a-x[4])                                # BC low types A
    h8 = -x[3] + lambda_a*pow(x[6]*x[1]*eta_h_a, 1-phi)+x[7]*pow(x[5],eta_h_a)+x[15] *(k_h_a-x[5])                              # BC high types A

    h9 = -kappa*pow(x[8],1/nu)+pow(x[10],-sigma)*x[14]*eta_l_b*lambda_b*(1-phi)*pow(x[14]*x[8]*eta_l_b,-phi)       # Euler low types country A
    h10 = -kappa*pow(x[9],1/nu)+pow(x[11],-sigma)*x[14]*eta_h_b*lambda_b*(1-phi)*pow(x[14]*x[9]*eta_h_b,-phi)       # Euler high types country A
    h11 = -x[15]+Z*(1-theta)*pow(x[12]+x[13]+(k_l_a-x[4])+(k_h_a-x[5]),-theta)*pow(x[8]*eta_l_b + x[9]*eta_h_b, theta)     # Firm FOC rate of return A
    h12 = -x[14]+Z*(theta)*pow(x[12]+x[13]+(k_l_a-x[4])+(k_h_a-x[5]),1-theta)*pow(x[8]*eta_l_b + x[9]*eta_h_b, theta-1)    # Firm FOC wage A
    h13 = -x[7]+x[15]*eta_l_b*pow(x[12],eta_l_b-1)                                                                           # Euler for low type capital A
    h14 = -x[7]+x[15]*eta_h_b*pow(x[13],eta_h_b-1)                                                                           # Euler for high type capital A
    h15 = -x[10] + lambda_b*pow(x[14]*x[8]*eta_l_b, 1-phi)+x[15]*pow(x[12],eta_l_b)+x[7]*(k_l_b-x[12])                                # BC low types A
    h16 = -x[11] + lambda_b*pow(x[14]*x[9]*eta_h_b, 1-phi)+x[15]*pow(x[13],eta_h_b)+x[7] *(k_h_b-x[13])                         # BC high types



    return[h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13,h14,h15,h16]

equilibrium_union = fsolve(h, [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
print(' ')
print('########## RESULTS FOR THE UNION ##########')
print('Union: The eq. consumption of low types in A is', round(equilibrium_union[2],2))
print('Union: The eq. consumption of high types in A is', round(equilibrium_union[3],2))
print('Union: The eq. labor supply of low types in A is', round(equilibrium_union[0],2))
print('Union: The eq. labor supply of high types in A is', round(equilibrium_union[1],2))
print('Union: The eq. domestic supply of capital of low types in A is', round(equilibrium_union[4],2))
print('Union: The eq. domestic supply of capital of high types in A is', round(equilibrium_union[5],2))
print('Union: The eq. rate of return in A is', round(equilibrium_union[7],2))
print('Union: The eq. wage in A is', round(equilibrium_union[6],2))
print('Union: The eq. consumption of low types in B is', round(equilibrium_union[10],2))
print('Union: The eq. consumption of high types in B is', round(equilibrium_union[11],2))
print('Union: The eq. labor supply of low types in B is', round(equilibrium_union[8],2))
print('Union: The eq. labor supply of high types in B is', round(equilibrium_union[9],2))
print('Union: The eq. domestic supply of capital of low types in B is', round(equilibrium_union[12],2))
print('Union: The eq. domestic supply of capital of high types in B is', round(equilibrium_union[13],2))
print('Union: The eq. rate of return in B is', round(equilibrium_union[15],2))
print('Union: The eq. wage in B is', round(equilibrium_union[14],2))
