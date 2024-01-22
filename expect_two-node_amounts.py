from mpmath import quad, mp
import numpy as np
import math
import time
from scipy.special import beta
import pandas as pd

mp.dps = 15

def min_v_0(t,f):
    func = - 0.4 - f + pow( 7-t, -1)
    return func

def min_v_1(t):
    kappa_ = 0.1
    func = - 0.2 - kappa_*pow(7-t, -1) 
    return func

def integrand(t,f):
    integrand = - min_v_1(t)+min_v_0(t,f)
    return integrand

def exp_min_v_0(f, t):
    int = quad(lambda s: min_v_0(s,f), [0,t])
    return math.exp(int)

def exp_min_v_1(t):
    int = quad(lambda s: min_v_1(s), [0,t])
    return math.exp(int)

def exp_1_0_integrand(r,f):
    int = quad(lambda s: integrand(s,f), [0,r] )
    return math.exp(int)*f

def integral_x_1(t,f):
    int = quad(lambda r: exp_1_0_integrand(r,f), [0,t])
    return int

def x_0(x_0_0, f, t):
    x_0 = x_0_0*exp_min_v_0(f,t)
    return x_0

def x_1(x_0_0, x_1_0, f, t):
    x_1 = exp_min_v_1(t)*(x_1_0+x_0_0*integral_x_1(t,f))
    return x_1

def expect_x_0(x_0_0, a, b, t):
    int = quad( lambda f: x_0(x_0_0, f, t)*pow(f, a-1)*pow(1-f, b-1), [0,1])
    return pow(beta(a, b, out=None),-1)*int
 
def expect_x_1(x_0_0, x_1_0, a, b, t):
    int = quad( lambda f: x_1(x_0_0, x_1_0, f, t)*pow(f, a-1)*pow(1-f, b-1), [0,1])
    return pow(beta(a, b, out=None),-1)*int 


a = 1
b = 3
x_0_0 = 1
x_1_0 = 1.5
t_max = 5
t_step = 0.25

list_exp_0 = []
list_exp_1 = []
list_t = []

t = 0
while t<=t_max+0.5*t_step:
    list_exp_0.append(expect_x_0(x_0_0, a, b, t))
    list_exp_1.append(expect_x_1(x_0_0, x_1_0, a, b, t))
    list_t.append(t)

    t += t_step

print(list_exp_0)    
print(list_exp_1)

dict = {'t': list_t, 'E[x_0(t)]': list_exp_0, 'E[x_1(t)]': list_exp_1}
   
df = pd.DataFrame(dict) 
    

df.to_csv('expected_two-node_kappa01.csv') 