from mpmath import quad, mp
import numpy as np
import math
import time
import pandas as pd

mp.dps = 15

def f_(a,b):
   arg = np.random.beta(a,b)
   return arg 

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


a=1.4
b=1.2

a = 1
b = 3
realiz = 1000
x_0_0 = 1
x_1_0 = 1.5
t_max = 5
t_step = 0.25

list_exp_0 = []
list_exp_1 = []
list_exp_tot = []
list_distr_0 = []
list_distr_1 = []
list_t = []

t=0
while t<=t_max+0.5*t_step:
    k=1
    list_x_0 =[]
    list_x_1 =[]
    total_amount=[]
    distr_0 = []
    distr_1 = []
    while k<=realiz:
        f = f_(a,b)
        x_0_t = x_0(x_0_0, f, t)
        list_x_0.append(x_0_t)
        x_1_t = x_1(x_0_0, x_1_0, f, t)
        list_x_1.append(x_1_t)
        tot = x_0_t + x_1_t
        total_amount.append(tot)
        distr_0.append(x_0_t/tot)
        distr_1.append(x_1_t/tot)    
        k+= 1

    list_exp_0.append(sum(list_x_0)/len(list_x_0))    
    list_exp_1.append(sum(list_x_1)/len(list_x_1)) 
    list_exp_tot.append(sum(total_amount)/len(total_amount)) 
    list_distr_0.append(sum(distr_0)/len(distr_0))
    list_distr_1.append(sum(distr_1)/len(distr_1))

    list_t.append(t)
    t += t_step

print(list_exp_0)    
print(list_exp_1)  
print(list_t)

dict = {'t': list_t, 'E[x_0(t)]': list_exp_0, 'E[x_1(t)]': list_exp_1, 'E[distr0]':list_distr_0, 'E[distr1]':list_distr_1}
   
df = pd.DataFrame(dict) 
    

df.to_csv('multip_realiz1000_two-node_kappa01.csv') 