from mpmath import quad
import numpy as np
import math
import time
from scipy.special import beta
import pandas as pd
from operator import add

def f_(a,b):
   arg = np.random.beta(a,b)
   return arg 

def g(n,r):
    if n==r:
        g = 0.1*(n)*(n) + 0.1
    else:
        g = 0    
    return g 


a=[1.2,1.5,1.1,1.7,4.1,2.1,3.1,1.6,0.82]
b=[0.4,1.2,3,1,2,2.2,1.8,0.7,0.1]  #[1,1,1,1,1,1,1,1,1]
delta_list = [0.2, 0.3, 0.2, 0.4, 0.1, 0.7, 0.2, 0.3, 0.5] #  [2,3,2,4,3,5,2.3,1.9,3.1]#


N =10
r_ = 3
r_ = r_-1
exp_x_n_r = [0]*N
exp_distr_x_n_r = [0]*N
realiz = 50000
r = 1
while r<= realiz:
    x_n = 1
    n=1
    n_list = [1]
    x_n_list = [x_n]
    total_amount = x_n
    f_n_min_1 = f_(a[0],b[0])
    while n<N:

        if n<N-1:
            f_n = f_(a[n],b[n])
            x_n = (f_n_min_1+g(n-1,r_))*x_n/(f_n+delta_list[n-1]+g(n,r_))
            f_n_min_1 = f_n 
        if n==N-1:
            x_n = (f_n_min_1+g(n-1,r_))*x_n/(delta_list[n-1]+g(n,r_))   
   
        n_list.append(n+1)
        x_n_list.append(x_n)
        total_amount += x_n
        n += 1

    exp_x_n_r = list(map(add, exp_x_n_r, x_n_list))
    newList = [x/ total_amount for x in x_n_list] 
    exp_distr_x_n_r = list(map(add, exp_distr_x_n_r,newList))  

    r += 1

exp_x_n = [s/realiz for s in exp_x_n_r]    
exp_distr_x_n = [s/realiz for s in exp_distr_x_n_r]    

print(exp_x_n)
print(n_list)    
print(exp_distr_x_n)    
dict = {'n': n_list, 'E[p_n(t)]': exp_distr_x_n}
   
df = pd.DataFrame(dict) 
    

df.to_csv('beta_distr_stationary_stim_1_at3_distr_mean.csv') 

