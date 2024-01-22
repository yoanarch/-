from mpmath import quad
import numpy as np
import math
import time
from scipy.special import beta
import pandas as pd


# интеграл, съответстващ на множителя, свързан с номера на разглежданата клетка

def g(n,r):
    if n==r:
        g = 0.1*(n)*(n) + 0.1
    else:
        g = 0    
    return g 
def integrand_n(x,delta_,a,b,n,r):
   return pow((x),a[n]-1)*pow((1-x),b[n]-1)*pow((x+delta_[n-1]+g(n,r)),-1)

def integrand_less_than_n(x,delta_,a,b,k,r):
    return (x+g(k,r))*pow(x,a[k]-1)*pow((1-x),b[k]-1)*pow((x+delta_[k-1]+g(k,r)),-1)

def multiplier_n(delta_,a, b,n,r):
    int_n = quad(lambda x: integrand_n(x,delta_,a,b,n,r), [0, 1]) 
    return pow(beta(a[n], b[n], out=None),-1)*int_n

def multiplier_less_than_n(delta_,a,b,k,r):
    int_k = quad(lambda x: integrand_less_than_n(x,delta_,a,b,k,r), [0, 1]) 
    return pow(beta(a[k], b[k], out=None),-1)*int_k


start_time = time.time()

l = 1
N = 10

a=[1.2,1.5,1.1,1.7,4.1,2.1,3.1,1.6,0.82]
b=[0.4,1.2,3,1,2,2.2,1.8,0.7,0.1]  #[1,1,1,1,1,1,1,1,1]
delta_list = [0.2, 0.3, 0.2, 0.4, 0.1, 0.7, 0.2, 0.3, 0.5] #  [2,3,2,4,3,5,2.3,1.9,3.1]#
r = 2
r = r-1
x_l_0 =1

mean_x_in_n = []
n_list = []
analytical = []

prod = a[0]*x_l_0/(a[0]+b[0])
mean_x_in_n.append(x_l_0)
analytical.append(prod)
n_list.append(l)

# пресмятане на E[x_{l+1}]
last_ = multiplier_n(delta_list, a, b, 1,r)
prod_ = last_*prod
mean_x_in_n.append(prod_)
n_list.append(l+1)

k = l+2
while k<N:
    s = k-l

    prod_ = prod_/ last_
    last_ = multiplier_n(delta_list, a, b, s,r)
    
    prod_ = prod_*multiplier_less_than_n(delta_list, a, b, s-1,r)*last_
    # аналитичен резултат
    #prod = prod * multiplier_less_than_n(a,b,k)

    mean_x_in_n.append(prod_)
    n_list.append(k)
    #analytical.append(((a[k]+b[k]-1)/(a[k]-1))*(a[0]/(a[0]+b[0])))
    k += 1

# пресмятане на E[x_{N}]    
prod_ = prod_/ last_
prod_ = prod_*multiplier_less_than_n(delta_list, a, b, len(a)-2,r)*pow(delta_list[len(delta_list)-1],-1)
mean_x_in_n.append(prod_)
n_list.append(N)

print(n_list)
print(mean_x_in_n)
#print(analytical)

dict = {'n': n_list, 'E[x_n(t)]': mean_x_in_n}
   
df = pd.DataFrame(dict) 
    

#df.to_csv('beta_distr_stationary_stimulus2_at_3.csv') 
# print(t_list)


print("--- %s seconds ---" % (time.time() - start_time))