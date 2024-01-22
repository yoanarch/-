from mpmath import quad, mp
import numpy as np
import math
import time
import pandas as pd

# def minus_v_n(x,n):
#     # sigma_n_t = (1/(2-n))*np.cos(0.5*np.pi*x)
#     # mu_n_t = (1/(2-n))*np.cos((2/3)*np.pi*x)
#     # epsilon_n_t = (n+1)*np.log((x+1)*2)
#     # gamma_n_t = (1/(2-n))*np.log(x+1)*0.3
#     # f_n_t = (1/(2-n))*0.3*(x+1)
    
#     # arg = sigma_n_t-mu_n_t+epsilon_n_t-gamma_n_t-f_n_t
#     # arg = np.cos(np.pi*x) +n*3.1
#     arg = (n+1)#*x**2
    
#     return arg

def minus_v_n(x,n):
    
    if n == 0:
        arg = -0.7* math.sin(math.pi*x) - 0.2 + 0.3

    if n == 1:
        arg = - 0.3* math.sin(math.pi*x) - 0.25 - 0.3

    if n == 2:
        arg = - 0.35* math.sin(math.pi*x) - 0.3    

    return arg

# def minus_v_n(x,n):
    
#     if n == 0:
#         arg = - 0.9 + pow(7-x, -1)

#     kappa_ = -1.2
#     if n == 1:
#         arg = - 0.2 - kappa_ * pow(7-x, -1)  

#     lambda_ = -1.4
#     if n == 2:
#         arg = - 0.1 - lambda_ * pow(7-x, -1)  
#     return arg        

   

# пресмятане на \exp( \int_0^s {(v_1(k) - v_0(k))} dk) * f_0(s)    
def integrand_1_0(s):
    n = 0 
    in_int_0 = quad(lambda x:minus_v_n(x,n), [0, s])#, args=n) 
    n = 1
    in_int_1 = quad(lambda x:minus_v_n(x,n), [0, s])#, args=n) 
    delta = - in_int_1 + in_int_0
    return math.exp(delta) * 0.2 #0.05

# пресмятане на \int_0^r {\exp( \int_0^s {(v_1(k) - v_0(k))} dk) * f_0(s)} ds
def ins_integr(r):
    ins_int = quad(lambda s: integrand_1_0(s), [0, r]) 
    return ins_int

# пресмятане на \exp( \int_0^r {(v_2(s) - v_1(s))} ds) * f_1(r)
def integrand_2_1(r):
    n = 2
    in_int_2 = quad(lambda x:minus_v_n(x,n), [0, r])#, args=n) 
    n = 1
    in_int_1 = quad(lambda x:minus_v_n(x,n), [0, r])#, args=n) 
    delta = - in_int_2 + in_int_1
    return math.exp(delta) *0.25 # 0.03

# задаване на  на \exp( \int_0^r {(v_2(s) - v_1(s))} ds) * f_1(r) * \int_0^r {\exp( \int_0^s {(v_1(k) - v_0(k))} dk) * f_0(s)} ds
def final_integrand(r):
    return integrand_2_1(r)*ins_integr(r)
    
mp.dps = 15

start_time = time.time()

step_traject = 0.2
time_int = 5

# начални условия
x_0_0 = 1
x_1_0 = 0.8
x_2_0 = 1.1

# списъци за стойностите при еволюцията във времето
t_list = []
x_2_in_t = []


# пресмятане на интегралите от уравнения (3.2.15)
###  в началния момент
integral_2 = 0
integral_exp21_f1 = 0
integral_21_f1_10_f0 = 0
t_list.append(0)
x_2_in_t.append(x_2_0)


###  в следващи моменти
k = 1
while (k*step_traject)<=time_int:
    
    t_list.append(k*step_traject)

    # интеграл в exp(- \int_0^t {v_2(s)} ds)
    n = 2
    integr_2 = quad( lambda x: minus_v_n(x,n), [(k-1)*step_traject, k*step_traject])#, args=n) 
    integral_2 = integral_2 + integr_2
    
    # \int_^t {\exp( \int_0^r {(v_2(s) - v_1(s))} ds) * f_1(r)} dr
    integr_21 = quad(lambda r: integrand_2_1(r), [(k-1)*step_traject, k*step_traject]) 
    integral_exp21_f1 = integral_exp21_f1 + integr_21
    
    # интеграл в третия член
    integr_21_f1_10_f0 = quad(lambda r: final_integrand(r), [(k-1)*step_traject, k*step_traject]) 
    integral_21_f1_10_f0 = integral_21_f1_10_f0 + integr_21_f1_10_f0
    
    x_2_t = math.exp(integral_2) * (x_2_0 + x_1_0 * integral_exp21_f1 + x_0_0 * integral_21_f1_10_f0)
    
    x_2_in_t.append(x_2_t)

    k += 1


dict = {'time': t_list, 'x_2(t)': x_2_in_t}
   
df = pd.DataFrame(dict) 
    

df.to_csv('x_2_num_gen_func_03.csv') 
# print(t_list)
print(x_2_in_t)


print("--- %s seconds ---" % (time.time() - start_time))