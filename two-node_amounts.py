from mpmath import quad, mp
import numpy as np
import math
import time
import pandas as pd

mp.dps = 15

# def minus_v_n(x,n):

    
#     if n == 0:
#         arg = -0.7* math.sin(math.pi*x) - 0.2 + 0.0

#     if n == 1:
#         arg = - 0.3* math.sin(math.pi*x) - 0.25 + 0.0

#     return arg

def minus_v_n( x, n):
    if n == 0:
        arg = - 0.4-0.25 + pow( 7-x, -1)

    kappa_ = 1.1
    if n == 1:
        arg = - 0.2 - kappa_*pow(7-x, -1)  

    return arg      



# пресмятане на \exp( \int_0^r {(v_1(s) - v_0(s))} ds) * f_0(r)
def integrand(r):
    n = 0 
    in_int_0 = quad(lambda x: minus_v_n(x,n), [0, r]) 
    n = 1
    in_int_1 = quad(lambda x: minus_v_n(x,n), [0, r]) 
    delta = - in_int_1 + in_int_0

    return math.exp(delta) * 0.25

start_time = time.time()


step_traject = 0.25
time_max = 5

# начални условия 
x_0_0 = 1
x_1_0 = 1.5


# списъци за стойностите при еволюцията във времето

t_list = []

x_0_in_t = []


x_1_in_t = []



# пресмятане на интегралите от уравнения (3.2.1)

integral_0 = 0
t_list.append(0)
x_0_in_t.append(x_0_0 * math.exp(integral_0))
x_1_in_t.append(x_1_0)


integral_1 = 0
integral_1_exp_f0 = 0

k = 1
while (k*step_traject) <= time_max+0.1:

    t_list.append(k*step_traject)
    
    #  за x_0(t)

    n = 0
    integr = quad(lambda x: minus_v_n(x,n), [(k-1)*step_traject, k*step_traject]) 
    integral_0 = integral_0 + integr
    x_0_in_t.append(x_0_0 * math.exp(integral_0))



    # за x_1(t)
  
    # exp(- \int_0^t {v_1(s)} ds)
    n = 1
    integr_1 = quad(lambda x: minus_v_n(x,n), [(k-1)*step_traject, k*step_traject]) 
    integral_1 = integral_1 + integr_1

    # интеграл от експонента
    integr_1_exp_f0 = quad(lambda x: integrand(x), [(k-1)*step_traject, k*step_traject]) 
    integral_1_exp_f0 = integral_1_exp_f0 + integr_1_exp_f0

    x_1_k = math.exp(integral_1)*(x_1_0 + x_0_0 * integral_1_exp_f0)
 #    print(math.exp(integral_1))

    x_1_in_t.append(x_1_k)

    k += 1


#print(t_list)
print(x_0_in_t)
print(x_1_in_t)

dict = {'time': t_list, 'x_0(t)':x_0_in_t,'x_2(t)': x_1_in_t}
   
df = pd.DataFrame(dict) 
    

df.to_csv('x_0_x_1_num_typical_plus.csv') 


print("--- %s seconds ---" % (time.time() - start_time))