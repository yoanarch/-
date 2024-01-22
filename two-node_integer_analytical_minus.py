from mpmath import nsum, mp
import math
import pandas as pd

mp.dps = 15

kappa = - 3
k = - kappa - 1

v_0 = 0.1
v_1 = 0.4
f_0 = 0.3

x_0_0 = 1
x_1_0 = 1.5

alpha = 7
t_max = 5
t_step = 0.2

const_series = nsum(lambda j: pow(-1, k-j)*math.factorial(k)*pow(math.factorial(int(j)), -1)*pow(-(v_1-v_0)*alpha,j), [0,k])
const = x_1_0*pow(alpha,(k+1)) + x_0_0*f_0*alpha*pow(-(v_1-v_0), -k-1)*const_series

x_1_t_list = []
t_list = []
t = 0
while t <= t_max+0.1:

    t_list.append(t)

    time_dependent_series = nsum(lambda j: pow(-1, k-j)*math.factorial(int(k))*pow(math.factorial(int(j)), -1)*pow(-(v_1-v_0)*(alpha-t),j), [0,k])
    x_1_plus_t = pow((alpha-t),-k-1)*math.exp(-v_1*t)*(const - x_0_0*f_0*alpha*math.exp(t*(v_1-v_0))*pow(-(v_1-v_0),-k-1)*time_dependent_series)
    x_1_t_list.append(x_1_plus_t)

    t += t_step

print(x_1_t_list)

dict = {'time': t_list, 'x_2(t)': x_1_t_list}
   
df = pd.DataFrame(dict) 
    

df.to_csv('x_1_min_analyt_kappa-3.csv') 