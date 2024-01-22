from mpmath import nsum, inf , mp
import math
import pandas as pd

mp.dps = 15

kappa = -3
k = - kappa - 1
lambda_ = - 4
l = kappa - lambda_

v_0 = 0.9
v_1 = 0.2
v_2 = 0.1

f_0 = 0.05
f_1 = 0.03

x_0_0 = 1
x_1_0 = 0.8
x_2_0 = 1.1

alpha = 7
t_max = 5
t_step = 0.2

const_series_1 = nsum(lambda j: math.factorial(l)*pow(math.factorial(int(j)),-1)*pow((v_2-v_1)*alpha, j-l-1), [0, l])
const_series_2 = nsum(lambda j: nsum(lambda i: math.factorial(k)*math.factorial(int(l+j))*pow(math.factorial(int(j))*math.factorial(int(i)),-1)\
                * pow((v_1-v_0), j-k-1)*pow(v_2-v_0, i-l-j-1)* pow(alpha,i) ,[0, l+j]), [0, k])
const_series_3 = nsum(lambda j: math.factorial(l)*pow(math.factorial(int(j)),-1)*pow((v_2-v_1), j-l-1)*pow(alpha,j), [0, l])
const_series_4 = nsum(lambda j: math.factorial(k)*pow(math.factorial(int(j)),-1)*pow((v_1-v_0), j-k-1)*pow(alpha,j), [0, k])

constant_ = x_2_0*pow(alpha, k+l+1) - x_1_0*f_1*pow(alpha, k+l+2)* const_series_1 - x_0_0*f_0*f_1*alpha*const_series_2 \
     + x_0_0*f_0*f_1*alpha*const_series_3*const_series_4


x_2_t_list = []
t_list = []
t = 0
while t <= t_max+0.1:

    t_list.append(t)

    time_dependent_series_1 = nsum(lambda j: math.factorial(l)*pow(math.factorial(int(j)),-1)*pow((v_2-v_1), j-l-1)*pow(alpha-t,j), [0, l])
    time_dependent_series_2 = nsum(lambda j: nsum(lambda i: math.factorial(k)*math.factorial(int(l+j))*pow(math.factorial(int(j))*math.factorial(int(i)),-1)\
                * pow((v_1-v_0), j-k-1)*pow(v_2-v_0, i-l-j-1)* pow(alpha-t,i) ,[0, l+j]), [0, k])
    time_dependent_series_3 = nsum(lambda j: math.factorial(l)*pow(math.factorial(int(j)),-1)*pow((v_2-v_1), j-l-1)*pow(alpha-t,j), [0, l])
    
    x_2_plus_t = constant_ + x_1_0*f_1*pow(alpha, k+1)* math.exp((v_2-v_1)*t)*time_dependent_series_1 + x_0_0*f_0*f_1*alpha*math.exp((v_2-v_0)*t)*time_dependent_series_2 \
     - x_0_0*f_0*f_1*alpha*math.exp((v_2-v_1)*t)*time_dependent_series_3*const_series_4
    
    x_2_t_list.append(pow(alpha-t,-l-k-1)*math.exp(-v_2*t)*x_2_plus_t)

    t += t_step

print(x_2_t_list)

dict = {'time': t_list, 'x_2(t)': x_2_t_list}
   
df = pd.DataFrame(dict) 
    

df.to_csv('x_2_min_analyt_kappa-3_lambda-4.csv') 