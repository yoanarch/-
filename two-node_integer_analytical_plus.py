from mpmath import nsum, inf, mp
import math
import pandas as pd

mp.dps = 15

kappa = 3
k = kappa + 1

v_0 = 0.1
v_1 = 0.4
f_0 = 0.3

x_0_0 = 1
x_1_0 = 1.5

alpha = 7
t_max = 5
t_step = 0.2

const_series_1 = nsum(lambda s: -math.exp(-(v_1-v_0)*alpha)*math.factorial(int(k-1-s))*pow((math.factorial(int(k-1))),-1)*pow((-(v_1-v_0)*(alpha)),s-k), [1, k-1])
const_series_2 = nsum(lambda s: pow((-alpha*(v_1-v_0)), s)/(s*math.factorial(int(s))), [1, inf])
const = x_1_0*pow(alpha, 1-k) + x_0_0*f_0*alpha*math.exp(alpha*(v_1-v_0))*pow((-(v_1-v_0)),k-1)*(const_series_1+pow(math.factorial(k-1),-1)*(math.log(abs((v_1-v_0)*alpha))+const_series_2))
# print(const_series_1)
# print(const_series_2)

x_1_t_list = []
t_list = []
t = 0
while t <= t_max+0.1:

    t_list.append(t)

    time_dependent_series_1 = nsum(lambda j: -math.exp(-(v_1-v_0)*(alpha-t))*math.factorial(int(k-1-j))*pow((math.factorial(int(k-1))),-1)*pow((-(v_1-v_0)*(alpha-t)),j-k), [1, k-1])
    time_dependent_series_2  = nsum(lambda j: pow((-(v_1-v_0)*(alpha-t)),j)/(j*math.factorial(int(j))), [1, inf])
    # print(time_dependent_series_1)
    # print(time_dependent_series_2)
    x_1_plus_t = pow((alpha-t),k-1)*math.exp(-v_1*t)*(const - x_0_0*f_0*alpha*math.exp(alpha*(v_1-v_0))*pow(-(v_1-v_0),k-1)*(time_dependent_series_1 + pow(math.factorial(k-1),-1)*(math.log(abs((v_1-v_0)*(alpha-t))) + time_dependent_series_2)))
    x_1_t_list.append(x_1_plus_t)

    t += t_step

print(x_1_t_list)
# print(nsum(lambda k: pow(k*math.factorial(int(k)),-1), [1, inf]))
dict = {'time': t_list, 'x_2(t)': x_1_t_list}
   
df = pd.DataFrame(dict) 
    

df.to_csv('x_1_pl_analyt_kappa3.csv') 