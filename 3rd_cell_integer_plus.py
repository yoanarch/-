from mpmath import nsum, inf, mp 
import math
import pandas as pd

mp.dps = 15

kappa = -3
k = - kappa - 1
lambda_ = 2
l = lambda_ - kappa

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


const_series_1 = nsum(lambda s: -math.exp(-(v_2-v_1)*alpha)*math.factorial(int(l-1-s))*pow((math.factorial(int(l-1))),-1)*pow((-(v_2-v_1)*(alpha)),s-l), [1, l-1])
const_series_2 = nsum(lambda s: pow((-alpha*(v_2-v_1)), s)/(s*math.factorial(int(s))), [1, inf])
const_series_3 = nsum(lambda s: pow((-alpha*(v_2-v_0)), s)/(s*math.factorial(int(s))), [1, inf])
const_series_4 = nsum(lambda j: math.factorial(k)*pow(math.factorial(int(j)),-1)*pow((v_0-v_2),l-j-1)*pow(v_1-v_0, j-k-1)\
                      * (nsum(lambda r: -math.exp(-(v_2-v_0)*alpha)*math.factorial(int(l-j-1-r))*pow((math.factorial(int(l-j-1))),-1)*pow((-(v_2-v_0)*(alpha)),r-l+j), [1, l-j-1] )\
                         + pow(math.factorial(int(l-j-1)),-1)*(math.log(abs((v_2-v_0)*alpha))+const_series_3)), [0, k])
const_series_5 = nsum(lambda j: math.factorial(k)*pow(math.factorial(int(j)),-1)*pow(alpha, j+1)*pow(v_1-v_0, j-k-1), [0, k])

constant_ = x_2_0 * pow(alpha, 1+k-l) + x_1_0 * f_1 * pow(alpha, k+1) * math.exp(alpha*(v_2-v_1))* pow(v_1-v_2, l-1) \
                * (const_series_1 + pow(math.factorial(l-1),-1)*(math.log(abs(alpha*(v_2-v_1)))+const_series_2))
constant_ += x_0_0 * f_0 * f_1 * alpha * math.exp(alpha*(v_2-v_0)) * const_series_4
constant_ += -x_0_0 * f_0 * f_1 * pow(v_1 - v_2, l-1) * math.exp(alpha*(v_2-v_1)) * const_series_5 *  (const_series_1 \
                    + pow(math.factorial(l-1),-1)*(math.log(abs((v_2-v_1)*(alpha)))+const_series_2))

#print(const_series_5)

x_2_t_list = []
t_list = []
t = 0
while t <= t_max+0.1:

    t_list.append(t)

    time_dependent_series_1 = nsum(lambda s: -math.exp(-(v_2-v_1)*(alpha-t))*math.factorial(int(l-1-s))*pow((math.factorial(int(l-1))),-1)*pow((-(v_2-v_1)*(alpha-t)),s-l), [1, l-1])
    time_dependent_series_2 = nsum(lambda s: pow((-(alpha-t)*(v_2-v_1)), s)/(s*math.factorial(int(s))), [1, inf])
    time_dependent_series_3 = nsum(lambda s: pow((-(alpha-t)*(v_2-v_0)), s)/(s*math.factorial(int(s))), [1, inf])
    time_dependent_series_4 = nsum(lambda j: math.factorial(k)*pow(math.factorial(int(j)),-1)*pow((v_0-v_2),l-j-1)*pow(v_1-v_0, j-k-1)\
                        * (nsum(lambda r: -math.exp(-(v_2-v_0)*(alpha-t))*math.factorial(int(l-j-1-r))*pow((math.factorial(int(l-j-1))),-1)*pow((-(v_2-v_0)*(alpha-t)),r-l+j), [1, l-j-1] )\
                            + pow(math.factorial(int(l-j-1)),-1)*(math.log(abs((v_2-v_0)*(alpha-t)))+time_dependent_series_3)), [0, k])
    
    #print(time_dependent_series_4)
    x_2_plus_t = constant_ - x_1_0*f_1*pow(alpha,k+1)*math.exp(alpha*(v_2-v_1))*pow(v_1-v_2,l-1)\
                    * (time_dependent_series_1 + pow(math.factorial(l-1),-1)*(math.log(abs((alpha-t)*(v_2-v_1)))+time_dependent_series_2))
    x_2_plus_t += -1*x_0_0 * f_0 * f_1 * alpha * math.exp(alpha*(v_2-v_0)) * time_dependent_series_4
    x_2_plus_t += x_0_0 * f_0 * f_1 * pow(v_1 - v_2, l-1) * math.exp(alpha*(v_2-v_1)) * const_series_5 * (time_dependent_series_1 \
                    + pow(math.factorial(l-1),-1)*(math.log(abs((v_2-v_1)*(alpha-t)))+time_dependent_series_2))
    
    x_2_t_list.append(pow(alpha-t,l-k-1)*math.exp(-v_2*t)*x_2_plus_t)

    t += t_step


dict = {'time': t_list, 'x_2(t)': x_2_t_list}
   
df = pd.DataFrame(dict) 
    

df.to_csv('x_2_pl_analyt_kappa-3_lambda2.csv') 
print(x_2_t_list)
print(t_list)
