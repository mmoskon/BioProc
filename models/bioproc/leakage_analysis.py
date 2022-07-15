import numpy as np
from hill_functions import *   
from scipy.integrate import odeint
import matplotlib.pyplot as plt
 
"""
FLIP-FLOP MODELS
"""
# MASTER-SLAVE D FLIP-FLOP MODEL
def ff_ode_model(Y, T, params): 
    
    a, not_a, q, not_q, d, clk = Y
    alpha1, alpha2, alpha3, alpha4, delta1, delta2, Kd, n, leak = params

    #delta1 += 2*delta1*leak
    #delta2 += 2*delta2*leak

    da_dt     = alpha1*(leak + (pow(d/Kd, n)/(1 + pow(d/Kd, n) + pow(clk/Kd, n) + pow(d/Kd, n)*pow(clk/Kd, n)))) + alpha2*(leak + (1/(1 + pow(not_a/Kd, n)))) - delta1 *a 
    dnot_a_dt = alpha1*(leak + (1/(1 + pow(d/Kd, n) + pow(clk/Kd, n) + pow(d/Kd, n)*pow(clk/Kd, n)))) + alpha2*(leak + (1/(1 + pow(a/Kd, n)))) - delta1*not_a   
    dq_dt     = alpha3*(leak + ((pow(a/Kd, n)*pow(clk/Kd, n))/(1 + pow(a/Kd, n) + pow(clk/Kd, n) + pow(a/Kd, n)*pow(clk/Kd, n)))) + alpha4*(leak + (1/(1 + pow(not_q/Kd, n)))) - delta2*q  
    dnot_q_dt = alpha3*(leak + ((pow(not_a/Kd, n)*pow(clk/Kd, n))/(1 + pow(not_a/Kd, n) + pow(clk/Kd, n) + pow(not_a/Kd, n)*pow(clk/Kd, n)))) + alpha4*(leak + (1/(1 + pow(q/Kd, n)))) - delta2*not_q   


    return np.array([da_dt, dnot_a_dt, dq_dt, dnot_q_dt]) 

"""
JOHSON COUNTER MODELS 
"""
	
def one_bit_model_switching(Y, T, params):
    a, not_a, q, not_q= Y

    clk = get_clock(T) 

    d = not_q
    Y_FF1 = [a, not_a, q, not_q, d, clk]

    dY = ff_ode_model(Y_FF1, T, params)

    return dY

def one_bit_model_holding(Y, T, params):
    a, not_a, q, not_q= Y

    clk = get_clock(T) 

    d = q
    Y_FF1 = [a, not_a, q, not_q, d, clk]

    dY = ff_ode_model(Y_FF1, T, params)

    return dY

##########
def run_MODEL(y0, T, params, dt = 0.001, leak = 0, switching=False, plot=False):
    if switching:
        model = one_bit_model_switching
    else:
        model = one_bit_model_holding

    N = int(T/dt)    
    ts = np.linspace(0, T, N)   
    #clk = get_clock(ts)

    params2 = params + [leak]

    Y = odeint(model, y0, ts, args=(params2,))

    q = Y[:,2]
    not_q = Y[:,3]

    if plot:       
        plt.plot(ts, q)
        plt.plot(ts, not_q)
        plt.show()

    q = np.array(q[len(q)//2:])
    not_q = np.array(not_q[len(not_q)//2:])

    return max(q-not_q)

# alpha1, alpha2, alpha3, alpha4, delta1, delta2, Kd, n = params

params = [1.558143736000000068e+01, 2.521915470000000159e+00, 2.435646121999999991e+01, 3.612683424000000088e+01, 2.022238700000000000e-01, 4.961420200000000169e-01, 1.701689959000000130e+01, 5.000000000000000000e+00]
#params = [2.226719939999999909e+01, 1.189967680000000083e+00, 1.963390148999999951e+01, 5.000000000000000000e+01, 1.363039000000000056e-01, 6.367345500000000103e-01, 1.557853451999999983e+01, 5.000000000000000000e+00]
#params = [1.674008136000000135e+01, 1.629092289999999998e+00, 2.185285782000000054e+01, 5.000000000000000000e+01, 1.363039000000000056e-01, 6.367345500000000103e-01, 2.060897531000000171e+01, 5.000000000000000000e+00]

#alpha1_s = 

#alpha1 = 0.8508*3600
#alpha2 = 1.5299*3600
#alpha3 = 0.3431*3600
#alpha4 = 1.5299*3600
#delta1=0.0036*3600
#delta2=0.0036*3600
#Kd=0.05
#n = 4
#params = (alpha1, alpha2, alpha3, alpha4, delta1, delta2, Kd, n)

#d = run_MODEL([10,0,10,0], 250, params, switching=True, leak = 0.2, plot=True)
#print(d)

def analyse_leak(y0, T, params, dt=0.001, switching=False):
    L = np.arange(0, 0.26, 0.01)
    DQ = []

    for l in L:
        d = run_MODEL([10,0,10,0], 250, params, switching=False, leak = l, plot=False)
        DQ.append(d)

    DQ = np.array(DQ)
    DQ /= DQ[0]
    plt.plot(L, DQ)
    plt.xlabel("Relative leakage rate")
    plt.ylabel("Relative state distance")
    plt.savefig("leakage.pdf", bbox_inches = 'tight')
    plt.show()

analyse_leak([10,0,10,0], 250, params, switching=False)