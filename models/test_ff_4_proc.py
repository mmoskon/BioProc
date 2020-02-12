from bioproc.proc_models import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt
    

"""
    TESTING
"""

# simulation parameters
t_end = 200
N = 1000


# model parameters
alpha1 = 34.73
alpha2 = 49.36
alpha3 = 32.73
alpha4 = 49.54
delta1 = 1.93
delta2 = 0.69
Kd = 10.44
n = 4.35
params_ff = (alpha1, alpha2, alpha3, alpha4, delta1, delta2, Kd, n)


# addressing params
alpha = 20
delta = 1
Kd = 80
n = 2

params_addr = (alpha, delta, Kd, n)

points = np.loadtxt('selected_points.txt')
params = points[0]
params_ff = list(params[:8])
params_addr = list(params[8:])

# four-bit counter with external clock
# a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, a3, not_a3, q3, not_q3, a4, not_a4, q4, not_q4, i1, i2, i3, i4, i5, i6, i7, i8
Y0 = np.array([0]*24)
T = np.linspace(0, t_end, N)

Y = odeint(four_bit_processor_ext, Y0, T, args=(params_ff, params_addr))

Y_reshaped = np.split(Y, Y.shape[1], 1)
"""
Q = Y_reshaped[2::4]

for q in Q:
    plt.plot(T, q)
plt.show()
"""
Q1 = Y_reshaped[2]
not_Q1 = Y_reshaped[3]
Q2 = Y_reshaped[6]
not_Q2 = Y_reshaped[7]
Q3 = Y_reshaped[10]
not_Q3 = Y_reshaped[11]
Q4 = Y_reshaped[14]
not_Q4 = Y_reshaped[15]



i1 = Y_reshaped[-8]
i2 = Y_reshaped[-7]
i3 = Y_reshaped[-6]
i4 = Y_reshaped[-5]
i5 = Y_reshaped[-4]
i6 = Y_reshaped[-3]
i7 = Y_reshaped[-2]
i8 = Y_reshaped[-1]



plt.plot(T, Q1, label='q1')
plt.plot(T, Q2, label='q2')
plt.plot(T, Q3, label='q3')
plt.plot(T, Q4, label='q4')

#plt.plot(T, not_Q1, label='not q1')
#plt.plot(T, not_Q2, label='not q2')
plt.plot(T, i1, label='i1')
plt.plot(T, i2, label='i2')
plt.plot(T, i3, label='i3')
plt.plot(T, i4, label='i4')
plt.plot(T, i5, label='i5')
plt.plot(T, i6, label='i6')
plt.plot(T, i7, label='i7')
plt.plot(T, i8, label='i8')



plt.plot(T, get_clock(T),  '--', linewidth=2, label="CLK", color='black', alpha=0.25)

plt.legend()
plt.savefig('figs\\proc4.pdf')
plt.show()
