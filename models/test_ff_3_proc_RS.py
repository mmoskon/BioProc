from bioproc.proc_models import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def from_addr_to_i(addr):
    if addr == [0,0,0]:
        return "i1"
    elif addr == [1,0,0]:
        return "i2"
    elif addr == [1,1,0]:
        return "i3"
    elif addr == [1,1,1]:
        return "i4"
    elif addr == [0,1,1]:
        return "i5"
    else:
        return "i6"
    

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

# proteolytic degradation parameters
KM = 5
deltaE = 500
params_ff = params_ff + (deltaE, KM)


# addressing params
alpha = 20
delta = 1
Kd = 80
n = 2

params_addr = (alpha, delta, Kd, n)

# three-bit counter with external clock
# a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, a3, not_a3, q3, not_q3, i1, i2, i3, i4, i5, i6
Y0 = np.array([0]*18)
T = np.linspace(0, t_end, N)
jump_src=[1,1,1]
jump_dst = [0,1,1]

i_src = from_addr_to_i(jump_src)
i_dst = from_addr_to_i(jump_dst)
print(i_src)
print(i_dst)

Y = odeint(three_bit_processor_ext_RS, Y0, T, args=(params_ff, params_addr, jump_src, jump_dst, i_src, i_dst))

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


i1 = Y_reshaped[-6]
i2 = Y_reshaped[-5]
i3 = Y_reshaped[-4]
i4 = Y_reshaped[-3]
i5 = Y_reshaped[-2]
i6 = Y_reshaped[-1]


plt.plot(T, Q1, label='q1')
plt.plot(T, Q2, label='q2')
plt.plot(T, Q3, label='q3')
#plt.plot(T, not_Q1, label='not q1')
#plt.plot(T, not_Q2, label='not q2')
plt.plot(T, i1, label='i1')
plt.plot(T, i2, label='i2')
plt.plot(T, i3, label='i3')
plt.plot(T, i4, label='i4')
plt.plot(T, i5, label='i5')
plt.plot(T, i6, label='i6')

plt.plot(T, get_clock(T),  '--', linewidth=2, label="CLK", color='black', alpha=0.25)

plt.legend()
plt.savefig('figs\\proc3_RS.pdf')
plt.show()
