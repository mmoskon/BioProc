import numpy as np
from bioproc.hill_functions import *   
 
"""
FLIP-FLOP MODELS
"""
# MASTER-SLAVE D FLIP-FLOP QSSA MODEL
def ff_stochastic_model(Y, T, params, omega):
	p = np.zeros(12)  

	a, not_a, q, not_q, d, clk = Y
	alpha1, alpha2, alpha3, alpha4, delta1, delta2, Kd, n = params
    
	p[0] = alpha1*(pow(d/(Kd*omega), n)/(1 + pow(d/(Kd*omega), n) + pow(clk/(Kd*omega), n) + pow(d/(Kd*omega), n)*pow(clk/(Kd*omega), n)))*omega   
	p[1] = alpha2*(1/(1 + pow(not_a/(Kd*omega), n)))*omega    
	p[2] = delta1*a  
	p[3] = alpha1*(1/(1 + pow(d/(Kd*omega), n) + pow(clk/(Kd*omega), n) + pow(d/(Kd*omega), n)*pow(clk/(Kd*omega), n)))*omega   
	p[4] = alpha2*(1/(1 + pow(a/(Kd*omega), n)))*omega   
	p[5] = delta1*not_a 
	p[6] = alpha3*((pow(a/(Kd*omega), n)*pow(clk/(Kd*omega), n))/(1 + pow(a/(Kd*omega), n) + pow(clk/(Kd*omega), n) + pow(a/(Kd*omega), n)*pow(clk/(Kd*omega), n)))*omega
	p[7] = alpha4*(1/(1 + pow(not_q/(Kd*omega), n)))*omega   
	p[8] = delta2*q 
	p[9] = alpha3*((pow(not_a/(Kd*omega), n)*pow(clk/(Kd*omega), n))/(1 + pow(not_a/(Kd*omega), n) + pow(clk/(Kd*omega), n) + pow(not_a/(Kd*omega), n)*pow(clk/(Kd*omega), n)))*omega   
	p[10] = alpha4*(1/(1 + pow(q/(Kd*omega), n)))*omega 
	p[11] = delta2*not_q 

	#propensities     
	return p   
	
	
# MASTER-SLAVE D FLIP-FLOP MODEL
def ff_ode_model(Y, T, params): 
    
    a, not_a, q, not_q, d, clk = Y
    alpha1, alpha2, alpha3, alpha4, delta1, delta2, Kd, n = params

    da_dt     = alpha1*(pow(d/Kd, n)/(1 + pow(d/Kd, n) + pow(clk/Kd, n) + pow(d/Kd, n)*pow(clk/Kd, n))) + alpha2*(1/(1 + pow(not_a/Kd, n))) - delta1 *a 
    dnot_a_dt = alpha1*(1/(1 + pow(d/Kd, n) + pow(clk/Kd, n) + pow(d/Kd, n)*pow(clk/Kd, n))) + alpha2*(1/(1 + pow(a/Kd, n))) - delta1*not_a   
    dq_dt     = alpha3*((pow(a/Kd, n)*pow(clk/Kd, n))/(1 + pow(a/Kd, n) + pow(clk/Kd, n) + pow(a/Kd, n)*pow(clk/Kd, n))) + alpha4*(1/(1 + pow(not_q/Kd, n))) - delta2*q  
    dnot_q_dt = alpha3*((pow(not_a/Kd, n)*pow(clk/Kd, n))/(1 + pow(not_a/Kd, n) + pow(clk/Kd, n) + pow(not_a/Kd, n)*pow(clk/Kd, n))) + alpha4*(1/(1 + pow(q/Kd, n))) - delta2*not_q   


    return np.array([da_dt, dnot_a_dt, dq_dt, dnot_q_dt]) 

# FF MODEL WITH ASYNCHRONOUS RESET AND SET
# dodana parametra deltaE, KM
# dodani vhodni spremenljivki RESET in SET
# dodano 23. 1. 2020
def ff_ode_model_RS(Y, T, params): 
    
    a, not_a, q, not_q, d, clk, RESET, SET = Y

    repress_both = True

    if repress_both:
            sum_one = a + q
            sum_zero = not_a + not_q
    
    alpha1, alpha2, alpha3, alpha4, delta1, delta2, Kd, n, deltaE, KM = params


    da_dt     = alpha1*(pow(d/Kd, n)/(1 + pow(d/Kd, n) + pow(clk/Kd, n) + pow(d/Kd, n)*pow(clk/Kd, n))) + alpha2*(1/(1 + pow(not_a/Kd, n))) - delta1 *a 

    #deltaE = delta1
    if repress_both:
        da_dt += -a*(deltaE*RESET/(KM+sum_one))
    else:
        da_dt += -a*(deltaE*RESET/(KM+a))


    dnot_a_dt = alpha1*(1/(1 + pow(d/Kd, n) + pow(clk/Kd, n) + pow(d/Kd, n)*pow(clk/Kd, n))) + alpha2*(1/(1 + pow(a/Kd, n))) - delta1*not_a
    if repress_both:
        dnot_a_dt += -not_a*(deltaE*SET/(KM+sum_zero))
    else:
        dnot_a_dt += -not_a*(deltaE*SET/(KM+not_a))    


    #deltaE = delta2
    dq_dt     = alpha3*((pow(a/Kd, n)*pow(clk/Kd, n))/(1 + pow(a/Kd, n) + pow(clk/Kd, n) + pow(a/Kd, n)*pow(clk/Kd, n))) + alpha4*(1/(1 + pow(not_q/Kd, n))) - delta2*q
    if repress_both:
        dq_dt += -q*(deltaE*RESET/(KM+sum_one))
    
    dnot_q_dt = alpha3*((pow(not_a/Kd, n)*pow(clk/Kd, n))/(1 + pow(not_a/Kd, n) + pow(clk/Kd, n) + pow(not_a/Kd, n)*pow(clk/Kd, n))) + alpha4*(1/(1 + pow(q/Kd, n))) - delta2*not_q   
    if repress_both:
        dnot_q_dt += -not_q*(deltaE*SET/(KM+sum_zero))
   
    return np.array([da_dt, dnot_a_dt, dq_dt, dnot_q_dt]) 


"""
ADRESSING MODELS
"""

# ADDRESSING 1-BIT QSSA MODEL
def addressing_stochastic_one_bit_model(Y, T, params, omega):   
    alpha, delta, Kd, n = params
    _,_, q1, not_q1, i1, i2 = Y  
    p = np.zeros(4) 
	
    p[0] = alpha*activate_1(not_q1, Kd*omega, n)*omega  
    p[1] = delta*i1  
    p[2] = alpha*activate_1(q1, Kd*omega, n)*omega 
    p[3] = delta*i2
	
    #propensities    
    return p

# ADDRESSING 2-BIT QSSA MODEL 
def addressing_stochastic_two_bit_model(Y, T, params, omega):   
    alpha, delta, Kd, n = params
    _, _, q1, not_q1, _, _, q2, not_q2, i1, i2, i3, i4 = Y
    p = np.zeros(8)  
	
    p[0] = alpha * activate_2(not_q1, not_q2, Kd*omega, n)*omega 
    p[1] = delta * i1
    p[2] = alpha * activate_2(q1, not_q2, Kd*omega, n)*omega  
    p[3] = delta * i2 
    p[4] = alpha * activate_2(q1, q2, Kd*omega, n)*omega    
    p[5] = delta * i3  
    p[6] = alpha * activate_2(not_q1, q2, Kd*omega, n)*omega
    p[7] = delta * i4    
			
    #propensities    
    return p 

# ADDRESSING 3-BIT QSSA MODEL 
def addressing_stochastic_three_bit_model(Y, T, params, omega):   
    alpha, delta, Kd, n = params
    _, _, q1, not_q1, _, _, q2, not_q2, _, _, q3, not_q3, i1, i2, i3, i4, i5, i6 = Y 
    p = np.zeros(12)  
	
    p[0] = alpha * activate_2(not_q1, not_q3, Kd*omega, n)*omega
    p[1] = delta * i1
    p[2] = alpha * activate_2(q1, not_q2, Kd*omega, n)*omega 
    p[3] = delta * i2
    p[4] = alpha * activate_2(q2, not_q3, Kd*omega, n)*omega
    p[5] = delta * i3
    p[6] = alpha * activate_2(q1, q3, Kd*omega, n)*omega
    p[7] = delta * i4
    p[8] = alpha * activate_2(not_q1, q2, Kd*omega, n)*omega
    p[9] = delta * i5  
    p[10] = alpha * activate_2(not_q2, q3, Kd*omega, n)*omega  
    p[11] = delta * i6    	
	
	#propensities      
    return p   	 	

# ONE BIT ADDRESSING MODEL SIMPLE
def one_bit_simple_addressing_ode_model(Y, T, params):
    alpha, delta, Kd, n = params
    
    q1, not_q1, i1, i2 = Y

    di1_dt = alpha * activate_1(not_q1, Kd, n) - delta * i1
    di2_dt = alpha * activate_1(q1, Kd, n) - delta * i2
    
    return np.array([di1_dt, di2_dt])

   
# TWO BIT ADDRESSING MODEL SIMPLE
def two_bit_simple_addressing_ode_model(Y, T, params):
    alpha, delta, Kd, n = params
    
    q1, not_q1, q2, not_q2, i1, i2, i3, i4 = Y

    di1_dt = alpha * activate_2(not_q1, not_q2, Kd, n) - delta * i1
    di2_dt = alpha * activate_2(q1, not_q2, Kd, n) - delta * i2
    di3_dt = alpha * activate_2(q1, q2, Kd, n) - delta * i3
    di4_dt = alpha * activate_2(not_q1, q2, Kd, n) - delta * i4

    return np.array([di1_dt, di2_dt, di3_dt, di4_dt])

# THREE BIT ADDRESSING MODEL SIMPLE
def three_bit_simple_addressing_ode_model(Y, T, params):
    alpha, delta, Kd, n = params
    
    q1, not_q1, q2, not_q2, q3, not_q3, i1, i2, i3, i4, i5, i6 = Y

    di1_dt = alpha * activate_2(not_q1, not_q3, Kd, n) - delta * i1
    di2_dt = alpha * activate_2(q1, not_q2, Kd, n) - delta * i2
    di3_dt = alpha * activate_2(q2, not_q3, Kd, n) - delta * i3
    di4_dt = alpha * activate_2(q1, q3, Kd, n) - delta * i4
    di5_dt = alpha * activate_2(not_q1, q2, Kd, n) - delta * i5
    di6_dt = alpha * activate_2(not_q2, q3, Kd, n) - delta * i6

    return np.array([di1_dt, di2_dt, di3_dt, di4_dt, di5_dt, di6_dt])

# FOUR BIT ADDRESSING MODEL SIMPLE
def four_bit_simple_addressing_ode_model(Y, T, params):
    alpha, delta, Kd, n = params
    
    q1, not_q1, q2, not_q2, q3, not_q3, q4, not_q4, i1, i2, i3, i4, i5, i6, i7, i8 = Y

    di1_dt = alpha * activate_2(not_q1, not_q4, Kd, n) - delta * i1
    di2_dt = alpha * activate_2(q1, not_q2, Kd, n) - delta * i2
    di3_dt = alpha * activate_2(q2, not_q3, Kd, n) - delta * i3
    di4_dt = alpha * activate_2(q3, not_q4, Kd, n) - delta * i4
    
    di5_dt = alpha * activate_2(q1, q4, Kd, n) - delta * i5
    di6_dt = alpha * activate_2(not_q1, q2, Kd, n) - delta * i6
    di7_dt = alpha * activate_2(not_q2, q3, Kd, n) - delta * i7
    di8_dt = alpha * activate_2(not_q3, q4, Kd, n) - delta * i8


    return np.array([di1_dt, di2_dt, di3_dt, di4_dt, di5_dt, di6_dt, di7_dt, di8_dt])


# FIVE BIT ADDRESSING MODEL SIMPLE
def five_bit_simple_addressing_ode_model(Y, T, params):
    alpha, delta, Kd, n = params
    
    q1, not_q1, q2, not_q2, q3, not_q3, q4, not_q4, q5, not_q5, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10 = Y

    di1_dt = alpha * activate_2(not_q1, not_q5, Kd, n) - delta * i1
    di2_dt = alpha * activate_2(q1, not_q2, Kd, n) - delta * i2
    di3_dt = alpha * activate_2(q2, not_q3, Kd, n) - delta * i3
    di4_dt = alpha * activate_2(q3, not_q4, Kd, n) - delta * i4
    di5_dt = alpha * activate_2(q4, not_q5, Kd, n) - delta * i5

    di6_dt = alpha * activate_2(q1, q5, Kd, n) - delta * i6
    
    di7_dt = alpha * activate_2(not_q1, q2, Kd, n) - delta * i7
    di8_dt = alpha * activate_2(not_q2, q3, Kd, n) - delta * i8
    di9_dt = alpha * activate_2(not_q3, q4, Kd, n) - delta * i9
    di10_dt = alpha * activate_2(not_q4, q5, Kd, n) - delta * i10

    return np.array([di1_dt, di2_dt, di3_dt, di4_dt, di5_dt, di6_dt, di7_dt, di8_dt, di9_dt, di10_dt])



"""
JOHSON COUNTER MODELS 
"""
	
# TOP MODEL (JOHNSON): ONE BIT MODEL WITH EXTERNAL CLOCK
def one_bit_model(Y, T, params):
    a, not_a, q, not_q= Y

    clk = get_clock(T) 

    d = not_q
    Y_FF1 = [a, not_a, q, not_q, d, clk]

    dY = ff_ode_model(Y_FF1, T, params)

    return dY

# TOP MODEL (JOHNSON): TWO BIT MODEL WITH EXTERNAL CLOCK    
def two_bit_model(Y, T, params): 
    a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2 = Y

    clk = get_clock(T) 

    d1 = not_q2
    d2 = q1
    
    Y_FF1 = [a1, not_a1, q1, not_q1, d1, clk]
    Y_FF2 = [a2, not_a2, q2, not_q2, d2, clk]

    dY1 = ff_ode_model(Y_FF1, T, params)
    dY2 = ff_ode_model(Y_FF2, T, params)

    dY = np.append(dY1, dY2)

    return dY

# TOP MODEL (JOHNSON): THREE BIT MODEL WITH EXTERNAL CLOCK    
def three_bit_model(Y, T, params):
    a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, a3, not_a3, q3, not_q3 = Y

    clk = get_clock(T) 

    d1 = not_q3
    d2 = q1
    d3 = q2
    
    Y_FF1 = [a1, not_a1, q1, not_q1, d1, clk]
    Y_FF2 = [a2, not_a2, q2, not_q2, d2, clk]
    Y_FF3 = [a3, not_a3, q3, not_q3, d3, clk]

    dY1 = ff_ode_model(Y_FF1, T, params)
    dY2 = ff_ode_model(Y_FF2, T, params)
    dY3 = ff_ode_model(Y_FF3, T, params)

    dY = np.append(np.append(dY1, dY2), dY3)

    return dY

# TOP MODEL (JOHNSON): FOUR BIT MODEL WITH EXTERNAL CLOCK    
def four_bit_model(Y, T, params):
    a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, a3, not_a3, q3, not_q3, a4, not_a4, q4, not_q4 = Y

    clk = get_clock(T) 

    d1 = not_q4
    d2 = q1
    d3 = q2
    d4 = q3

    Y_FF1 = [a1, not_a1, q1, not_q1, d1, clk]
    Y_FF2 = [a2, not_a2, q2, not_q2, d2, clk]
    Y_FF3 = [a3, not_a3, q3, not_q3, d3, clk]
    Y_FF4 = [a4, not_a4, q4, not_q4, d4, clk]

    dY1 = ff_ode_model(Y_FF1, T, params)
    dY2 = ff_ode_model(Y_FF2, T, params)
    dY3 = ff_ode_model(Y_FF3, T, params)
    dY4 = ff_ode_model(Y_FF4, T, params)

    dY = np.append(np.append(np.append(dY1, dY2), dY3), dY4)

    return dY


# TOP MODEL (JOHNSON): FIVE BIT MODEL WITH EXTERNAL CLOCK    
def five_bit_model(Y, T, params):
    a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, a3, not_a3, q3, not_q3, a4, not_a4, q4, not_q4, a5, not_a5, q5, not_q5 = Y

    clk = get_clock(T) 

    d1 = not_q5
    d2 = q1
    d3 = q2
    d4 = q3
    d5 = q4

    Y_FF1 = [a1, not_a1, q1, not_q1, d1, clk]
    Y_FF2 = [a2, not_a2, q2, not_q2, d2, clk]
    Y_FF3 = [a3, not_a3, q3, not_q3, d3, clk]
    Y_FF4 = [a4, not_a4, q4, not_q4, d4, clk]
    Y_FF5 = [a5, not_a5, q5, not_q5, d5, clk]

    dY1 = ff_ode_model(Y_FF1, T, params)
    dY2 = ff_ode_model(Y_FF2, T, params)
    dY3 = ff_ode_model(Y_FF3, T, params)
    dY4 = ff_ode_model(Y_FF4, T, params)
    dY5 = ff_ode_model(Y_FF5, T, params)

    dY = np.append(np.append(np.append(np.append(dY1, dY2), dY3), dY4), dY5)

    return dY


"""
JOHSON COUNTER MODELS THAT USE FLIP-FLOPS WITH ASYNCRHONOUS SET/RESET
dodano 23. 1. 2020
"""
	
# TOP MODEL (JOHNSON): ONE BIT MODEL WITH EXTERNAL CLOCK AND FLIP-FLOPS WITH ASYNCRHONOUS SET/RESET
def one_bit_model_RS(Y, T, params):
    a, not_a, q, not_q, R, S = Y

    clk = get_clock(T) 

    d = not_q
    Y_FF1 = [a, not_a, q, not_q, d, clk, R, S]

    dY = ff_ode_model_RS(Y_FF1, T, params)
    
    return dY

# TOP MODEL (JOHNSON): TWO BIT MODEL WITH EXTERNAL CLOCK AND FLIP-FLOPS WITH ASYNCRHONOUS SET/RESET    
def two_bit_model_RS(Y, T, params): 
    a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, R1, S1, R2, S2 = Y

    clk = get_clock(T) 

    d1 = not_q2
    d2 = q1
    
    Y_FF1 = [a1, not_a1, q1, not_q1, d1, clk, R1, S1]
    Y_FF2 = [a2, not_a2, q2, not_q2, d2, clk, R2, S2]

    dY1 = ff_ode_model_RS(Y_FF1, T, params)
    dY2 = ff_ode_model_RS(Y_FF2, T, params)

    dY = np.append(dY1, dY2)

    return dY

# TOP MODEL (JOHNSON): THREE BIT MODEL WITH EXTERNAL CLOCK AND FLIP-FLOPS WITH ASYNCRHONOUS SET/RESET    
def three_bit_model_RS(Y, T, params):
    a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, a3, not_a3, q3, not_q3, R1, S1, R2, S2, R3, S3 = Y

    clk = get_clock(T) 

    d1 = not_q3
    d2 = q1
    d3 = q2

       
    Y_FF1 = [a1, not_a1, q1, not_q1, d1, clk, R1, S1]
    Y_FF2 = [a2, not_a2, q2, not_q2, d2, clk, R2, S2]
    Y_FF3 = [a3, not_a3, q3, not_q3, d3, clk, R3, S3]

    dY1 = ff_ode_model_RS(Y_FF1, T, params)
    dY2 = ff_ode_model_RS(Y_FF2, T, params)
    dY3 = ff_ode_model_RS(Y_FF3, T, params)

    dY = np.append(np.append(dY1, dY2), dY3)

    return dY



"""
PROCESSOR MODEL


!!!OPTIMIZACIJA NAD TEMI MODELI!!!


"""

# TOP MODEL OF PROCESSOR WITH ONE BIT ADDRESSING
def one_bit_processor_ext(Y, T, params_johnson, params_addr):
    a1, not_a1, q1, not_q1, i1, i2  = Y

    Y_johnson = [a1, not_a1, q1, not_q1]
    Y_address = [q1, not_q1, i1, i2]
    
    
    dY_johnson = one_bit_model(Y_johnson, T, params_johnson)
    dY_addr = one_bit_simple_addressing_ode_model(Y_address, T, params_addr)

    dY = np.append(dY_johnson, dY_addr)
    return dY


# TOP MODEL OF PROCESSOR WITH TWO BIT ADDRESSING
def two_bit_processor_ext(Y, T, params_johnson, params_addr):
    a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, i1, i2, i3, i4  = Y

    Y_johnson = [a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2]
    Y_address = [q1, not_q1, q2, not_q2, i1, i2, i3, i4]
    
    
    dY_johnson = two_bit_model(Y_johnson, T, params_johnson)
    dY_addr = two_bit_simple_addressing_ode_model(Y_address, T, params_addr)

    dY = np.append(dY_johnson, dY_addr)
    return dY

# TOP MODEL OF PROCESSOR WITH THREE BIT ADDRESSING
def three_bit_processor_ext(Y, T, params_johnson, params_addr):
    a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, a3, not_a3, q3, not_q3, i1, i2, i3, i4, i5, i6  = Y

    Y_johnson = [a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, a3, not_a3, q3, not_q3]
    Y_address = [q1, not_q1, q2, not_q2, q3, not_q3, i1, i2, i3, i4, i5, i6]
    
    
    dY_johnson = three_bit_model(Y_johnson, T, params_johnson)
    dY_addr = three_bit_simple_addressing_ode_model(Y_address, T, params_addr)

    dY = np.append(dY_johnson, dY_addr)
    return dY

# TOP MODEL OF PROCESSOR WITH FOUR BIT ADDRESSING
def four_bit_processor_ext(Y, T, params_johnson, params_addr):
    a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, a3, not_a3, q3, not_q3, a4, not_a4, q4, not_q4, i1, i2, i3, i4, i5, i6, i7, i8  = Y

    Y_johnson = [a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, a3, not_a3, q3, not_q3, a4, not_a4, q4, not_q4]
    Y_address = [q1, not_q1, q2, not_q2, q3, not_q3, q4, not_q4, i1, i2, i3, i4, i5, i6, i7, i8]
    
    
    dY_johnson = four_bit_model(Y_johnson, T, params_johnson)
    dY_addr = four_bit_simple_addressing_ode_model(Y_address, T, params_addr)

    dY = np.append(dY_johnson, dY_addr)
    return dY


# TOP MODEL OF PROCESSOR WITH FIVE BIT ADDRESSING
def five_bit_processor_ext(Y, T, params_johnson, params_addr):
    a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, a3, not_a3, q3, not_q3, a4, not_a4, q4, not_q4, a5, not_a5, q5, not_q5, i1, i2, i3, i4, i5, i6, i7, i8 ,i9, i10 = Y

    Y_johnson = [a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, a3, not_a3, q3, not_q3, a4, not_a4, q4, not_q4, a5, not_a5, q5, not_q5]
    Y_address = [q1, not_q1, q2, not_q2, q3, not_q3, q4, not_q4, q5, not_q5, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10]
    
    
    dY_johnson = five_bit_model(Y_johnson, T, params_johnson)
    dY_addr = five_bit_simple_addressing_ode_model(Y_address, T, params_addr)

    dY = np.append(dY_johnson, dY_addr)
    return dY




"""
PROCESSOR MODEL WITH EXTERNAL CLOCK AND RS inputs
external clock is required, more robust
jumps allowed
dodano 23. 1. 2020
"""

# TOP MODEL OF PROCESSOR WITH ONE BIT ADDRESSING AND FLIP-FLOP WITH RS ASYNCHRONOUS INPUTS
def one_bit_processor_ext_RS(Y, T, params_johnson_RS, params_addr):
    a1, not_a1, q1, not_q1, i1, i2  = Y

    R1 = 0
    S1 = 0

    Y_johnson = [a1, not_a1, q1, not_q1, R1, S1]
    Y_address = [q1, not_q1, i1, i2]
    
    
    dY_johnson = one_bit_model_RS(Y_johnson, T, params_johnson_RS)
    dY_addr = one_bit_simple_addressing_ode_model(Y_address, T, params_addr)

    dY = np.append(dY_johnson, dY_addr)
    return dY


# TOP MODEL OF PROCESSOR WITH TWO BIT ADDRESSING AND FLIP-FLOPS WITH RS ASYNCHRONOUS INPUTS
def two_bit_processor_ext_RS(Y, T, params_johnson_RS, params_addr):
    a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, i1, i2, i3, i4  = Y

    R1 = 0
    S1 = 0
    R2 = 0
    S2 = 0

    Y_johnson = [a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, R1, S1, R2, S2]
    Y_address = [q1, not_q1, q2, not_q2, i1, i2, i3, i4]
    
    
    dY_johnson = two_bit_model_RS(Y_johnson, T, params_johnson_RS)
    dY_addr = two_bit_simple_addressing_ode_model(Y_address, T, params_addr)

    dY = np.append(dY_johnson, dY_addr)
    return dY

# TOP MODEL OF PROCESSOR WITH THREE BIT ADDRESSING AND FLIP-FLOPS WITH RS ASYNCHRONOUS INPUTS
def three_bit_processor_ext_RS(Y, T, params_johnson_RS, params_addr, jump_src, jump_dst, i_src, i_dst):
    a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, a3, not_a3, q3, not_q3, i1, i2, i3, i4, i5, i6  = Y

    i_src = eval(i_src)
    
    R = [0,0,0]
    S = [0,0,0]

    for i in range(len(jump_src)):
        if jump_src[i] > jump_dst[i]:
                R[i] = i_src
        elif jump_src[i] < jump_dst[i]:
                S[i] = i_src

    R1, R2, R3 = R if T > 1 else [100,100,100]
    S1, S2, S3 = S

    Y_johnson = [a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, a3, not_a3, q3, not_q3, R1, S1, R2, S2, R3, S3]
    Y_address = [q1, not_q1, q2, not_q2, q3, not_q3, i1, i2, i3, i4, i5, i6]
    
    
    dY_johnson = three_bit_model_RS(Y_johnson, T, params_johnson_RS)
    dY_addr = three_bit_simple_addressing_ode_model(Y_address, T, params_addr)

    dY = np.append(dY_johnson, dY_addr)
    return dY

"""
PROCESSOR MODEL WITH EXTERNAL CLOCK AND RS inputs AND JUMP CONDITIONS
dodano 24. 1. 2020
"""
def get_condition(x0, delta, t):
        return x0 * np.e**(-delta*t)

# TOP MODEL OF PROCESSOR WITH THREE BIT ADDRESSING AND CONDITIONAL JUMPS
def three_bit_processor_ext_RS_cond(Y, T, params_johnson_RS, params_addr, jump_src, jump_dst, i_src, i_dst, condition):
    a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, a3, not_a3, q3, not_q3, i1, i2, i3, i4, i5, i6  = Y

    x0_cond, delta_cond, KD_cond, condition_type = condition
    cond = get_condition(x0_cond, delta_cond, T)

    

    i_src = eval(i_src)   
    
    R = np.array([0,0,0])
    S = np.array([0,0,0])

    for i in range(len(jump_src)):
        if jump_src[i] > jump_dst[i]:
                R[i] = i_src
        elif jump_src[i] < jump_dst[i]:
                S[i] = i_src
    if condition_type == "induction":
        R = induction(R, cond, KD_cond)
        S = induction(S, cond, KD_cond)
    else:
        R = inhibition(R, cond, KD_cond)
        S = inhibition(S, cond, KD_cond)


    R1, R2, R3 = R
    S1, S2, S3 = S

    Y_johnson = [a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, a3, not_a3, q3, not_q3, R1, S1, R2, S2, R3, S3]
    Y_address = [q1, not_q1, q2, not_q2, q3, not_q3, i1, i2, i3, i4, i5, i6]
    
    
    dY_johnson = three_bit_model_RS(Y_johnson, T, params_johnson_RS)
    dY_addr = three_bit_simple_addressing_ode_model(Y_address, T, params_addr)

    dY = np.append(dY_johnson, dY_addr)
    return dY
