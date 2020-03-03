from bioproc.hill_functions import *
def model(Y, T, params):
	a0, not_a0, q0, not_q0,a1, not_a1, q1, not_q1,a2, not_a2, q2, not_q2,a3, not_a3, q3, not_q3,i1, i2, i3, i4, i5, i6, i7, i8,A,B,C,=Y
	alpha1, alpha2, alpha3, alpha4, delta1, delta2, Kd, n, deltaE, KM,alpha_a, delta_a, Kd_a, n_a,KD_cond,=params
	clk = get_clock(T)
	prog_Kd_B=10
	prog_delta_B=0.1
	prog_alpha_C=10
	prog_delta_C=0.1
	prog_n_B=2
	prog_delta_A=0.1
	prog_alpha_B=10
	prog_n_A=2
	prog_n_C=2
	prog_alpha_A=10
	prog_Kd_A=10
	prog_Kd_C=10
	dA_dt=+prog_alpha_A*activate_1(i1,prog_Kd_A,prog_n_A)-prog_delta_A*A
	dB_dt=+prog_alpha_B*activate_1(i2,prog_Kd_B,prog_n_B)-prog_delta_B*B
	dC_dt=+prog_alpha_C*activate_2(i3,A,prog_Kd_C,prog_n_C)+prog_alpha_C*activate_2(i3,B,prog_Kd_C,prog_n_C)-prog_delta_C*C
	RESET0=0 if T > 1 else 100
	SET0=0 
	RESET1=0 if T > 1 else 100
	SET1=0 
	RESET2=0 if T > 1 else 100
	SET2=0 
	RESET3=0 if T > 1 else 100
	SET3=0 
	d0=not_q3
	sum_one = a0 + q0
	sum_zero = not_a0 + not_q0
	da0_dt     = alpha1*(pow(d0/Kd, n)/(1 + pow(d0/Kd, n) + pow(clk/Kd, n) + pow(d0/Kd, n)*pow(clk/Kd, n))) + alpha2*(1/(1 + pow(not_a0/Kd, n))) - delta1 *a0 -a0*(deltaE*RESET0/(KM+sum_one))
	dnot_a0_dt = alpha1*(1/(1 + pow(d0/Kd, n) + pow(clk/Kd, n) + pow(d0/Kd, n)*pow(clk/Kd, n))) + alpha2*(1/(1 + pow(a0/Kd, n))) - delta1*not_a0-not_a0*(deltaE*SET0/(KM+sum_zero))
	dq0_dt     = alpha3*((pow(a0/Kd, n)*pow(clk/Kd, n))/(1 + pow(a0/Kd, n) + pow(clk/Kd, n) + pow(a0/Kd, n)*pow(clk/Kd, n))) + alpha4*(1/(1 + pow(not_q0/Kd, n))) - delta2*q0-q0*(deltaE*RESET0/(KM+sum_one))
	dnot_q0_dt = alpha3*((pow(not_a0/Kd, n)*pow(clk/Kd, n))/(1 + pow(not_a0/Kd, n) + pow(clk/Kd, n) + pow(not_a0/Kd, n)*pow(clk/Kd, n))) + alpha4*(1/(1 + pow(q0/Kd, n))) - delta2*not_q0 -not_q0*(deltaE*SET0/(KM+sum_zero))
	d1=q0
	sum_one = a1 + q1
	sum_zero = not_a1 + not_q1
	da1_dt     = alpha1*(pow(d1/Kd, n)/(1 + pow(d1/Kd, n) + pow(clk/Kd, n) + pow(d1/Kd, n)*pow(clk/Kd, n))) + alpha2*(1/(1 + pow(not_a1/Kd, n))) - delta1 *a1 -a1*(deltaE*RESET1/(KM+sum_one))
	dnot_a1_dt = alpha1*(1/(1 + pow(d1/Kd, n) + pow(clk/Kd, n) + pow(d1/Kd, n)*pow(clk/Kd, n))) + alpha2*(1/(1 + pow(a1/Kd, n))) - delta1*not_a1-not_a1*(deltaE*SET1/(KM+sum_zero))
	dq1_dt     = alpha3*((pow(a1/Kd, n)*pow(clk/Kd, n))/(1 + pow(a1/Kd, n) + pow(clk/Kd, n) + pow(a1/Kd, n)*pow(clk/Kd, n))) + alpha4*(1/(1 + pow(not_q1/Kd, n))) - delta2*q1-q1*(deltaE*RESET1/(KM+sum_one))
	dnot_q1_dt = alpha3*((pow(not_a1/Kd, n)*pow(clk/Kd, n))/(1 + pow(not_a1/Kd, n) + pow(clk/Kd, n) + pow(not_a1/Kd, n)*pow(clk/Kd, n))) + alpha4*(1/(1 + pow(q1/Kd, n))) - delta2*not_q1 -not_q1*(deltaE*SET1/(KM+sum_zero))
	d2=q1
	sum_one = a2 + q2
	sum_zero = not_a2 + not_q2
	da2_dt     = alpha1*(pow(d2/Kd, n)/(1 + pow(d2/Kd, n) + pow(clk/Kd, n) + pow(d2/Kd, n)*pow(clk/Kd, n))) + alpha2*(1/(1 + pow(not_a2/Kd, n))) - delta1 *a2 -a2*(deltaE*RESET2/(KM+sum_one))
	dnot_a2_dt = alpha1*(1/(1 + pow(d2/Kd, n) + pow(clk/Kd, n) + pow(d2/Kd, n)*pow(clk/Kd, n))) + alpha2*(1/(1 + pow(a2/Kd, n))) - delta1*not_a2-not_a2*(deltaE*SET2/(KM+sum_zero))
	dq2_dt     = alpha3*((pow(a2/Kd, n)*pow(clk/Kd, n))/(1 + pow(a2/Kd, n) + pow(clk/Kd, n) + pow(a2/Kd, n)*pow(clk/Kd, n))) + alpha4*(1/(1 + pow(not_q2/Kd, n))) - delta2*q2-q2*(deltaE*RESET2/(KM+sum_one))
	dnot_q2_dt = alpha3*((pow(not_a2/Kd, n)*pow(clk/Kd, n))/(1 + pow(not_a2/Kd, n) + pow(clk/Kd, n) + pow(not_a2/Kd, n)*pow(clk/Kd, n))) + alpha4*(1/(1 + pow(q2/Kd, n))) - delta2*not_q2 -not_q2*(deltaE*SET2/(KM+sum_zero))
	d3=q2
	sum_one = a3 + q3
	sum_zero = not_a3 + not_q3
	da3_dt     = alpha1*(pow(d3/Kd, n)/(1 + pow(d3/Kd, n) + pow(clk/Kd, n) + pow(d3/Kd, n)*pow(clk/Kd, n))) + alpha2*(1/(1 + pow(not_a3/Kd, n))) - delta1 *a3 -a3*(deltaE*RESET3/(KM+sum_one))
	dnot_a3_dt = alpha1*(1/(1 + pow(d3/Kd, n) + pow(clk/Kd, n) + pow(d3/Kd, n)*pow(clk/Kd, n))) + alpha2*(1/(1 + pow(a3/Kd, n))) - delta1*not_a3-not_a3*(deltaE*SET3/(KM+sum_zero))
	dq3_dt     = alpha3*((pow(a3/Kd, n)*pow(clk/Kd, n))/(1 + pow(a3/Kd, n) + pow(clk/Kd, n) + pow(a3/Kd, n)*pow(clk/Kd, n))) + alpha4*(1/(1 + pow(not_q3/Kd, n))) - delta2*q3-q3*(deltaE*RESET3/(KM+sum_one))
	dnot_q3_dt = alpha3*((pow(not_a3/Kd, n)*pow(clk/Kd, n))/(1 + pow(not_a3/Kd, n) + pow(clk/Kd, n) + pow(not_a3/Kd, n)*pow(clk/Kd, n))) + alpha4*(1/(1 + pow(q3/Kd, n))) - delta2*not_q3 -not_q3*(deltaE*SET3/(KM+sum_zero))
	di1_dt = alpha_a * activate_2(not_q0, not_q3, Kd_a, n_a) - delta_a * i1
	di2_dt = alpha_a * activate_2(q0, not_q1, Kd_a, n_a) - delta_a * i2
	di3_dt = alpha_a * activate_2(q1, not_q2, Kd_a, n_a) - delta_a * i3
	di4_dt = alpha_a * activate_2(q2, not_q3, Kd_a, n_a) - delta_a * i4
	di5_dt = alpha_a * activate_2(q0, q3, Kd_a, n_a) - delta_a * i5
	di6_dt = alpha_a * activate_2(not_q0, q1, Kd_a, n_a) - delta_a * i6
	di7_dt = alpha_a * activate_2(not_q1, q2, Kd_a, n_a) - delta_a * i7
	di8_dt = alpha_a * activate_2(not_q2, q3, Kd_a, n_a) - delta_a * i8
	return [da0_dt,dnot_a0_dt,dq0_dt,dnot_q0_dt,da1_dt,dnot_a1_dt,dq1_dt,dnot_q1_dt,da2_dt,dnot_a2_dt,dq2_dt,dnot_q2_dt,da3_dt,dnot_a3_dt,dq3_dt,dnot_q3_dt,di1_dt,di2_dt,di3_dt,di4_dt,di5_dt,di6_dt,di7_dt,di8_dt,dA_dt,dB_dt,dC_dt]