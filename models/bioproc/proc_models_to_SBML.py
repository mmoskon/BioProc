import numpy as np
import simplesbml #https://simplesbml.readthedocs.io/en/latest/
import os
 
"""
FLIP-FLOP MODELS
"""
	
# MASTER-SLAVE D FLIP-FLOP MODEL

def ff_ode_model(Y, params): 
    
    species = {"$a$": Y[0],
               "$not_a$": Y[1],
               "$q$": Y[2],
               "$not_q$": Y[3], 
               "$d$": Y[4],
               "$clk$": Y[5]}

    parameters = {"$alpha1$": params[0], 
                  "$alpha2$": params[1],
                  "$alpha3$": params[2],
                  "$alpha4$": params[3],
                  "$delta1$": params[4],
                  "$delta2$": params[5],
                  "$Kd$": params[6],
                  "$n$": params[7]}

    #a, not_a, q, not_q, d, clk = Y
    #alpha1, alpha2, alpha3, alpha4, delta1, delta2, Kd, n = params

    da_dt     = "$alpha1$*(pow($d$/$Kd$, $n$)/(1 + pow($d$/$Kd$, $n$) + pow($clk$/$Kd$, $n$) + pow($d$/$Kd$, $n$)*pow($clk$/$Kd$, $n$))) + $alpha2$*(1/(1 + pow($not_a$/$Kd$, $n$))) - $delta1$ *$a$"
    dnot_a_dt = "$alpha1$*(1/(1 + pow($d$/$Kd$, $n$) + pow($clk$/$Kd$, $n$) + pow($d$/$Kd$, $n$)*pow($clk$/$Kd$, $n$))) + $alpha2$*(1/(1 + pow($a$/$Kd$, $n$))) - $delta1$*$not_a$"
    dq_dt     = "$alpha3$*((pow($a$/$Kd$, $n$)*pow($clk$/$Kd$, $n$))/(1 + pow($a$/$Kd$, $n$) + pow($clk$/$Kd$, $n$) + pow($a$/$Kd$, $n$)*pow($clk$/$Kd$, $n$))) + $alpha4$*(1/(1 + pow($not_q$/$Kd$, $n$))) - $delta2$*$q$"  
    dnot_q_dt = "$alpha3$*((pow($not_a$/$Kd$, $n$)*pow($clk$/$Kd$, $n$))/(1 + pow($not_a$/$Kd$, $n$) + pow($clk$/$Kd$, $n$) + pow($not_a$/$Kd$, $n$)*pow($clk$/$Kd$, $n$))) + $alpha4$*(1/(1 + pow($q$/$Kd$, $n$))) - $delta2$*$not_q$"   

    for old_s, new_s in species.items():
        da_dt = da_dt.replace(old_s, new_s)
        dnot_a_dt = dnot_a_dt.replace(old_s, new_s)
        dq_dt = dq_dt.replace(old_s, new_s)
        dnot_q_dt = dnot_q_dt.replace(old_s, new_s)
    
    for old_s, new_s in parameters.items():
        da_dt = da_dt.replace(old_s, new_s)
        dnot_a_dt = dnot_a_dt.replace(old_s, new_s)
        dq_dt = dq_dt.replace(old_s, new_s)
        dnot_q_dt = dnot_q_dt.replace(old_s, new_s)

    return [da_dt, dnot_a_dt, dq_dt, dnot_q_dt]

#f = open("test.txt", "w")
#X = ff_ode_model(["a", "not_a", "q", "not_q", "d", "clk"], ["alpha1", "alpha2", "alpha3", "alpha4", "delta1", "delta2", "Kd", "n"])
#for x in X:
#    print(x, file=f)
#f.close()

# MASTER-SLAVE D FLIP-FLOP MODEL with external set and reset inputs
def ff_ode_model_RS(Y, params): 
        
    species = {"$a$": Y[0],
               "$not_a$": Y[1],
               "$q$": Y[2],
               "$not_q$": Y[3], 
               "$d$": Y[4],
               "$clk$": Y[5],
               "$RESET$": Y[6],
               "$SET$": Y[7]}

    parameters = {"$alpha1$": params[0], 
                  "$alpha2$": params[1],
                  "$alpha3$": params[2],
                  "$alpha4$": params[3],
                  "$delta1$": params[4],
                  "$delta2$": params[5],
                  "$Kd$": params[6],
                  "$n$": params[7], 
                  "$deltaE$": params[8],
                  "$KM$": params[9]}

    #a, not_a, q, not_q, d, clk, RESET, SET = Y
    #alpha1, alpha2, alpha3, alpha4, delta1, delta2, Kd, n, deltaE, KM = params

    repress_both = True

    
    

    da_dt     = "$alpha1$*(pow($d$/$Kd$, $n$)/(1 + pow($d$/$Kd$, $n$) + pow($clk$/$Kd$, $n$) + pow($d$/$Kd$, $n$)*pow($clk$/$Kd$, $n$))) + $alpha2$*(1/(1 + pow($not_a$/$Kd$, $n$))) - $delta1$ * $a$"

    #deltaE = delta1
    if repress_both:
        da_dt += "-$a$*($deltaE$*$RESET$/($KM$+$sum_one$))"
    else:
        da_dt += "-a*(deltaE*RESET/(KM+a))"


    dnot_a_dt = "$alpha1$*(1/(1 + pow($d$/$Kd$, $n$) + pow($clk$/$Kd$, $n$) + pow($d$/$Kd$, $n$)*pow($clk$/$Kd$, $n$))) + $alpha2$*(1/(1 + pow($a$/$Kd$, $n$))) - $delta1$*$not_a$"
    if repress_both:
        dnot_a_dt += "-$not_a$*($deltaE$*$SET$/($KM$+$sum_zero$))"
    else:
        dnot_a_dt += "-$not_a$*($deltaE$*$SET$/($KM$+$not_a$))"

    #deltaE = delta2
    dq_dt     = "$alpha3$*((pow($a$/$Kd$, $n$)*pow($clk$/$Kd$, $n$))/(1 + pow($a$/$Kd$, $n$) + pow($clk$/$Kd$, $n$) + pow($a$/$Kd$, $n$)*pow($clk$/$Kd$, $n$))) + $alpha4$*(1/(1 + pow($not_q$/$Kd$, $n$))) - $delta2$*$q$"
    if repress_both:
        dq_dt += "-$q$*($deltaE$*$RESET$/($KM$+$sum_one$))"
    
    dnot_q_dt = "$alpha3$*((pow($not_a$/$Kd$, $n$)*pow($clk$/$Kd$, $n$))/(1 + pow($not_a$/$Kd$, $n$) + pow($clk$/$Kd$, $n$) + pow($not_a$/$Kd$, $n$)*pow($clk$/$Kd$, $n$))) + $alpha4$*(1/(1 + pow($q$/$Kd$, $n$))) - $delta2$*$not_q$"
    if repress_both:
        dnot_q_dt += "-$not_q$*($deltaE$*$SET$/($KM$+$sum_zero$))"
   

    if repress_both:                 
            species2 = {}
            species2["$sum_one$"] = "$a$ + $q$"
            species2["$sum_zero$"] = "$not_a$ + $not_q$"

            for old_s, new_s in species2.items():
                da_dt = da_dt.replace(old_s, new_s)
                dnot_a_dt = dnot_a_dt.replace(old_s, new_s)
                dq_dt = dq_dt.replace(old_s, new_s)
                dnot_q_dt = dnot_q_dt.replace(old_s, new_s)
    
    for old_s, new_s in species.items():
        da_dt = da_dt.replace(old_s, new_s)
        dnot_a_dt = dnot_a_dt.replace(old_s, new_s)
        dq_dt = dq_dt.replace(old_s, new_s)
        dnot_q_dt = dnot_q_dt.replace(old_s, new_s)
    
    for old_s, new_s in parameters.items():
        da_dt = da_dt.replace(old_s, new_s)
        dnot_a_dt = dnot_a_dt.replace(old_s, new_s)
        dq_dt = dq_dt.replace(old_s, new_s)
        dnot_q_dt = dnot_q_dt.replace(old_s, new_s)



    return [da_dt, dnot_a_dt, dq_dt, dnot_q_dt]

#f = open("test.txt", "w")
#print(ff_ode_model_RS(["a", "not_a", "q", "not_q", "d", "clk","RESET", "SET"], ["alpha1", "alpha2", "alpha3", "alpha4", "delta1", "delta2", "Kd", "n", "deltaE", "KM"]), file=f)
#f.close()




"""
ADRESSING MODELS
"""

def activate_1(A, Kd, n):
    d = {"$A$":A,
         "$Kd$":Kd,
         "$n$":n}
        
    expr = "pow($A$/$Kd$, $n$)/(1 + pow($A$/$Kd$, $n$))"
    for old_s, new_s in d.items():
        expr = expr.replace(old_s, new_s)

    return expr

def activate_2(A1, A2, Kd1, n1, Kd2="", n2=""):

    if not Kd2:
        Kd2 = Kd1
    if not n2:
        n2 = n1

    d = {"$A1$":A1,
         "$A2$":A2,
         "$Kd1$":Kd1,
         "$n1$":n1,
         "$Kd2$":Kd2,
         "$n2$":n2}
     
    expr = "pow($A1$/$Kd1$, $n1$) * pow($A2$/$Kd2$, $n2$)/(1 + pow($A1$/$Kd1$, $n1$) + pow($A2$/$Kd2$, $n2$) + pow($A1$/$Kd1$, $n1$) * pow($A2$/$Kd2$, $n2$))"
    
    for old_s, new_s in d.items():
        expr = expr.replace(old_s, new_s)

    return expr
    

# ONE BIT ADDRESSING MODEL SIMPLE
def one_bit_simple_addressing_ode_model(Y, params):
    #q1, not_q1, i1, i2 = Y
    #alpha, delta, Kd, n = params

    d = {"$q1$": Y[0], 
         "$not_q1$": Y[1], 
         "$i1$": Y[2], 
         "$i2$": Y[3],
         "$alpha$": params[0], 
         "$delta$": params[1], 
         "$Kd$": params[2], 
         "$n$": params[3]
        }

    di1_dt = "$alpha$ * " + activate_1("$not_q1$", "$Kd$", "$n$") + " - $delta$ * $i1$"
    di2_dt = "$alpha$ * " + activate_1("$q1$", "$Kd$", "$n$") + " - $delta$ * $i2$"
    
    for old_s, new_s in d.items():
        di1_dt = di1_dt.replace(old_s, new_s)
        di2_dt = di2_dt.replace(old_s, new_s)

    return [di1_dt, di2_dt]

#f = open("test.txt", "w")
#X = one_bit_simple_addressing_ode_model(["q1", "not_q1", "i1", "i2"], ["alpha", "delta", "Kd",  "n"])
#for x in X:
#    print(x, file=f)
#f.close()


   
# TWO BIT ADDRESSING MODEL SIMPLE
def two_bit_simple_addressing_ode_model(Y, params):
    #alpha, delta, Kd, n = params
    #q1, not_q1, q2, not_q2, i1, i2, i3, i4 = Y
    
    d = {"$q1$"     : Y[0], 
         "$not_q1$" : Y[1], 
         "$q2$"     : Y[2], 
         "$not_q2$" : Y[3], 
         "$i1$"     : Y[4], 
         "$i2$"     : Y[5],
         "$i3$"     : Y[6],
         "$i4$"     : Y[7],         
         "$alpha$"  : params[0], 
         "$delta$"  : params[1], 
         "$Kd$"     : params[2], 
         "$n$"      : params[3]
        }

    di1_dt = "$alpha$ * " + activate_2("$not_q1$", "$not_q2$", "$Kd$", "$n$") + " - $delta$ * $i1$"
    di2_dt = "$alpha$ * " + activate_2("$q1$", "$not_q2$", "$Kd$", "$n$") + "- $delta$ * $i2$"
    di3_dt = "$alpha$ * " + activate_2("$q1$", "$q2$", "$Kd$", "$n$") + "- $delta$ * $i3$"
    di4_dt = "$alpha$ * " + activate_2("$not_q1$", "$q2$", "$Kd$", "$n$") + "- $delta$ * $i4$"

    for old_s, new_s in d.items():
        di1_dt = di1_dt.replace(old_s, new_s)
        di2_dt = di2_dt.replace(old_s, new_s)
        di3_dt = di3_dt.replace(old_s, new_s)
        di4_dt = di4_dt.replace(old_s, new_s)
        
    return [di1_dt, di2_dt, di3_dt, di4_dt]

# THREE BIT ADDRESSING MODEL SIMPLE
def three_bit_simple_addressing_ode_model(Y, params):
    #alpha, delta, Kd, n = params    
    #q1, not_q1, q2, not_q2, q3, not_q3, i1, i2, i3, i4, i5, i6 = Y

    d = {"$q1$"     : Y[0], 
         "$not_q1$" : Y[1], 
         "$q2$"     : Y[2], 
         "$not_q2$" : Y[3], 
         "$q3$"     : Y[4], 
         "$not_q3$" : Y[5], 
         "$i1$"     : Y[6], 
         "$i2$"     : Y[7],
         "$i3$"     : Y[8],
         "$i4$"     : Y[9],
         "$i5$"     : Y[10],
         "$i6$"     : Y[11],         
         "$alpha$"  : params[0], 
         "$delta$"  : params[1], 
         "$Kd$"     : params[2], 
         "$n$"      : params[3]
        }


    di1_dt = "$alpha$ * " + activate_2("$not_q1$", "$not_q3$", "$Kd$", "$n$") + "- $delta$ * $i1$"
    di2_dt = "$alpha$ * " + activate_2("$q1$", "$not_q2$", "$Kd$", "$n$") + "- $delta$ * $i2$"
    di3_dt = "$alpha$ * " + activate_2("$q2$", "$not_q3$", "$Kd$", "$n$") + "- $delta$ * $i3$"
    di4_dt = "$alpha$ * " + activate_2("$q1$", "$q3$", "$Kd$", "$n$") + "- $delta$ * $i4$"
    di5_dt = "$alpha$ * " + activate_2("$not_q1$", "$q2$", "$Kd$", "$n$") + "- $delta$ * $i5$"
    di6_dt = "$alpha$ * " + activate_2("$not_q2$", "$q3$", "$Kd$", "$n$") + "- $delta$ * $i6$"

    for old_s, new_s in d.items():
        di1_dt = di1_dt.replace(old_s, new_s)
        di2_dt = di2_dt.replace(old_s, new_s)
        di3_dt = di3_dt.replace(old_s, new_s)
        di4_dt = di4_dt.replace(old_s, new_s)
        di5_dt = di5_dt.replace(old_s, new_s)
        di6_dt = di6_dt.replace(old_s, new_s)        

    return [di1_dt, di2_dt, di3_dt, di4_dt, di5_dt, di6_dt]

# FOUR BIT ADDRESSING MODEL SIMPLE
def four_bit_simple_addressing_ode_model(Y, params):
    #alpha, delta, Kd, n = params    
    #q1, not_q1, q2, not_q2, q3, not_q3, q4, not_q4, i1, i2, i3, i4, i5, i6, i7, i8 = Y

    d = {"$q1$"     : Y[0], 
         "$not_q1$" : Y[1], 
         "$q2$"     : Y[2], 
         "$not_q2$" : Y[3], 
         "$q3$"     : Y[4], 
         "$not_q3$" : Y[5], 
         "$q4$"     : Y[6], 
         "$not_q4$" : Y[7], 
         "$i1$"     : Y[8], 
         "$i2$"     : Y[9],
         "$i3$"     : Y[10],
         "$i4$"     : Y[11],
         "$i5$"     : Y[12],
         "$i6$"     : Y[13],
         "$i7$"     : Y[14],
         "$i8$"     : Y[15],
         "$alpha$"  : params[0], 
         "$delta$"  : params[1], 
         "$Kd$"     : params[2], 
         "$n$"      : params[3]
        }

    di1_dt = "$alpha$ * " + activate_2("$not_q1$", "$not_q4$", "$Kd$", "$n$") + "- $delta$ * $i1$"
    di2_dt = "$alpha$ * " + activate_2("$q1$", "$not_q2$", "$Kd$", "$n$") + " - $delta$ * $i2$"
    di3_dt = "$alpha$ * " + activate_2("$q2$", "$not_q3$", "$Kd$", "$n$") + "- $delta$ * $i3$"
    di4_dt = "$alpha$ * " + activate_2("$q3$", "$not_q4$", "$Kd$", "$n$") + "- $delta$ * $i4$"
    
    di5_dt = "$alpha$ * " + activate_2("$q1$", "$q4$", "$Kd$", "$n$") + "- $delta$ * $i5$"
    di6_dt = "$alpha$ * " + activate_2("$not_q1$", "$q2$", "$Kd$", "$n$") + "- $delta$ * $i6$"
    di7_dt = "$alpha$ * " + activate_2("$not_q2$", "$q3$", "$Kd$", "$n$") + "- $delta$ * $i7$"
    di8_dt = "$alpha$ * " + activate_2("$not_q3$", "$q4$", "$Kd$", "$n$") + "- $delta$ * $i8$"

    for old_s, new_s in d.items():
        di1_dt = di1_dt.replace(old_s, new_s)
        di2_dt = di2_dt.replace(old_s, new_s)
        di3_dt = di3_dt.replace(old_s, new_s)
        di4_dt = di4_dt.replace(old_s, new_s)
        di5_dt = di5_dt.replace(old_s, new_s)
        di6_dt = di6_dt.replace(old_s, new_s)
        di7_dt = di7_dt.replace(old_s, new_s)
        di8_dt = di8_dt.replace(old_s, new_s)

    return [di1_dt, di2_dt, di3_dt, di4_dt, di5_dt, di6_dt, di7_dt, di8_dt]



"""
JOHSON COUNTER MODELS 
"""
	
# TOP MODEL (JOHNSON): ONE BIT MODEL WITH EXTERNAL CLOCK
def one_bit_model(Y, params):
    a, not_a, q, not_q, clk= Y

    #clk = get_clock(T) 

    d = not_q
    Y_FF1 = [a, not_a, q, not_q, d, clk]

    dY = ff_ode_model(Y_FF1, params)

    return dY

# TOP MODEL (JOHNSON): TWO BIT MODEL WITH EXTERNAL CLOCK    
def two_bit_model(Y, params): 
    a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, clk = Y

    #clk = get_clock(T) 

    d1 = not_q2
    d2 = q1
    
    Y_FF1 = [a1, not_a1, q1, not_q1, d1, clk]
    Y_FF2 = [a2, not_a2, q2, not_q2, d2, clk]

    dY1 = ff_ode_model(Y_FF1, params)
    dY2 = ff_ode_model(Y_FF2, params)

    dY = dY1 + dY2

    return dY

# TOP MODEL (JOHNSON): THREE BIT MODEL WITH EXTERNAL CLOCK    
def three_bit_model(Y, params):
    a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, a3, not_a3, q3, not_q3, clk = Y

    #clk = get_clock(T) 

    d1 = not_q3
    d2 = q1
    d3 = q2
    
    Y_FF1 = [a1, not_a1, q1, not_q1, d1, clk]
    Y_FF2 = [a2, not_a2, q2, not_q2, d2, clk]
    Y_FF3 = [a3, not_a3, q3, not_q3, d3, clk]

    dY1 = ff_ode_model(Y_FF1, params)
    dY2 = ff_ode_model(Y_FF2, params)
    dY3 = ff_ode_model(Y_FF3, params)

    dY = dY1 + dY2 + dY3

    return dY

# TOP MODEL (JOHNSON): FOUR BIT MODEL WITH EXTERNAL CLOCK    
def four_bit_model(Y, params):
    a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, a3, not_a3, q3, not_q3, a4, not_a4, q4, not_q4, clk = Y

    #clk = get_clock(T) 

    d1 = not_q4
    d2 = q1
    d3 = q2
    d4 = q3

    Y_FF1 = [a1, not_a1, q1, not_q1, d1, clk]
    Y_FF2 = [a2, not_a2, q2, not_q2, d2, clk]
    Y_FF3 = [a3, not_a3, q3, not_q3, d3, clk]
    Y_FF4 = [a4, not_a4, q4, not_q4, d4, clk]

    dY1 = ff_ode_model(Y_FF1, params)
    dY2 = ff_ode_model(Y_FF2, params)
    dY3 = ff_ode_model(Y_FF3, params)
    dY4 = ff_ode_model(Y_FF4, params)

    dY = dY1 + dY2 + dY3 + dY4

    return dY


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

"""
PROCESSOR MODEL
"""

# TOP MODEL OF PROCESSOR WITH ONE BIT ADDRESSING
def one_bit_processor_ext(Y, params_johnson, params_addr):
    a1, not_a1, q1, not_q1, i1, i2, clk  = Y

    Y_johnson = [a1, not_a1, q1, not_q1, clk]
    Y_address = [q1, not_q1, i1, i2]
    
    
    dY_johnson = one_bit_model(Y_johnson, params_johnson)
    dY_addr = one_bit_simple_addressing_ode_model(Y_address, params_addr)

    dY = dY_johnson + dY_addr

    
    return dY

# TOP MODEL OF PROCESSOR WITH TWO BIT ADDRESSING
def two_bit_processor_ext(Y, params_johnson, params_addr):
    a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, i1, i2, i3, i4, clk  = Y

    Y_johnson = [a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, clk]
    Y_address = [q1, not_q1, q2, not_q2, i1, i2, i3, i4]
    
    
    dY_johnson = two_bit_model(Y_johnson, params_johnson)
    dY_addr = two_bit_simple_addressing_ode_model(Y_address, params_addr)

    dY = dY_johnson + dY_addr
    return dY

# TOP MODEL OF PROCESSOR WITH THREE BIT ADDRESSING
def three_bit_processor_ext(Y, params_johnson, params_addr):
    a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, a3, not_a3, q3, not_q3, i1, i2, i3, i4, i5, i6, clk  = Y

    Y_johnson = [a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, a3, not_a3, q3, not_q3, clk]
    Y_address = [q1, not_q1, q2, not_q2, q3, not_q3, i1, i2, i3, i4, i5, i6]
    
    
    dY_johnson = three_bit_model(Y_johnson, params_johnson)
    dY_addr = three_bit_simple_addressing_ode_model(Y_address, params_addr)

    dY = dY_johnson + dY_addr
    return dY

# TOP MODEL OF PROCESSOR WITH FOUR BIT ADDRESSING
def four_bit_processor_ext(Y, params_johnson, params_addr):
    a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, a3, not_a3, q3, not_q3, a4, not_a4, q4, not_q4, i1, i2, i3, i4, i5, i6, i7, i8, clk  = Y

    Y_johnson = [a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, a3, not_a3, q3, not_q3, a4, not_a4, q4, not_q4, clk]
    Y_address = [q1, not_q1, q2, not_q2, q3, not_q3, q4, not_q4, i1, i2, i3, i4, i5, i6, i7, i8]
    
    
    dY_johnson = four_bit_model(Y_johnson, params_johnson)
    dY_addr = four_bit_simple_addressing_ode_model(Y_address, params_addr)

    dY = dY_johnson + dY_addr
    return dY



def get_reaction(tmp, product):
    tmp=tmp.strip()    
    
    products_to_remove = set()
    product = product.replace("@","")
    products = {product}

    if tmp.startswith("-"):        
        products_to_remove = {product}
        tmp = tmp[1:]
    elif tmp.startswith("+"):
        tmp = tmp[1:]   

    

    reactants = set(tmp.split("@")[1::2])
 
    params = set(tmp.split("&")[1::2])
    
    
    products |= reactants
    products -= products_to_remove
    

    d = {'products': products,
         'reactants': reactants,
         'equation': tmp.replace("@","").replace("&",""),
         'params': params}

    return d

def to_reactions(Y,X):
    """
    reactions = {}

    
    for y,x in zip(Y,X):
        tmp = ""
        i = 0
        for s in y:                  
            if (s == "+" or s == "-") and (i == 0):
                if x in reactions:                         
                    reactions[x].append(tmp)
                else:
                    reactions[x]=[tmp]
                tmp = s
                continue
            elif s == "(":
                i += 1
            elif s == ")":
                i -= 1
            tmp += s
        if tmp: # if the end is reached
            if x in reactions:                         
                reactions[x].append(tmp)
            else:
                reactions[x]=[tmp]
    """
    reactions = []

    for y,x in zip(Y,X):
        tmp = ""
        i = 0
        for s in y:                  
            if (s == "+" or s == "-") and (i == 0):               
                reactions.append(get_reaction(tmp, x))
                tmp = s
                continue
            elif s == "(":
                i += 1
            elif s == ")":
                i -= 1
            tmp += s
        if tmp: # if the end is reached            
            reactions.append(get_reaction(tmp, x))
     
            
            

    return reactions

"""
f = open("test.txt", "w")

Y = ["@a1@", "@not_a1@", "@q1@", "@not_q1@", "@i1@", "@i2@", "@clk@"]
params_ff = ["&alpha1&", "&alpha2&", "&alpha3&", "&alpha4&", "&delta1&", "&delta2&", "&Kd&", "&n&"]
params_addr = ["&alpha_I&", "&delta_I&", "&Kd_I&", "&n_I&"]
#params_ff = ["$alpha1$", "$alpha2$", "$alpha3$", "$alpha4$", "$delta1$", "$delta2$", "$Kd$", "$n$"]
#params_addr = ["$alpha_I$", "$delta_I$", "$Kd_I$", "$n_I$"]

X = one_bit_processor_ext(Y, params_ff, params_addr)
for y,x in zip(Y,X):
    print(y + ": " + x, file=f)
f.close()

f = open("test2.txt", "w")
reacts = to_reactions(X,Y)
for r in reacts:
    print(r, file=f)
f.close()

#print(get_reaction("- &delta1& *@a1@", "@a1@"))
"""

def to_sbml(reacts, filename):
    model = simplesbml.SbmlModel()
    model.addCompartment(1, comp_id='comp')

    species = set()
    params = set()
    for react in reacts:
        species |= react['products']
        species |= react['reactants']
        params |= react['params']

    for s in species:
        model.addSpecies(s, 0, comp='comp')
    
    for p in params:
        model.addParameter(p, 0)

    for i,r in enumerate(reacts):
        model.addReaction(list(r['reactants']), list(r['products']), r['equation'], rxn_id=f"r{i}")
           

    f = open(filename, 'w')
    print(model.toSBML(), file=f)
    f.close()

def one_bit_processor_ext_to_sbml(filename):
    Y = ["@a1@", "@not_a1@", "@q1@", "@not_q1@", "@i1@", "@i2@", "@clk@"]
    params_ff = ["&alpha1&", "&alpha2&", "&alpha3&", "&alpha4&", "&delta1&", "&delta2&", "&Kd&", "&n&"]
    params_addr = ["&alpha_I&", "&delta_I&", "&Kd_I&", "&n_I&"]
    X = one_bit_processor_ext(Y, params_ff, params_addr)
    reacts = to_reactions(X,Y)

    to_sbml(reacts, filename)

def two_bit_processor_ext_to_sbml(filename):
    #Y = ["@a1@", "@not_a1@", "@q1@", "@not_q1@", "@i1@", "@i2@", "@clk@"]
    Y = ["@a1@", "@not_a1@", "@q1@", "@not_q1@", "@a2@", "@not_a2@", "@q2@", "@not_q2@", "@i1@", "@i2@", "@i3@", "@i4@", "@clk@"]
    params_ff = ["&alpha1&", "&alpha2&", "&alpha3&", "&alpha4&", "&delta1&", "&delta2&", "&Kd&", "&n&"]
    params_addr = ["&alpha_I&", "&delta_I&", "&Kd_I&", "&n_I&"]
    X = two_bit_processor_ext(Y, params_ff, params_addr)
    reacts = to_reactions(X,Y)

    to_sbml(reacts, filename)

def three_bit_processor_ext_to_sbml(filename):
    #Y = ["@a1@", "@not_a1@", "@q1@", "@not_q1@", "@i1@", "@i2@", "@clk@"]
    Y = ["@a1@", "@not_a1@", "@q1@", "@not_q1@", "@a2@", "@not_a2@", "@q2@", "@not_q2@", "@a3@", "@not_a3@", "@q3@", "@not_q3@", "@i1@", "@i2@", "@i3@", "@i4@", "@i5@", "@i6@", "@clk@"]
    params_ff = ["&alpha1&", "&alpha2&", "&alpha3&", "&alpha4&", "&delta1&", "&delta2&", "&Kd&", "&n&"]
    params_addr = ["&alpha_I&", "&delta_I&", "&Kd_I&", "&n_I&"]
    X = three_bit_processor_ext(Y, params_ff, params_addr)
    reacts = to_reactions(X,Y)

    #f = open("test.txt", "w")
    #for r in reacts:
    #    print(r, file=f)
    #f.close()

    to_sbml(reacts, filename)

def four_bit_processor_ext_to_sbml(filename):
    #Y = ["@a1@", "@not_a1@", "@q1@", "@not_q1@", "@i1@", "@i2@", "@clk@"]
    Y = ["@a1@", "@not_a1@", "@q1@", "@not_q1@", "@a2@", "@not_a2@", "@q2@", "@not_q2@", "@a3@", "@not_a3@", "@q3@", "@not_q3@", "@a4@", "@not_a4@", "@q4@", "@not_q4@", "@i1@", "@i2@", "@i3@", "@i4@", "@i5@", "@i6@", "@i7@", "@i8@", "@clk@"]
    params_ff = ["&alpha1&", "&alpha2&", "&alpha3&", "&alpha4&", "&delta1&", "&delta2&", "&Kd&", "&n&"]
    params_addr = ["&alpha_I&", "&delta_I&", "&Kd_I&", "&n_I&"]
    X = four_bit_processor_ext(Y, params_ff, params_addr)
    reacts = to_reactions(X,Y)

    to_sbml(reacts, filename)



one_bit_processor_ext_to_sbml(os.path.join("SBML", "one_bit_proc.sbml"))
two_bit_processor_ext_to_sbml(os.path.join("SBML", "two_bit_proc.sbml"))
three_bit_processor_ext_to_sbml(os.path.join("SBML", "three_bit_proc.sbml"))
four_bit_processor_ext_to_sbml(os.path.join("SBML", "four_bit_proc.sbml"))