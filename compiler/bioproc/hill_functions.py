import numpy as np

def repress_1(R, Kd, n):
    return 1/(1 + pow(R/Kd, n))

def repress_2(R1, R2, Kd1, n1, Kd2=0, n2=0, comp = 0):

    if not Kd2:
        Kd2 = Kd1
    if not n2:
        n2 = n1

    if comp:
        return 1/(1 + pow(R1/Kd1, n1) + pow(R2/Kd2, n2))        
    else:
        return 1/(1 + pow(R1/Kd1, n1) + pow(R2/Kd2, n2) + pow(R1/Kd1, n1) * pow(R2/Kd2, n2))


def activate_1(A, Kd, n):
    return pow(A/Kd, n)/(1 + pow(A/Kd, n))

def activate_2(A1, A2, Kd1, n1, Kd2=0, n2=0):

    if not Kd2:
        Kd2 = Kd1
    if not n2:
        n2 = n1
     
    return pow(A1/Kd1, n1) * pow(A2/Kd2, n2)/(1 + pow(A1/Kd1, n1) + pow(A2/Kd2, n2) + pow(A1/Kd1, n1) * pow(A2/Kd2, n2))
  
def activate_3(A1, A2, A3, Kd, n):
    return activate_1(A1, Kd, n) * activate_1(A2, Kd, n) * activate_1(A3, Kd, n)


def hybrid(A, R, Kd_A, n_A, Kd_R, n_R):
    return activate_1(A, Kd_A, n_A) * repress_1(R, Kd_R, n_R)
    
  
def hybrid_AAR(A1, A2, A3, Kd_A, n_A, Kd_R=0, n_R=0):
    if not Kd_R:
        Kd_R = Kd_A
    if not n_R:
        n_R = n_A
    return activate_1(A1, Kd_A, n_A) * activate_1(A2, Kd_A, n_A) * repress_1(A3, Kd_R, n_R)
    
def induction(x, i, Kd, n=1):
    # Michelis-Menten equation - see Alon p. 244, equation A.2.4
    # induction of protease by external condition / inducer - conditional jumps
    return x*i**n/(i**n+Kd**n)   

def inhibition(x, i, Kd, n=1):
    # Michelis-Menten equation - see Alon p. 245, equation A.2.5
    # inhibition of protease by external condition / inducer - conditional jumps
     return x/(1+(i/Kd)**n)
     

# CLOCK GENERATOR 
def get_clock(t, amp=100, per=24, phase = 0):                   
    return amp*(np.sin(2*np.pi*(t)/per + phase) + 1)/2    