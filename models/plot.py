from solver import Solver, Region
from numpy import random
import pickle
import os.path as path
import numpy as np 
from bioproc.flip_flops import *   
import matplotlib.pyplot as plt 
from scipy.integrate import odeint 
import os.path 
from deap import creator, base, tools, algorithms    

#def get_clock(t, amp=100, per=24, phase = 0): 
#    return amp*(np.sin(2*np.pi*(t)/per + phase) + 1)/2

creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Candidate", list, fitness=creator.FitnessMax)		
toolbox = base.Toolbox()	 
toolbox.register("candidate", Solver.generateCandidate) 

file =  os.path.join(".", "bioproc", "external_clock2", "bioprocViableSet_IterGA.p")    
viablePoints = np.array(pickle.load(open(file, "rb")))     
number = np.size(viablePoints, 0) 
rndPoints = np.array(np.random.randint(number, size=20)) 
#print(viablePoints)   
print(viablePoints[rndPoints]) 
points = viablePoints[rndPoints]
print(points.shape) 


dt = 0.001 
T = 192 #hours 
N = int(T/dt)    
ts = np.linspace(0, T, N)   
y0 = np.array([0]*12)     

clk = get_clock(ts)

for point in points:  
	print(point)
	params_addr = point 
		
	params_ff = point[0:8]        
	params_addr = point[8:]   
	print(params_addr)
	
	Y = odeint(two_bit_model_full_with_instructions, y0, ts, args=(params_ff, params_addr, True ))       
	print(Y.shape) 
	#Y = np.split(Y, Y.shape[1], 1)
	
	#clk = Y[:,-5] 
	i1 = Y[:,-4]         
	i2 = Y[:,-3]        
	i3 = Y[:,-2]  
	i4 = Y[:,-1]   
	
	q1 = Y[:,2]       
	q2 = Y[:,6]        
	
	
	plt.plot(ts, i1, label='i1')
	plt.plot(ts, i2, label='i2')
	plt.plot(ts, i3, label='i3')
	plt.plot(ts, i4, label='i4')
	plt.plot(ts, q1, label='q1') 
	plt.plot(ts, q2, label='q2') 
	plt.plot(ts, clk, label='clk') 
	plt.legend()       
	plt.show()    