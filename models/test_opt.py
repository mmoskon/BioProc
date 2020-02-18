from solver import Solver, Region
from numpy import random
import pickle
import os.path as path
import numpy as np 
from bioproc.proc_models import *   
from bioproc.proc_opt import *   
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

file =  os.path.join(".", "bioproc", "one_bit_model_new_new", "bioprocViableSet_IterGA.p")    
viablePoints = np.array(pickle.load(open(file, "rb")))     
number = np.size(viablePoints, 0) 
rndPoints = np.array(np.random.randint(number, size=20)) 
#print(viablePoints)   
print(viablePoints[rndPoints]) 
points = viablePoints[rndPoints]
print(points.shape)    
   	

model = BioProc(np.array(["protein_production", "protein_production", "protein_production", "protein_production", "protein_degradation", "protein_degradation", "Kd","hill", "protein_production", "protein_degradation", "Kd", "hill"]), model_mode=one_bit_processor_ext, avg_dev=30)          

for point in points:     
	print(point)  
	#params_addr = point  
	
	Y = model.simulate(point) 
	print(model.eval(point))   
 	
	#clk = Y[:,-5] 
	i1 = Y[:,-1]          
	i2 = Y[:,-2] 	 
	#i3 = Y[:,-3]         
	#i4 = Y[:,-4]  
	#i5 = Y[:,-5]          
	#i6 = Y[:,-6] 	    	
	
	q1 = Y[:,2]  
	#q2 = Y[:,6]   
	#q3 = Y[:,10]     	
	
	plt.plot(model.ts, i1, label='i1') 
	plt.plot(model.ts, i2, label='i2') 
	#plt.plot(model.ts, i3, label='i3') 
	#plt.plot(model.ts, i4, label='i4') 
	#plt.plot(model.ts, i5, label='i5') 
	#plt.plot(model.ts, i6, label='i6') 
	plt.plot(model.ts, q1, label='q1')  
	#plt.plot(model.ts, q2, label='q2') 
	#plt.plot(model.ts, q3, label='q3')  	
	
	plt.legend(loc="upper left")       
	plt.show()    
