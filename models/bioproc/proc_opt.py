import numpy as np
import math
import peakutils 
import numpy.fft as fft
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from bioproc.hill_functions import *  

class BioProc:
	def __init__(self, parameter_values, params, initial_conditions, avg_dev=10, fake_clock=False):    
		self.nParams = len(params)  
		self.params = params #model parameters		
		self.parameter_values = parameter_values #allowed parameter ranges  
		self.y0 = initial_conditions   
		self.modes = [self.eval]   		        
		self.fake_clock = fake_clock #if true synchronization signal clk is external, else clk is internal    
		self.dt = 0.001
		self.dt_eval = 1           		
		self.T =  192 #hours                 
		self.N = int(self.T/self.dt) 
		self.N_eval = int(self.T/self.dt_eval)   
		self.jump = int(self.dt_eval/self.dt)     		   
		self.ts = np.linspace(0, self.T, self.N)  
		self.ts_eval = np.linspace(0, self.T, self.N_eval) 		
		self.threshold = -(avg_dev*avg_dev)*self.N_eval  #20 nM from ideal signal on average                                  
		s_width = int(self.N_eval/8.0)       
		self.amp = 50 
		a = 0.5   
		b = 12    
		c =	0.5	
		#print(self.N)
		#print(self.N_eval)  
		snId = (np.sin(2*np.pi*np.linspace(-0.25, 0.75, s_width)) + 1)*self.amp
		snFF = 1.0/(1 + pow(abs((np.linspace(0, 1, 2*s_width) - c)/a), 2*b))*2*self.amp     	
		
		self.id1 = np.array([0]*self.N_eval, dtype=np.float64) 
		self.id2 = np.array([0]*self.N_eval, dtype=np.float64)  
		self.id3 = np.array([0]*self.N_eval, dtype=np.float64)  
		self.id4 = np.array([0]*self.N_eval, dtype=np.float64) 
		
		self.idFF1 = np.array([0]*self.N_eval, dtype=np.float64) 
		self.idFF2 = np.array([0]*self.N_eval, dtype=np.float64)  

		self.idFF1[1*s_width:3*s_width] = snFF  
		self.idFF1[5*s_width:7*s_width] = snFF   

		self.idFF2[2*s_width:4*s_width] = snFF     
		self.idFF2[6*s_width:8*s_width] = snFF   			
 		
		self.id1[0:s_width] = snId
		self.id1[4*s_width:5*s_width] = snId			

		self.id2[1*s_width:2*s_width] = snId
		self.id2[5*s_width:6*s_width] = snId		

		self.id3[2*s_width:3*s_width] = snId 
		self.id3[6*s_width:7*s_width] = snId 			

		self.id4[3*s_width:4*s_width] = snId
		self.id4[7*s_width:8*s_width] = snId   		
		 
		print(self.threshold) 
	
		"""
		CLK = get_clock(self.ts)
		plt.plot(self.ts, CLK)  
		plt.plot(self.ts_eval, self.id1, "r")
		plt.plot(self.ts_eval, self.id2, "g")
		plt.plot(self.ts_eval, self.id3, "b")		 
		plt.plot(self.ts_eval, self.id4, "o") 
		plt.plot(self.ts_eval, self.idFF1)     
		plt.plot(self.ts_eval, self.idFF2)     
		plt.show()      
		"""     
		
	def getFitness(self, signal, ideal):  
		#print(signal.shape) 
		diff = signal[0::self.jump] - ideal        
		fitness = -np.dot(diff, diff)  
		return fitness     
	
	def eval(self, candidate):    
		Y = np.array(self.simulate(candidate)) 
		i1 = Y[:,-4]           
		i2 = Y[:,-3]          
		i3 = Y[:,-2]    
		i4 = Y[:,-1]       		
		
		ins = [i1, i2, i3, i4]          
		ideals = [self.id1, self.id2, self.id3, self.id4]  
		weights = [0.25, 0.25, 0.25, 0.25]       		
				
		t_fitness = 0      
		for i, ideal, w in zip(ins, ideals, weights):       
			t_fitness += w*self.getFitness(i, ideal) 
			
			
		weights_flip_flop = [0.5, 0.5] 
		idealFlipFlops = [self.idFF1, self.idFF2]   
		ffs = [Y[:,2], Y[:,6]] #q1 and q2 signals 
		for ff, idealFlipFlop, w in zip(ffs, idealFlipFlops, weights_flip_flop):
			t_fitness += w*self.getFitness(ff, idealFlipFlop)           
			
		t_fitness = t_fitness/2.0     
		
		print(t_fitness) 
		return t_fitness,               		
		
	def isViable(self, point):  
		fitness = self.eval(point)   
		#print(fitness[0])   
		return fitness[0] >= self.threshold          		
	
	def simulate(self, candidate):   
		params_ff = candidate[0:8]        
		params_addr = candidate[8:]          
		
		Y = odeint(two_bit_model_full_with_instructions, self.y0, self.ts, args=(params_ff, params_addr, self.fake_clock)) 
        # MORALO BI BITI
        #Y = odeint(two_bit_processor_ext, self.y0, self.ts, args=(params_ff, params_addr)) 
        
		#print("simulating")    
		return Y    