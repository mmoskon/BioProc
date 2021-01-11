import numpy as np
import math
import peakutils 
import numpy.fft as fft
import matplotlib.pyplot as plt
import sys  
from scipy.integrate import odeint
from bioproc.hill_functions import *
from bioproc.proc_models import *     

def_parameter_values = {  "transcription": {"min": 0.01, "max": 50},   
            "translation": {"min": 0.01, "max": 50},  
            "protein_production": {"min": 0.1, "max": 50},                          
            "rna_degradation": {"min": 0.1, "max": 100},        
            "protein_degradation": {"min": 0.001, "max": 50},         
            "hill": {"min": 1, "max": 5},         
            "Kd": {"min": 0.01, "max": 250}, 
            "protease_concentration": {"min": 10, "max":1000}       
            }      

class BioProc:
    def __init__(self, params, model_mode, parameter_values=def_parameter_values, avg_dev=30, plot_fitness=False, plot_devs=True, load_init=True):       
        self.nParams = len(params)  
        self.params = params #model parameters      
        self.parameter_values = parameter_values #allowed parameter ranges  
        self.modes = [self.eval]                 
        self.dt = 0.001
        self.dt_eval = 1    
        self.model_mode = model_mode    
        self.multiplier = 1  
        self.omega = 1#nm^-1
    
        if self.model_mode == one_bit_processor_ext:
            self.multiplier = 1
            if load_init:
                self.y0 = np.loadtxt('init1bit.txt')
            else:                
                self.y0 = [0]*(self.multiplier*6) 
        elif self.model_mode == two_bit_processor_ext:
            self.multiplier = 2 
            if load_init:
                self.y0 = np.loadtxt('init2bit.txt')
            else:                
                self.y0 = [0]*(self.multiplier*6) 
        elif self.model_mode == three_bit_processor_ext:
            self.multiplier = 3         
            if load_init:
                self.y0 = np.loadtxt('init3bit.txt')
            else:                
                self.y0 = [0]*(self.multiplier*6) 
        elif self.model_mode == four_bit_processor_ext:
            self.multiplier = 4         
            if load_init:
                self.y0 = np.loadtxt('init4bit.txt')
            else:                
                self.y0 = [0]*(self.multiplier*6) 
        elif self.model_mode == five_bit_processor_ext:
            self.multiplier = 5         
            if load_init:
                self.y0 = np.loadtxt('init5bit.txt')
            else:                
                self.y0 = [0]*(self.multiplier*6) 
        else: 
            sys.exit('Error: unvalid model mode')       
            
                    
            
        self.T =  48*(2*self.multiplier) #hours * 2n                  
        self.N = int(self.T/self.dt) 
        self.N_eval = int(self.T/self.dt_eval)   
        self.jump = int(self.dt_eval/self.dt)              
        self.ts = np.linspace(0, self.T, self.N)  
        self.ts_eval = np.linspace(0, self.T, self.N_eval)      
                
        s_width = int(self.N_eval/(self.multiplier*4.0))          
        
        self.amp = 50 
        a = 0.5   
        b = 12    
        c = 0.5 

        snId = 1.0/(1 + pow(abs((np.linspace(0, 1, s_width) - c)/a), 2*b))*2*self.amp   #fitness bell shaped instructions #s_width*[self.amp*2]             
        snFF = 1.0/(1 + pow(abs((np.linspace(0, 1, self.multiplier*s_width) - c)/a), 2*b*self.multiplier))*2*self.amp       # self.multiplier*s_width*[self.amp*2]
        
        self.idealIns = []  
        self.idealsFF = []   
        
        for i in range(self.multiplier*2):      
            ins = np.array([0]*self.N_eval, dtype=np.float64) 
            ins[i*s_width:(i+1)*s_width] = snId 
            ins[(i + self.multiplier*2)*s_width:(i + self.multiplier*2 + 1)*s_width] = snId 
            self.idealIns.append(ins) 
            
        for i in range(self.multiplier): 
            idFF = np.array([0]*self.N_eval, dtype=np.float64)              
            idFF[(i+1)*s_width:(i+1+self.multiplier)*s_width] = snFF  
            idFF[(i+1 + self.multiplier*2)*s_width:(i+1+self.multiplier + self.multiplier*2)*s_width] = snFF                                                 
   
            self.idealsFF.append(idFF)              
                                
        self.threshold = -30 if self.model_mode == one_bit_processor_ext else -20 if self.model_mode == two_bit_processor_ext else -18 if self.model_mode == three_bit_processor_ext else -17 #self.getThreshold(avg_dev, plot_devs=plot_devs)                                             
        print(self.threshold) 
            
            
        if plot_fitness: 
            CLK = get_clock(self.ts)
            plt.plot(self.ts, CLK, alpha=0.5)  

            for ins in self.idealIns:
                plt.plot(self.ts_eval, ins, alpha=0.5)   
            for ff in self.idealsFF:
                plt.plot(self.ts_eval, ff, alpha=0.5)      
            plt.show()                  

    def getThreshold(self, avg_dev, plot_devs=False):       
        devIns = [] 
        devFFs = [] 

        for idealIn in self.idealIns:  
            loc = idealIn > np.max(idealIn) * 0.5
            devIn = idealIn.copy()
            devIn[loc] += avg_dev    
            #devIn += avg_dev    
            devIns.append(devIn)
            
        for idealFF in self.idealsFF:       
            loc = idealFF > np.max(idealFF) * 0.5
            devFF = idealFF.copy()
            devFF[loc] += avg_dev
            #devFF += avg_dev    
            devFFs.append(devFF)
                
        if plot_devs:  
            CLK = get_clock(self.ts)
            plt.plot(self.ts, CLK, alpha=0.5)  

            for ins in devIns: 
                plt.plot(self.ts_eval, ins, alpha=0.5)   
            for ff in devFFs:  
                plt.plot(self.ts_eval, ff, alpha=0.5)      
            plt.show() 
        return self.eval(None, devIns, devFFs)       

        
    def getTotalVolume(self):
        vol = 1.0
        for param in self.params:       
            vol = vol*(self.parameter_values[param]["max"] - self.parameter_values[param]["min"])
        return vol    
            
    def getFitness(self, signal, ideal):#, instruction = True):  
        if self.jump > len(signal):
            ideal = ideal[10:]
            signal = signal[10:]
        else:        
            ideal = ideal[10:]
            signal = signal[10*self.jump::self.jump]
        diff =  signal - ideal

         
        #diff = signal[0::self.jump] - ideal

        T1 = ideal > np.max(ideal)*0.9
        diff[T1] /= sum(T1)

        T0 = ideal <= np.max(ideal)*0.1
        diff[T0] /= sum(T0)
        #if instruction:
        #diff[T0] /= self.multiplier#*2
        #diff[T0] = 0
        
        T01 = ~(T0 | T1)
        diff[T01] = 0

            
        #fitness = -np.dot(diff, diff) 
        fitness = -np.sum(np.abs(diff))             
        #fitness /= np.exp(-1/4 * self.multiplier)
        
        #print((np.max(signal) - np.min(signal)/2 - self.amp))
        
        return fitness     
        
    def eval(self, candidate, ins=None, ffs=None):

        if ins == None and ffs == None:   
            Y = np.array(self.simulate(candidate)) 

            ins = []    
            for i in range(self.multiplier*2, 0,-1):   
                ins.append(Y[:,-i])
                #print(i)                

            ffs = []
            for i in range(self.multiplier):
                ff = Y[:,2 + i*4]  
                ffs.append(ff)          
            
        
        weight = 1.0/(self.multiplier*2)                        
        t_fitness = 0      
        
        for instruction, ideal in zip(ins, self.idealIns):       
            #print(len(instruction), len(ideal))
            t_fitness += weight*self.getFitness(instruction, ideal) 
            
        weight_flip_flop = 1.0/(self.multiplier)   
        
        for ff, idealFlipFlop in zip(ffs, self.idealsFF):  
            t_fitness += weight_flip_flop*self.getFitness(ff, idealFlipFlop)#, instruction=False)             
            
        t_fitness = t_fitness/2.0      
        
        print(t_fitness) 
        return t_fitness,                       
        
    def isViable(self, point, fitness=None): 
        if fitness == None:  
            fitness = self.eval(point)   
        #print(fitness[0])   
        return fitness[0] >= self.threshold       

    def getHotStart(self, params_ff, params_addr):
        y0 = self.y0  
        for i in range(3):                    
            ts_small = self.ts[self.ts < 5]
            Y = odeint(self.model_mode, y0, ts_small, args=(params_ff, params_addr))                     
            y0 = Y[-1,:].flatten()
        return y0
    
    def simulate(self, candidate, hot_start=True):   
        params_ff = candidate[0:8]        
        params_addr = candidate[8:]          
        
        y0 = self.y0       
        if hot_start:
            y0 = self.getHotStart(params_ff,params_addr)  
        
        Y = odeint(self.model_mode, y0, self.ts, args=(params_ff, params_addr))  
        """
        #Y = odeint(one_bit_processor_ext, self.y0, self.ts, args=(params_ff, params_addr))  
        #Y = odeint(two_bit_processor_ext, self.y0, self.ts, args=(params_ff, params_addr))  
        #Y = odeint(three_bit_processor_ext, self.y0, self.ts, args=(params_ff, params_addr))  
        """

        #print("simulating")     
        return Y

    def simulateStochastic(self, candidate,hot_start=True): 
        multiplier = self.multiplier   
        omega = self.omega 
        
                
        params_ff = candidate[0:8]          
        params_addr = candidate[8:] 
        
        y0 = self.y0 
        if hot_start:
            y0 = self.getHotStart(params_ff,params_addr)        
            
        y_conc = np.array(y0*omega).astype(int) 
        print(y_conc) 

        Y_total = []
        Y_total.append(y_conc)
        t = 0 
        T = [] 
        T.append(t) 
   

        N = np.zeros((6*multiplier, 16*multiplier)) #6*multiplier species, 16*multiplier reactions  
        for i in range(multiplier):
            #flip flops 
            #a
            N[i*4 + 0,i*12 + 0] = 1  
            N[i*4 + 0,i*12 + 1] = 1  
            N[i*4 + 0,i*12 + 2] = -1  
            #not a
            N[i*4 + 1,i*12 + 3] = 1  
            N[i*4 + 1,i*12 + 4] = 1  
            N[i*4 + 1,i*12 + 5] = -1  
            #q 
            N[i*4 + 2,i*12 + 6] = 1  
            N[i*4 + 2,i*12 + 7] = 1   
            N[i*4 + 2,i*12 + 8] = -1            
            #not q 
            N[i*4 + 3,i*12 + 9] = 1    
            N[i*4 + 3,i*12 + 10] = 1    
            N[i*4 + 3,i*12 + 11] = -1 
            #instructions               
            N[multiplier*4 + i*2+ 0, multiplier*12 + i*4 + 0] = 1  
            N[multiplier*4 + i*2 + 0, multiplier*12 + i*4 + 1] = -1     
            N[multiplier*4 + i*2 + 1, multiplier*12 + i*4 + 2] = 1  
            N[multiplier*4 + i*2 + 1, multiplier*12 + i*4 + 3] = -1             
            
        
        
        while t < self.T: 
            #choose two random numbers 
            r = np.random.uniform(size=2)
            r1 = r[0]
            r2 = r[1]           

            #get clk    
            clk = int(get_clock(t)*omega)       
            #get propensities               
            a = np.zeros(16*multiplier) 
            
            if self.model_mode == one_bit_processor_ext:
                ds = [y_conc[3]] #not q1   
                a[multiplier*12:] = addressing_stochastic_one_bit_model(y_conc, t, params_addr, omega)
                
            elif self.model_mode == two_bit_processor_ext:
                ds = [y_conc[7],y_conc[2]] #not q2, q1     
                a[multiplier*12:] = addressing_stochastic_two_bit_model(y_conc, t, params_addr, omega)
                
            else: #three_bit_processor_ext  
                ds = [y_conc[11],y_conc[2],y_conc[6]] #not q3, q1, q2       
                a[multiplier*12:] = addressing_stochastic_three_bit_model(y_conc, t, params_addr, omega) 

            for i in range(multiplier):
                y = y_conc[i*4:i*4+4] 
                y = np.append(y, ds[i]) #to do 
                y = np.append(y, clk) 
                a[i*12:i*12+12] = ff_stochastic_model(y, t, params_ff, omega)    
                    
            asum = np.cumsum(a)
            a0 = np.sum(a)  
            #get tau
            tau = (1.0/a0)*np.log(1.0/r1)    
        
            #print(t)       
            #select reaction 
            reaction_number = np.argwhere(asum > r2*a0)[0,0] #get first element         
        
            #update concentrations
            y_conc = y_conc + N[:,reaction_number]      
            Y_total.append(y_conc) 
            #update time
            t = t + tau  
            T.append(t)

            
        T = np.array(T)
        Y_total = np.array(Y_total) 
        return T, Y_total       