from bioproc.proc_opt import BioProc  
from bioproc.proc_models import *  
from bioproc.hill_functions import * 
import matplotlib.pyplot as plt 
from solver import * 
import numpy as np  
import os.path      
import pickle 
import time   

import seaborn as sns

sns.set_style("white")   

	
models = [one_bit_processor_ext, two_bit_processor_ext, three_bit_processor_ext]        
base_path = os.path.join(".", "bioproc")   
folders = [os.path.join(base_path, "01_lidija"), os.path.join(base_path, "02_ziga"), os.path.join(base_path, "03_miha")]                       
model_regions = []        
	
for model_index in range(len(models)):      				
	folder = folders[model_index]        		
	model = BioProc(np.array(["protein_production", "protein_production", "protein_production", "protein_production", "protein_degradation", "protein_degradation", "Kd","hill", "protein_production", "protein_degradation", "Kd", "hill"]), model_mode=models[model_index])                                       
	solver = Solver(model)    

	region_files = []
	region_files.append(os.path.join(folder, "bioprocViableSet_IterGA.p"))   
	
	for i in range(10):
		region_files.append(os.path.join(folder, "bioproc_Region0ViableSet_Iter" + str(i+1) + ".p"))

	viablePoints = [] 	
	for region_file in region_files: 
		viablePointsRegion = pickle.load(open(region_file, "rb"))   
		print(len(viablePointsRegion))   
		viablePoints.extend(viablePointsRegion)
 	                      
	region = Region(viablePoints, model, "region")          	
	model_regions.append(region)        
		
def calculateVolumes(model_index=0): 	
	model = BioProc(np.array(["protein_production", "protein_production", "protein_production", "protein_production", "protein_degradation", "protein_degradation", "Kd","hill", "protein_production", "protein_degradation", "Kd", "hill"]), model_mode=models[model_index], avg_dev=30)    
		
	print(model.threshold)  
	solver = Solver(model)
	folder = folders[model_index]
	print(folder)    
	
	_, vol_desc = solver.getViableVolume([model_regions[model_index]])           
	
	f = open(os.path.join(folder, "viable_volume.txt"), "w")      
	f.write(vol_desc)   
	f.close()      

def plotBoxPlots():	 
	number_points = int(1e4)     
	
	all_costs = [] 
	for box_plot in range(len(models)):   	
		plt.subplot(1,3,box_plot+1)
		
		model_evals = []     
		for model_index in range(len(models)):     			
			evals = []
			region = model_regions[model_index]   
			print(type(region.points))   
			samples = region.points[np.random.choice(region.points.shape[0], number_points, replace=False), :]            
			
			model = BioProc(np.array(["protein_production", "protein_production", "protein_production", "protein_production", "protein_degradation", "protein_degradation", "Kd","hill", "protein_production", "protein_degradation", "Kd", "hill"]), model_mode=models[box_plot])                                      	
			for sample in samples:
				evals.append(-model.eval(sample)[0]) 
				
			model_evals.append(evals)      
		
		all_costs.append(model_evals)      
		print(model_evals[0])     	
		plt.boxplot(model_evals)            	

	#dump data to file  		
	f = open(os.path.join(base_path, "boxplot_data.p"), "wb")        
	pickle.dump(all_costs, f)   
	f.close()
	
	plt.show()  	
	
def plotStochasticSimulations(from_file = True): 
	number_points = 3 
	print("Plotting stochastic simulations") 
	
	fig, axs = plt.subplots(len(models), number_points)		

	plotFlipflops = False  
	
	for model_index in range(len(models)):  
		print(model_index) 
		region = model_regions[model_index] 
		model = BioProc(np.array(["protein_production", "protein_production", "protein_production", "protein_production", "protein_degradation", "protein_degradation", "Kd","hill", "protein_production", "protein_degradation", "Kd", "hill"]), model_mode=models[model_index], plot_fitness=False, plot_devs=False)   
		
		samples = []
		if not from_file:
			samples = region.points[np.random.choice(region.points.shape[0], number_points, replace=False), :]    
		else:
			for i in range(number_points):
				samples.append(pickle.load(open("model" + str(model_index + 1) + "sample" + str(i + 1) + ".p", "rb")))  	
			
		
		sample_num = 0
		for sample in samples:
			print(sample) 
			T, Y = model.simulateStochastic(sample)  
			Y_ode = model.simulate(sample)   

			colors = plt.rcParams["axes.prop_cycle"]() #get color iterator    					
			for i in range(model_index + 1):  
				
				if plotFlipflops: 				
					c = next(colors)["color"]
					axs[model_index, sample_num].plot(T, Y[:, i*4 + 2], label='q'+str(i+1), color=c, alpha=0.7) 
					axs[model_index, sample_num].plot(model.ts, Y_ode[:, i*4 + 2], color=c, linestyle="--")  
					#axs[model_index, sample_num].get_xaxis().set_visible(False)
			
			
			for i in range(2*(model_index + 1)):
				c = next(colors)["color"]
				axs[model_index,sample_num].plot(T, Y[:, -(2*(model_index + 1) - i)], label="i"+str(i + 1), color=c, alpha=0.7)  
				axs[model_index,sample_num].plot(model.ts, Y_ode[:, -(2*(model_index + 1) - i)], color=c, linestyle="--")  
			
			if not from_file:
				pickle.dump(sample, open("model" + str(model_index + 1) + "sample" + str(sample_num + 1) + ".p", "wb+")) 
			
			sample_num += 1 
			
	for i in range(len(samples)):
		axs[len(models) - 1, i].set(xlabel='Time [h]') 
	for i in range(len(models)):
		axs[i, 0].set(ylabel='Concentration [nM]')   
		#plt.legend(loc="upper left")
		#plt.legend()  			
	plt.show()   

if __name__ == '__main__':
			
	tic = time.perf_counter()			
	plotStochasticSimulations(from_file = True)       
	toc = time.perf_counter()
	print("Elapsed time: " + str(toc - tic) + "s")   


   


