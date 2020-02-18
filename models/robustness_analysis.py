from bioproc.proc_opt import BioProc  
from bioproc.proc_models import *  
from bioproc.hill_functions import * 
import matplotlib.pyplot as plt 
from solver import * 
import numpy as np  
import os.path      
import pickle     

	
models = [one_bit_processor_ext, two_bit_processor_ext, three_bit_processor_ext]   
base_path = os.path.join(".", "bioproc") 
folders = [os.path.join(base_path, "one_bit_model"), os.path.join(base_path, "two_bit_model"), os.path.join(base_path, "three_bit_model")]   
model_regions = []      
	
for model_index in range(3):     				
	folder = folders[model_index]        		
	model = BioProc(np.array(["protein_production", "protein_production", "protein_production", "protein_production", "protein_degradation", "protein_degradation", "Kd","hill", "protein_production", "protein_degradation", "Kd", "hill"]), model_mode=models[model_index])                                       
	solver = Solver(model)    

	region_files = [os.path.join(folder, "bioprocViableSet_IterGA.p")] 
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
	for box_plot in range(3):  	
		plt.subplot(1,3,box_plot+1)
		
		model_evals = []    
		for model_index in range(3):    			
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
	
def plotStochasticSimulations():
	number_points = 5
	print("Plotting stochastic simulations") 
	for model_index in range(3):
		print(model_index) 
		region = model_regions[model_index] 
		model = BioProc(np.array(["protein_production", "protein_production", "protein_production", "protein_production", "protein_degradation", "protein_degradation", "Kd","hill", "protein_production", "protein_degradation", "Kd", "hill"]), model_mode=models[model_index], plot_fitness=True)  
		samples = region.points[np.random.choice(region.points.shape[0], number_points, replace=False), :]   
		
		for sample in samples:
			T, Y = model.simulateStochastic(sample) 
			
			#for i in range(model_index + 1):  
			#	plt.plot(T, Y[:, i*4 + 2], label='q'+str(i+1))     

			for i in range(2*(model_index + 1)):
				plt.plot(T, Y[:, -(i+1)], label="i"+str(2*model_index + 2 - i))  
			
			#plt.legend(loc="upper left")
			plt.legend()  			
			plt.show()  
			
			
#calculateVolumes(model_index=0)  
#calculateVolumes(model_index=1)   	
#calculateVolumes(model_index=2)      

#plotBoxPlots()  

plotStochasticSimulations()   


