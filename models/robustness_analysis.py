from bioproc.proc_opt import BioProc  
from bioproc.proc_models import *  
from bioproc.hill_functions import * 
import matplotlib.pyplot as plt 
from solver import * 
import numpy as np  
import os.path      
import pickle     
import seaborn as sns

sns.set_style("white")


models = [one_bit_processor_ext, two_bit_processor_ext, three_bit_processor_ext]   
base_path_opt = os.path.join(".", "results_opt") 
base_path_robustness = os.path.join(".", "results_robustness") 
#folders = [os.path.join(base_path, "one_bit_model"), os.path.join(base_path, "two_bit_model"), os.path.join(base_path, "three_bit_model")]   
model_regions = []

num_models_fitness = 3
num_models_regions = 1

    
for model_index in range(num_models_regions):                    
    #folder = folders[model_index]               
    model = BioProc(np.array(["protein_production", "protein_production", "protein_production", "protein_production", "protein_degradation", "protein_degradation", "Kd","hill", "protein_production", "protein_degradation", "Kd", "hill"]), model_mode=models[model_index])                                       
    solver = Solver(model)    

    model_str = '0'+str(model_index+1)+'_'

    region_files = [os.path.join(base_path_opt, model_str+"bioprocViableSet_IterGA.p")] 

    for i in range(10):
        region_files.append(os.path.join(base_path_opt, model_str+"bioproc_Region0ViableSet_Iter" + str(i+1) + ".p"))

    viablePoints = []   
    for region_file in region_files: 
        viablePointsRegion = pickle.load(open(region_file, "rb"))   
        print(len(viablePointsRegion))   
        viablePoints.extend(viablePointsRegion)
    print("Number of points ("+str(model_index+1)+"-bit):",len(viablePoints))
    region = Region(viablePoints, model, "region")              
    model_regions.append(region)        
        
def calculateVolumes(model_index=0):    
    model = BioProc(np.array(["protein_production", "protein_production", "protein_production", "protein_production", "protein_degradation", "protein_degradation", "Kd","hill", "protein_production", "protein_degradation", "Kd", "hill"]), model_mode=models[model_index], avg_dev=30)    
        
    print(model.threshold)  
    solver = Solver(model)
    #folder = folders[model_index]
    #print(folder)    
    
    _, vol_desc = solver.getViableVolume([model_regions[model_index]])           
    
    model_str = '0'+str(model_index+1)+'_'
    
    f = open(os.path.join(base_path_robustness, model_str+"viable_volume.txt"), "w")      
    f.write(vol_desc)   
    f.close()      
"""
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
    f = open(os.path.join(base_path_robustness, "boxplot_data.p"), "wb")        
    pickle.dump(all_costs, f)   
    f.close()
    
    plt.show()      
""" 
def plotCost(number_points = 0):     
    
    rand_samples = []
    for model_index in range(num_models_regions):
        region = model_regions[model_index]   
        if number_points:
            samples = region.points[np.random.choice(region.points.shape[0], number_points, replace=False), :]            
        else:
            samples = region.points

        rand_samples.append(samples)
    
    all_costs = [] 

    _, axes = plt.subplots(1,3)
      
    for box_plot in range(num_models_fitness):   
        ax = axes.flat[box_plot]

        #plt.subplot(1,3,box_plot+1)
        ax.set_title(str(box_plot+1)+'-bit model')
        ax.set_xlabel('region id')
        if box_plot == 0:
            ax.set_ylabel('cost')
    
        model = BioProc(np.array(["protein_production", "protein_production", "protein_production", "protein_production", "protein_degradation", "protein_degradation", "Kd","hill", "protein_production", "protein_degradation", "Kd", "hill"]), model_mode=models[box_plot])                                          
        
        model_evals = []    
        for model_index in range(num_models_regions):
            evals = []
            for sample in rand_samples[model_index]:
                evals.append(-model.eval(sample)[0]) 
            
            model_evals.append(evals)      
        
        #plt.boxplot(model_evals)
        ax.violinplot(model_evals)                
        all_costs.append(model_evals)      
        
        

    #dump data to file          
    f = open(os.path.join(base_path_robustness, "cost_distrib.p"), "wb")        
    pickle.dump(all_costs, f)   
    f.close()

    plt.savefig(os.path.join(base_path_robustness, 'cost_distrib.pdf'), bbox_inches = 'tight')
    plt.show()      

def plotParamDistrib(number_points = 0):     
    
    rand_samples = []    
    for model_index in range(num_models_regions):
        region = model_regions[model_index]   
        if number_points:
            samples = region.points[np.random.choice(region.points.shape[0], number_points, replace=False), :]
        else:
            samples = region.points

        rand_samples.append(samples)

    rand_samples = np.array(rand_samples)
        
    #["protein_production", "protein_production", "protein_production", "protein_production", "protein_degradation", "protein_degradation", "Kd","hill", "protein_production", "protein_degradation", "Kd", "hill"]
    param_names = [r"$\alpha_1$", r"$\alpha_2$", r"$\alpha_3$", r"$\alpha_4$", r"$\delta_1$", r"$\delta_2$", r"$K_d$", r"$n$", r"$\alpha_I$", r"$\delta_I$", r"$K_{d_I}$", r"$n_I$"]
    
    _, axes = plt.subplots(3,4)

    for param_id in range(len(param_names)):
        ax = axes.flat[param_id]
        
        #plt.subplot(1,len(param_names),param_id+1)
        ax.set_ylabel(param_names[param_id])
        ax.set_xlabel('model id')

        #plt.boxplot(rand_samples[0,:,i],rand_samples[1,:,i], rand_samples[2,:,i])
        ax.violinplot(rand_samples[0,:,param_id])
  
    plt.savefig(os.path.join(base_path_robustness, 'params_distrib.pdf'), bbox_inches = 'tight')
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
            #   plt.plot(T, Y[:, i*4 + 2], label='q'+str(i+1))     

            for i in range(2*(model_index + 1)):
                plt.plot(T, Y[:, -(i+1)], label="i"+str(2*model_index + 2 - i))  
            
            #plt.legend(loc="upper left")
            plt.legend()            
            plt.show()  
            
            
#calculateVolumes(model_index=0)  
#calculateVolumes(model_index=1)    
#calculateVolumes(model_index=2)      

plotCost(number_points = 100)  

#plotParamDistrib()#number_points = 10)

#plotStochasticSimulations()   


