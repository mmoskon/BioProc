from bioproc.proc_opt import BioProc  
from bioproc.proc_models import *  
from bioproc.hill_functions import * 
import matplotlib.pyplot as plt 
from solver import * 
import numpy as np  
import os.path      
import pickle     
import seaborn as sns
import pandas as pd

sns.set_style("white")
 
        
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

def getCosts(number_points = 0, file_name = ""):     
    rand_samples = []
    for model_index in range(num_models_regions):
        region = model_regions[model_index]   
        if number_points:
            samples = region.points[np.random.choice(region.points.shape[0], number_points, replace=False), :]            
        else:
            samples = region.points

        rand_samples.append(samples)
    
    df = pd.DataFrame(columns=['Model id', 'Region id', 'cost'])

    for model_id in range(num_models_fitness):
        model = BioProc(np.array(["protein_production", "protein_production", "protein_production", "protein_production", "protein_degradation", "protein_degradation", "Kd","hill", "protein_production", "protein_degradation", "Kd", "hill"]), model_mode=models[model_id])                                          
        for region_id in range(num_models_regions):
            for sample in rand_samples[region_id]:
                cost = -model.eval(sample)[0]
                df = df.append({'Model id': model_id+1, 'Region id': region_id+1, 'cost':cost}, ignore_index=True)

    if file_name:
        df.to_csv(file_name, index=False)

    return df

def plotCostdf(df=None, number_points = 0):
    if not type(df):
        df = getCosts(number_points)
    
    df['Model id'] = df['Model id'].astype(int)
    df['Region id'] = df['Region id'].astype(int)

    sns.violinplot(x="Model id", y="cost", hue="Region id", data=df, palette="Pastel1")
    plt.savefig(os.path.join(base_path_robustness, 'cost_distrib_sns.pdf'), bbox_inches = 'tight')
    plt.show()





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
        ax.boxplot(model_evals) 
        ax.set_xticks([1,2,3])
                      
        all_costs.append(model_evals)      
        
        

    #dump data to file          
    f = open(os.path.join(base_path_robustness, "cost_distrib.p"), "wb")        
    pickle.dump(all_costs, f)   
    f.close()

    plt.savefig(os.path.join(base_path_robustness, 'cost_distrib.pdf'), bbox_inches = 'tight')
    plt.show()      




def getParamDistrib(number_points = 0, file_name = ""):    
    rand_samples = []    
    for model_index in range(num_models_regions):
        region = model_regions[model_index]   
        if number_points:
            samples = region.points[np.random.choice(region.points.shape[0], number_points, replace=False), :]
        else:
            samples = region.points

        rand_samples.append(samples)

    rand_samples = np.array(rand_samples)

    param_names = [r"$\alpha_1$", r"$\alpha_2$", r"$\alpha_3$", r"$\alpha_4$", r"$\delta_1$", r"$\delta_2$", r"$K_d$", r"$n$", r"$\alpha_I$", r"$\delta_I$", r"$K_{d_I}$", r"$n_I$"]
    
    df1 = pd.DataFrame(rand_samples[0])
    df1.columns = param_names
    df1["Model id"] = 1
    
    df = df1

    if num_models_regions >= 2:

        df2 = pd.DataFrame(rand_samples[1])
        df2.columns = param_names
        df2["Model id"] = 2

        df = pd.concat([df1, df2], ignore_index=True)

    if num_models_regions >= 3:

        df3 = pd.DataFrame(rand_samples[2])
        df3.columns = param_names
        df3["Model id"] = 3

        df = pd.concat([df1, df2, df3], ignore_index=True)

    
    
    """
    for model_id in range(num_models_regions):
        df_tmp = pd.DataFrame(columns = ['Parameter', 'Value'])
        for sample in rand_samples[model_id]:
            for name, value in zip(param_names, sample):
                df_tmp.append({'Parameter': name, 'Value': value}, ignore_index=True)
        
        df_tmp['Model id'] = model_id+1
        if model_id == 0:
            df = df_tmp
        else:
            df = pd.concat([df, df_tmp], ignore_index=True)
    """
    if file_name:
        df.to_csv(file_name, index=False)
    return df

def plotParamsdf(df=None, number_points = 0):
    if not type(df):
        df = getParamDistrib(number_points)
    
    df['Model id'] = df['Model id'].astype(int)


    param_names = [r"$\alpha_1$", r"$\alpha_2$", r"$\alpha_3$", r"$\alpha_4$", r"$\delta_1$", r"$\delta_2$", r"$K_d$", r"$n$", r"$\alpha_I$", r"$\delta_I$", r"$K_{d_I}$", r"$n_I$"]
    
    fig, axes = plt.subplots(3,4)

    for param_id in range(len(param_names)):
        ax = axes.flat[param_id]

        sns.violinplot(y = param_names[param_id], x="Model id", data=df[[param_names[param_id], "Model id"]], palette="Pastel1", ax = ax)
    
    fig.set_size_inches([20,12])
    plt.savefig(os.path.join(base_path_robustness, 'params_distrib_sns.pdf'), bbox_inches = 'tight')   
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
    
    fig, axes = plt.subplots(3,4)
    

    for param_id in range(len(param_names)):
        ax = axes.flat[param_id]
        #ax.set_yscale('log')
        #plt.subplot(1,len(param_names),param_id+1)
        ax.set_ylabel(param_names[param_id])

        if param_id > 7:
            ax.set_xlabel('model id')

        
        #plt.boxplot([rand_samples[0,:,i],rand_samples[1,:,i], rand_samples[2,:,i]])
        ax.violinplot([rand_samples[0][:,param_id]])#,rand_samples[1][:,param_id]])
        #ax.boxplot([rand_samples[0][:,param_id]])#,rand_samples[1][:,param_id]])
        ax.set_xticks([1,2,3])
        #ax.set_xticklabels([1,2,3])



    fig.set_size_inches([20,12])
    plt.savefig(os.path.join(base_path_robustness, 'params_distrib.pdf'), bbox_inches = 'tight')
    plt.show()      
    

    
def plotStochasticSimulations(number_points = 5):
 	
	print("Plotting stochastic simulations")
	for model_index in range(len(models)): 
		print(model_index) 
		region = model_regions[model_index] 
		model = BioProc(np.array(["protein_production", "protein_production", "protein_production", "protein_production", "protein_degradation", "protein_degradation", "Kd","hill", "protein_production", "protein_degradation", "Kd", "hill"]), model_mode=models[model_index], plot_fitness=True)  
		samples = region.points[np.random.choice(region.points.shape[0], number_points, replace=False), :]	 
		
		for sample in samples:
			T, Y = model.simulateStochastic(sample) 
			Y_ode = model.simulate(sample, hot_start=False)	 
			
			for i in range(model_index + 1):  
				p = plt.plot(T, Y[:, i*4 + 2], label='q'+str(i+1), alpha=0.7) 
				plt.plot(model.ts, Y_ode[:, i*4 + 2], p[0].get_color(), linestyle="--") 

			for i in range(2*(model_index + 1)):
				p = plt.plot(T, Y[:, -(i+1)], label="i"+str(2*model_index + 2 - i), alpha=0.7)	
				plt.plot(model.ts, Y_ode[:, -(i+1)],  p[0].get_color(), linestyle="--")	 
			
			#plt.legend(loc="upper left")
			plt.legend()			
			plt.show()	 
            
if __name__ == '__main__':  


    #
    # SETTINGS
    #
    ga_solutions = False
    local_solutions = True

    #base_paths_opt = [os.path.join(".", "results_opt"), os.path.join(".", "results_opt2")]#, os.path.join(".", "results_opt3")]
    base_paths_opt = [os.path.join(".", "results_opt")]
    num_models_fitness = 3
    num_models_regions = 3


    base_path_robustness = os.path.join(".", "results_robustness") 


    #
    # END OF SETTINGS
    #


    models = [one_bit_processor_ext, two_bit_processor_ext, three_bit_processor_ext]   

    #folders = [os.path.join(base_path, "one_bit_model"), os.path.join(base_path, "two_bit_model"), os.path.join(base_path, "three_bit_model")]   
    model_regions = []

        
    for model_index in range(num_models_regions):                    
        #folder = folders[model_index]               
        model = BioProc(np.array(["protein_production", "protein_production", "protein_production", "protein_production", "protein_degradation", "protein_degradation", "Kd","hill", "protein_production", "protein_degradation", "Kd", "hill"]), model_mode=models[model_index])                                       
        solver = Solver(model)    

        model_str = '0'+str(model_index+1)+'_'
        region_files =  []
        for base_path_opt in base_paths_opt:
            if ga_solutions:
                region_files.append(os.path.join(base_path_opt, model_str+"bioprocViableSet_IterGA.p"))
            if local_solutions:
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







    #calculateVolumes(model_index=0)  
    #calculateVolumes(model_index=1)    
    #calculateVolumes(model_index=2)      

    #plotCost(number_points = 5)  

    #plotParamDistrib()

    #plotStochasticSimulations()   


    #df = getCosts(number_points = 5, file_name = "costs.csv")


    #df = getCosts(number_points=5, file_name="results_robustness\\costs.csv")
    #df = pd.read_csv("results_robustness\\costs.csv")
    #plotCostdf(df)

    df = getParamDistrib(file_name="results_robustness\\params.csv")
    df = pd.read_csv("results_robustness\\params.csv")
    plotParamsdf(df)