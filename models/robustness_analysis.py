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

import multiprocessing


 
if __name__ == '__main__':  
    sns.set_style("white")
    flatui = ['#d9d9d9','#bdbdbd','#969696','#636363']
    #sns.palplot(sns.color_palette(flatui))
    sns.set_palette(flatui)

    #
    # SETTINGS
    #
    ga_solutions = False
    local_solutions = True

    
    base_paths_opt = [os.path.join(".", "results_opt")]
    #base_paths_opt = [os.path.join(".", "results_opt"), os.path.join(".", "results_opt_rep1"), os.path.join(".", "results_opt_rep2")]
    num_models_fitness = 4
    num_models_regions = 4
    
    num_models_fitness = 5
    num_models_regions = 5
    


    base_path_robustness = os.path.join(".", "results_robustness") 


    #
    # END OF SETTINGS
    #


    models = [one_bit_processor_ext, two_bit_processor_ext, three_bit_processor_ext, four_bit_processor_ext, five_bit_processor_ext]   

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
        print(region_files)  
        for region_file in region_files: 
            viablePointsRegion = pickle.load(open(region_file, "rb"))   
            print(len(viablePointsRegion))   
            viablePoints.extend(viablePointsRegion)
        print("Number of points ("+str(model_index+1)+"-bit):",len(viablePoints))
        region = Region(viablePoints, model, "region")              
        model_regions.append(region)                                                                        
        

def calculateVolumes(model_indexes=None):    
    if model_indexes == None:
        model_indexes = range(num_models_regions)
    elif type(model_indexes) == int:
        model_indexes = [model_indexes]

    for model_index in model_indexes:
        model = BioProc(np.array(["protein_production", "protein_production", "protein_production", "protein_production", "protein_degradation", "protein_degradation", "Kd","hill", "protein_production", "protein_degradation", "Kd", "hill"]), model_mode=models[model_index], avg_dev=30)    
        
        print(model.threshold)  
        solver = Solver(model)
        _, vol_desc = solver.getViableVolume([model_regions[model_index]])     
        
        model_str = '0'+str(model_index+1)+'_'
    
        f = open(os.path.join(base_path_robustness, model_str+"viable_volume.txt"), "w")      
        f.write(vol_desc)   
        f.close()    

def calculateVolumesRepeats(model_indexes=None, repeats = 10):    
    if model_indexes == None:
        model_indexes = range(num_models_regions)
    elif type(model_indexes) == int:
        model_indexes = [model_indexes]

    df = pd.DataFrame(columns=["Model id", "Volume", "Total", "Ratio"])

    for model_index in model_indexes:
        model = BioProc(np.array(["protein_production", "protein_production", "protein_production", "protein_production", "protein_degradation", "protein_degradation", "Kd","hill", "protein_production", "protein_degradation", "Kd", "hill"]), model_mode=models[model_index], avg_dev=30)    
        
        print(model.threshold)  
        solver = Solver(model)

        for _ in range(repeats):
            (volume, total, ratio), _ = solver.getViableVolume([model_regions[model_index]])        
            df = df.append({"Model id":model_index+1, "Volume": volume, "Total":total, "Ratio":ratio}, ignore_index=True)
    
    df.to_csv(os.path.join(base_path_robustness, 'volumes.csv'), index=False)

def plotVolumesFromCsv(file_name = ""):
    if not file_name:
        file_name = os.path.join(base_path_robustness, 'volumes.csv')
    df = pd.read_csv(file_name)
    model_indexes = list(df["Model id"].unique())
    model_indexes.sort()
    df_avg = pd.DataFrame(columns=["Model id", "Volume", "Total", "Ratio"])   

    for model_index in model_indexes:
        df_model = df[df["Model id"]==model_index]
        n = df_model.shape[0]
        df_avg = df_avg.append({"Model id":model_index, 
                                "Volume": sum(df_model.Volume)/n, 
                                "Total":sum(df_model.Total)/n, 
                                "Ratio":sum(df_model.Ratio)/n}, ignore_index=True)
    df_avg.to_csv(os.path.join(base_path_robustness, 'volumes_avg.csv'), index=False)    

    sns.barplot(x = 'Model id', y = 'Ratio', data = df_avg)#, palette="Pastel1")
    plt.ylabel('Volume [a.u.]')
    plt.yscale('log')
    #fig = plt.gcf()
    #fig.set_size_inches([20,8])
    plt.savefig(os.path.join(base_path_robustness, 'volumes_avg.pdf'), bbox_inches = 'tight')
    plt.close()

    #sns.violinplot(x="Model id", y="Ratio", data=df_avg, palette="Pastel1")
    #plt.ylabel('Volume [a.u.]')
    #plt.yscale('log')
    #plt.show()



def plotVolumesFromTxt(model_indexes=None):  
    if model_indexes == None:
        model_indexes = range(num_models_regions)
    elif type(model_indexes) == int:
        model_indexes = [model_indexes]

    df = pd.DataFrame(columns = ["Model id", "Total", "Ratio"])

    for model_index in model_indexes:
        model_str = '0'+str(model_index+1)+'_'
        f = open(os.path.join(base_path_robustness, model_str+"viable_volume.txt"))
        f.readline()
        volume = float(f.readline().strip().split(":")[-1])
        total = float(f.readline().strip().split(":")[-1])
        ratio = float(f.readline().strip().split(":")[-1])
        df = df.append({"Model id":model_index+1, "Volume": volume, "Total":total, "Ratio":ratio}, ignore_index=True)

    sns.barplot(x = 'Model id', y = 'Ratio', data = df)#, palette="Pastel1")
    plt.ylabel('Volume [a.u.]')
    plt.yscale('log')
    #fig = plt.gcf()
    #fig.set_size_inches([20,8])
    plt.savefig(os.path.join(base_path_robustness, 'volumes.pdf'), bbox_inches = 'tight')
    plt.show()

    
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

def getCostsParallel(number_points = 0, file_name = ""):     
    rand_samples = []
    for model_index in range(num_models_regions):
        region = model_regions[model_index]   
        if number_points:
            samples = region.points[np.random.choice(region.points.shape[0], number_points, replace=False), :]            
        else:
            samples = region.points

        rand_samples.append(samples)
    
    df = pd.DataFrame(columns=['Model id', 'Region id', 'cost'])
    
    pool = multiprocessing.Pool()
            

    for model_id in range(num_models_fitness):
        model = BioProc(np.array(["protein_production", "protein_production", "protein_production", "protein_production", "protein_degradation", "protein_degradation", "Kd","hill", "protein_production", "protein_degradation", "Kd", "hill"]), model_mode=models[model_id])                                          
        for region_id in range(num_models_regions):
            costs = pool.map(model.eval, rand_samples[region_id])
            
            for costA in costs:
                cost = -costA[0]
                df = df.append({'Model id': model_id+1, 'Region id': region_id+1, 'cost':cost}, ignore_index=True)

    pool.close()

    if file_name:
        df.to_csv(file_name, index=False)

    return df



def plotCostdf(df=None, number_points = 0, normalize=True):
    if not type(df):
        df = getCosts(number_points)
    
    df['Model id'] = df['Model id'].astype(int)
    df['Region id'] = df['Region id'].astype(int)

    
    thresholds = [30, 20, 18, 17]
    #thresholds = [30, 20, 20, 20]

   

    _, axes = plt.subplots(1,2, gridspec_kw={'width_ratios': [2, 1]})
    
    df_comp = pd.DataFrame(columns =['Model id', 'Region id', 'compatible'])


    for model_id in range(1,num_models_fitness+1):        
        for region_id in range(1,num_models_regions+1):
            threshold = max(thresholds[region_id-1], thresholds[model_id-1])
            n = df[(df['Model id']== model_id) & (df['Region id']== region_id)].shape[0]
            comp = df[(df['Model id']== model_id) & (df['Region id']== region_id) & (df['cost'] <= threshold)].shape[0]
            df_comp = df_comp.append({'Model id': model_id, 'Region id': region_id, 'compatible': comp/n}, ignore_index=True)

    


    g = sns.barplot(x = 'Model id', y = 'compatible', data = df_comp, hue="Region id", ax = axes[1]) #, palette="Pastel1")
    g.legend_.remove()
    axes[1].set_ylabel('Fractions')
    #l = axes[1].get_legend()
    #l.set_bbox_to_anchor((1, 0.75))

    if normalize:
        for model_id in range(1,num_models_fitness+1):
            locs = df['Model id']== model_id
            df.loc[locs,'cost'] /= thresholds[model_id-1]

    g=sns.violinplot(x="Model id", y="cost", hue="Region id", data=df, ax = axes[0]) #, palette="Pastel1")
    #g.legend_.remove()
    g.legend(ncol=10, loc='upper center', bbox_to_anchor=(0.5, 0.95), title="Region id")
    if normalize:
        axes[0].set_ylabel('Normalized costs [a.u.]')    
    else:
        axes[0].set_ylabel('Costs')


    """
    plt.legend(ncol=10, 
          loc='upper center',
          bbox_to_anchor=(0.5, 0.95),
          bbox_transform=plt.gcf().transFigure)
    """    
    fig = plt.gcf()
    fig.set_size_inches([20,8])
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

    if num_models_regions >= 4:

        df4 = pd.DataFrame(rand_samples[3])
        df4.columns = param_names
        df4["Model id"] = 4

        df = pd.concat([df1, df2, df3, df4], ignore_index=True)

    
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

        sns.violinplot(y = param_names[param_id], x="Model id", data=df[[param_names[param_id], "Model id"]], ax = ax) #,palette="Pastel1")
    
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
    

"""   
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
"""           
def plotStochasticSimulations(from_file = True, number_points = 3, plotFlipflops = False, pickle_dump=False): 

    sns.set_style("white")
    #flatui = ['#f7f7f7','#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525']
    #flatui = ['#bdbdbd','#969696','#737373','#525252','#252525','#000000']
    flatui = ['#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525']
    #sns.palplot(sns.color_palette(flatui))
    sns.set_palette(flatui)
    
     
    print("Plotting stochastic simulations") 
    
    fig, axs = plt.subplots(3, number_points, sharey='row')     

    for model_index in range(3):  
        print(model_index) 
        region = model_regions[model_index] 
        model = BioProc(np.array(["protein_production", "protein_production", "protein_production", "protein_production", "protein_degradation", "protein_degradation", "Kd","hill", "protein_production", "protein_degradation", "Kd", "hill"]), model_mode=models[model_index], plot_fitness=False, plot_devs=False)   
        
        samples = []
        if not from_file:
            samples = region.points[np.random.choice(region.points.shape[0], number_points, replace=False), :]    
        else:
            for i in range(number_points):
                samples.append(pickle.load(open("selected_points/model" + str(model_index + 1) + "sample" + str(i + 1) + ".p", "rb")))      
            
        
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
                pickle.dump(sample, open("selected_points/test_model" + str(model_index + 1) + "sample" + str(sample_num + 1) + ".p", "wb+")) 
            
            sample_num += 1 
            
    for i in range(len(samples)):
        axs[2, i].set(xlabel='Time [h]') 
    for i in range(3):
        axs[i, 0].set(ylabel='Concentration [nM]')   
        #plt.legend(loc="upper left")
        #plt.legend()           
    

    plt.legend(ncol=10, 
          loc='upper center',
          bbox_to_anchor=(0.5, 0.95),
          bbox_transform=plt.gcf().transFigure)

    plt.gcf().set_size_inches(15,12) 
    plt.savefig(os.path.join(base_path_robustness, 'ssa.pdf'), bbox_inches = 'tight')
    if pickle_dump:
        pickle.dump(fig, open("selected_points/plot_SSA.pickle", "wb"))

    plt.show()  

    

if __name__ == "__main__":
    print('------------------------------------------------------------- ok 0')
    df = getCostsParallel(file_name="results_robustness/costs.csv")  
    print('------------------------------------------------------------- ok 1')
    df = pd.read_csv("results_robustness/costs.csv")
    print('------------------------------------------------------------- ok 2')
    plotCostdf(df)
    print('------------------------------------------------------------- ok 3')


    df = getParamDistrib(file_name="results_robustness/params.csv")
    print('------------------------------------------------------------- ok 4')
    df = pd.read_csv("results_robustness/params.csv")
    print('------------------------------------------------------------- ok 5')

    plotParamsdf(df)
    print('------------------------------------------------------------- ok 6')
    plotStochasticSimulations(pickle_dump=False)   
