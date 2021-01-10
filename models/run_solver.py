import solver
from bioproc.proc_opt import BioProc 
from bioproc.proc_models import * 
import os
import numpy as np

if __name__ == "__main__":

    param_values = {  "transcription": {"min": 0.01, "max": 50},   
                    "translation": {"min": 0.01, "max": 50},  
                    "protein_production": {"min": 0.1, "max": 50},                          
                    "rna_degradation": {"min": 0.1, "max": 100},        
                    "protein_degradation": {"min": 0.001, "max": 50},         
                    "hill": {"min": 1, "max": 5},         
                    "Kd": {"min": 0.01, "max": 250}, 
                    "protease_concentration": {"min": 10, "max":1000}       
                    }         

    #model_mode = one_bit_processor_ext
    #model_mode = two_bit_processor_ext
    #model_mode = three_bit_processor_ext
    model_mode = four_bit_processor_ext
    
    folder_name = "results_opt"
    if model_mode == one_bit_processor_ext:
        file_name = "01_bioproc"
    elif model_mode == two_bit_processor_ext:
        file_name = "02_bioproc"
    elif model_mode == three_bit_processor_ext:
        file_name = "03_bioproc"
    elif model_mode == three_bit_processor_ext:
        file_name = "04_bioproc"
    else:
        file_name = "05_bioproc"

    filename =  os.path.join(".", folder_name, file_name)                                             
    print(filename)        

    model = BioProc(np.array(["protein_production", "protein_production", "protein_production", "protein_production", "protein_degradation", "protein_degradation", "Kd","hill", "protein_production", "protein_degradation", "Kd", "hill"]), model_mode=model_mode, parameter_values=param_values, avg_dev=30)                                      
    sol = solver.Solver(model)                        
    sol.run(filename, maxDepth=1) #do not cluster    