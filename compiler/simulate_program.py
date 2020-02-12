from bioproc.hill_functions import *
from generate_model import generate_model
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np
import importlib


def simulate_program(program_name, t_end, N, params_ff, params_addr, params_proteolysis, params_condition, params_prog, n_bits, ax=plt, plot_clock = True):
    
    model_name = program_name.split(".")[0]
    
    prog_alpha, prog_delta, prog_n, prog_Kd = params_prog 


    generate_model(program_name, model_name, n_bits, prog_alpha, prog_delta, prog_n, prog_Kd)      
    model = importlib.import_module(model_name.replace("\\","."))

    f_description = open(model_name+"description.txt")
    ops = f_description.readline().strip().split(",")[:-1]
    f_description.close()

    # params
    params = (params_ff + params_proteolysis + params_addr + params_condition,)
    
    # solving
    if n_bits == 3:
        Y0 = np.array([0]*(18+len(ops)))
    else: #if n_bits == 4:
        Y0 = np.array([0]*(24+len(ops)))

    T = np.linspace(0, t_end, N)
    Y = odeint(model.model, Y0, T, args=params)

    i = -len(ops)
    for op in ops:
        o = Y[:,i]
        ax.plot(T,o)
        i += 1
        


    if plot_clock:
        clk = get_clock(T)
        ax.plot(T,clk, color="black", alpha=0.1)
        ax.legend(ops + ['clk'], loc='upper left')

    else:

        if n_bits == 3:
            i1 = Y[:,-6-len(ops)]
            i2 = Y[:,-5-len(ops)]
            i3 = Y[:,-4-len(ops)]
            i4 = Y[:,-3-len(ops)]
            i5 = Y[:,-2-len(ops)]
            i6 = Y[:,-1-len(ops)]
        else: #if n_bits == 4:
            i1 = Y[:,-8-len(ops)]
            i2 = Y[:,-7-len(ops)]
            i3 = Y[:,-6-len(ops)]
            i4 = Y[:,-5-len(ops)]
            i5 = Y[:,-4-len(ops)]
            i6 = Y[:,-3-len(ops)]
            i7 = Y[:,-2-len(ops)]
            i8 = Y[:,-1-len(ops)]


        ax.plot(T,i1)
        ax.plot(T,i2)
        ax.plot(T,i3)
        ax.plot(T,i4)
        ax.plot(T,i5)
        ax.plot(T,i6)
        if n_bits == 4:
            ax.plot(T,i7)
            ax.plot(T,i8)

       


        
        if n_bits == 3:
            ax.legend(['i1','i2','i3','i4','i5','i6']+ops, loc='upper left')
        else:
            ax.legend(['i1','i2','i3','i4','i5','i6','i7','i8']+ops, loc='upper left')
    
    if ax == plt:
        plt.savefig("figs\\"+program_name.split(".")[0]+".pdf", bbox_inches = 'tight')
        plt.savefig("figs\\"+program_name.split(".")[0]+".png")
        
        plt.show()




if __name__ == '__main__':

    # simulation parameters

    t_end = 200
    N = 1000

    n_bits = 4
    #program_name = "programs\\program_add.txt"
    #program_name = "programs\\program_add_nop.txt"
    #program_name = "programs\\program_sub.txt"
    #program_name = "programs\\program_sub_one.txt"
    #program_name = "programs\\program_sub_halt.txt"
    #program_name = "programs\\program_add_halt.txt"
    #program_name = "programs\\program_if.txt"
    #program_name = "programs\\program_if_false.txt"
    #program_name = "programs\\program_add_multi.txt"
    program_name = "programs\\program_while.txt"
    #
    # program_name = "programs\\program_if2.txt"


    plot_multi = False

    prog_alpha = 10
    prog_delta = 0.1#0.01
    prog_n = 2
    prog_Kd = 10
    params_prog  = prog_alpha, prog_delta, prog_n, prog_Kd


    # proteolysis and induction of protease (conditional jumps)
    deltaE = 250
    KM = 100
    KD_cond = 0.1
    params_proteolysis = [deltaE, KM]
    params_condition = [KD_cond]

    # other params - outputs of GA
    points = np.loadtxt('selected_points.txt')

    

    if plot_multi:

        ax1=plt.subplot(1, 3, 1)
        ax2=plt.subplot(1, 3, 2, sharey = ax1)
        ax3=plt.subplot(1, 3, 3, sharey = ax1)
        axes = [ax1, ax2, ax3]

        
        
        

        for i,ax in enumerate(axes):
            
            params = points[i]

            params_ff = list(params[:8])
            params_addr = list(params[8:])   

            simulate_program(program_name, t_end, N, params_ff, params_addr, params_proteolysis, params_condition, params_prog, n_bits, ax)
        plt.gcf().set_size_inches(15,5)
        plt.savefig("figs\\"+program_name+".pdf", bbox_inches = 'tight')
        plt.show()  
    
    else:
        #params = np.array([15.58143736,  2.52191547, 24.35646122, 36.12683424,  0.20222387, 0.49614202, 17.01689959,  5.        , 44.03567633,  0.73364135, 34.36351963,  5.        ])
        params = points[0]
        params_ff = list(params[:8])
        params_addr = list(params[8:])   

        simulate_program(program_name, t_end, N, params_ff, params_addr, params_proteolysis, params_condition, params_prog, n_bits)

    