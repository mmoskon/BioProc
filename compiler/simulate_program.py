from bioproc.hill_functions import *
from generate_model import generate_model
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np
import importlib
import seaborn as sns
import pickle


def simulate_program(program_name, t_end, N, params_ff, params_addr, params_proteolysis, params_condition, params_prog, n_bits, ax=plt, plot_clock = True, plot_instructions = False, plot_ops = True, alpha_plot = 0.75):
    
    model_name = program_name.split(".")[0]
    
    prog_alpha, prog_delta, prog_n, prog_Kd = params_prog 


    generate_model(program_name, model_name, n_bits, prog_alpha, prog_delta, prog_n, prog_Kd)      
    model = importlib.import_module(model_name.replace("/","."))

    f_description = open(model_name+"description.txt")
    ops = f_description.readline().strip().split(",")[:-1]
    f_description.close()

    # params
    params = (params_ff + params_proteolysis + params_addr + params_condition,)
    
    # solving
    if n_bits == 3:
        Y0 = np.array([0]*(18+len(ops)))
    elif n_bits == 4:
        Y0 = np.array([0]*(24+len(ops)))
    else:
        Y0 = np.array([0]*(30+len(ops)))

    T = np.linspace(0, t_end, N)
    Y = odeint(model.model, Y0, T, args=params)

    legend = []

        
    if plot_instructions:
        if n_bits == 3:
            i1 = Y[:,-6-len(ops)]
            i2 = Y[:,-5-len(ops)]
            i3 = Y[:,-4-len(ops)]
            i4 = Y[:,-3-len(ops)]
            i5 = Y[:,-2-len(ops)]
            i6 = Y[:,-1-len(ops)]
        elif n_bits == 4:
            i1 = Y[:,-8-len(ops)]
            i2 = Y[:,-7-len(ops)]
            i3 = Y[:,-6-len(ops)]
            i4 = Y[:,-5-len(ops)]
            i5 = Y[:,-4-len(ops)]
            i6 = Y[:,-3-len(ops)]
            i7 = Y[:,-2-len(ops)]
            i8 = Y[:,-1-len(ops)]
        else: #elif n_bits == 5:
            i1 = Y[:,-10-len(ops)]
            i2 = Y[:,-9-len(ops)]
            i3 = Y[:,-8-len(ops)]
            i4 = Y[:,-7-len(ops)]
            i5 = Y[:,-6-len(ops)]
            i6 = Y[:,-5-len(ops)]
            i7 = Y[:,-4-len(ops)]
            i8 = Y[:,-3-len(ops)]
            i9 = Y[:,-2-len(ops)]
            i10 = Y[:,-1-len(ops)]

        ax.plot(T,i1, alpha=alpha_plot)
        ax.plot(T,i2, alpha=alpha_plot)
        ax.plot(T,i3, alpha=alpha_plot)
        ax.plot(T,i4, alpha=alpha_plot)
        ax.plot(T,i5, alpha=alpha_plot)
        ax.plot(T,i6, alpha=alpha_plot)
        if n_bits == 4:
            ax.plot(T,i7, alpha=alpha_plot)
            ax.plot(T,i8, alpha=alpha_plot)
        elif n_bits == 5:
            ax.plot(T,i7, alpha=alpha_plot)
            ax.plot(T,i8, alpha=alpha_plot)
            ax.plot(T,i9, alpha=alpha_plot)
            ax.plot(T,i10, alpha=alpha_plot)
    
    
        
        if n_bits == 3:
            ax.legend(ops+['i1','i2','i3','i4','i5','i6'], loc='upper left')
            legend += ['i1','i2','i3','i4','i5','i6']
        elif n_bits == 4: #else:
            ax.legend(ops+['i1','i2','i3','i4','i5','i6','i7','i8'], loc='upper left')
            legend += ['i1','i2','i3','i4','i5','i6','i7','i8']
        elif n_bits == 5: #else:
            ax.legend(ops+['i1','i2','i3','i4','i5','i6','i7','i8','i9','i10'], loc='upper left')
            legend += ['i1','i2','i3','i4','i5','i6','i7','i8','i9','i10']


    if plot_ops:
        i = -len(ops)
        for op in ops:
            o = Y[:,i]
            ax.plot(T,o, alpha=alpha_plot)
            i += 1

        legend += ops
    
    
    if plot_clock:
        clk = get_clock(T)
        ax.plot(T,90*clk, color="black", alpha=0.1)
        ax.legend(ops + ['clk'], loc='upper left')
        legend += ['clk']
    
    
    #if ax != plt:
    #    ax.set_xlabel('Time [h]')
    #    ax.set_ylabel('Concentrations [nM]')

    ax.legend(legend, ncol=10, 
          loc='upper center',
          bbox_to_anchor=(0.5, 0.95),
          bbox_transform=plt.gcf().transFigure)

    if ax == plt:
        plt.savefig("figs/"+program_name.split(".")[0]+".pdf", bbox_inches = 'tight')
        plt.savefig("figs/"+program_name.split(".")[0]+".png")
        
        plt.xlabel('Time [h]')
        plt.ylabel('Concentrations [nM]')

        plt.show()




if __name__ == '__main__':

    # simulation parameters

    t_end = 290 #200
    N = 1200 #1000

    n_bits = 5

    plot_ops = True

    prog_alpha = 1000
    prog_delta = 0.1#0.01
    prog_n = 2
    prog_Kd = 10
    


    #program_name = "programs/program_add.txt"
    #program_name = "programs\\program_add_nop.txt"
    #program_name = "programs\\program_sub.txt"
    #program_name = "programs\\program_sub_one.txt"
    #program_name = "programs\\program_sub_halt.txt"
    #program_name = "programs\\program_add_halt.txt"
    #program_name = "programs\\program_if.txt"
    program_name = "programs/program_if_false.txt"
    #program_name = "programs\\program_add_multi.txt"
    #program_name = "programs\\program_while.txt"
    #
    #program_name = "programs\\program_if.txt"
    
    ### Figures for the paper ###
    #program_name, t_end = "programs\\Figure_nop.txt", 150
    #program_name = "programs\\Figure_halt.txt"
    #program_name, t_end = "programs\\Figure_jump_unconditional.txt", 120
    #program_name, t_end = "programs\\Figure_jump_conditional_false.txt", 120
    #program_name, t_end, prog_delta = "programs\\Figure_jump_conditional_true.txt", 120, 0.1
    #program_name, t_end, prog_delta = "programs\\Figure_if_true.txt", 120, 0.1
    #program_name, t_end, prog_delta = "programs\\Figure_if_false.txt", 120, 0.1
    #program_name, t_end = "programs/Figure_while.txt", 180

    plot_multi = False #True

    params_prog  = prog_alpha, prog_delta, prog_n, prog_Kd


    # proteolysis and induction of protease (conditional jumps)
    deltaE = 250
    KM = 100
    KD_cond = 0.1
    params_proteolysis = [deltaE, KM]
    params_condition = [KD_cond]

    # other params - outputs of GA

    #points = np.loadtxt('selected_points//selected_points_old.txt')[:3]
    p31 = pickle.load(open("selected_points/model3sample1.p", "rb"))      
    p32 = pickle.load(open("selected_points/model3sample2.p", "rb"))      
    p33 = pickle.load(open("selected_points/model3sample3.p", "rb"))      
    points = np.array([p31, p32, p33])
    
    #with plt.style.context('fivethirtyeight'):
    sns.set_style("white")

    if plot_multi:

        ax1=plt.subplot(1, 3, 1)
        ax2=plt.subplot(1, 3, 2, sharey = ax1)
        ax3=plt.subplot(1, 3, 3, sharey = ax1)
        axes = [ax1, ax2, ax3]

        
        
        

        for i,ax in enumerate(axes):
            
            params = points[i]

            params_ff = list(params[:8])
            params_addr = list(params[8:])   

            simulate_program(program_name, t_end, N, params_ff, params_addr, params_proteolysis, params_condition, params_prog, n_bits, ax, plot_instructions=True, plot_ops = plot_ops)

            if i == 0:
                ax.set_ylabel('Concentrations [nM]')  
            ax.set_xlabel('Time [h]')


        plt.gcf().set_size_inches(15,5)
        plt.savefig("figs/"+program_name.split(".")[0]+".pdf", bbox_inches = 'tight')
        plt.show()  
    
    else:
        #params = np.array([15.58143736,  2.52191547, 24.35646122, 36.12683424,  0.20222387, 0.49614202, 17.01689959,  5.        , 44.03567633,  0.73364135, 34.36351963,  5.        ])
        params = points[0]
        params_ff = list(params[:8])
        params_addr = list(params[8:])   

        simulate_program(program_name, t_end, N, params_ff, params_addr, params_proteolysis, params_condition, params_prog, n_bits)

    