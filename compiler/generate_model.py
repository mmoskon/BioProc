from collections import defaultdict
import numpy as np

def from_addr_to_i(addr, n_bits = 3):
    if n_bits == 3:
        if addr == [0,0,0]:
            return "i1"
        elif addr == [1,0,0]:
            return "i2"
        elif addr == [1,1,0]:
            return "i3"
        elif addr == [1,1,1]:
            return "i4"
        elif addr == [0,1,1]:
            return "i5"
        else:
            return "i6"
    elif n_bits == 4:
        if addr == [0,0,0,0]:
            return "i1"
        elif addr == [1,0,0,0]:
            return "i2"
        elif addr == [1,1,0,0]:
            return "i3"
        elif addr == [1,1,1,0]:
            return "i4"
        elif addr == [1,1,1,1]:
            return "i5"
        elif addr == [0,1,1,1]:
            return "i6"
        elif addr == [0,0,1,1]:
            return "i7"
        else:
            return "i8"
    elif n_bits == 5:
        if addr == [0,0,0,0,0]:
            return "i1"
        elif addr == [1,0,0,0,0]:
            return "i2"
        elif addr == [1,1,0,0,0]:
            return "i3"
        elif addr == [1,1,1,0,0]:
            return "i4"
        elif addr == [1,1,1,1,0]:
            return "i5"
        elif addr == [1,1,1,1,1]:
            return "i6"
        elif addr == [0,1,1,1,1]:
            return "i7"
        elif addr == [0,0,1,1,1]:
            return "i8"
        elif addr == [0,0,0,1,1]:
            return "i9"
        else:
            return "i10"

    

def from_i_to_addr(i, n_bits = 3):
    if n_bits == 3:
        if i == "i1":
            return [0,0,0]
        elif i == "i2":
            return [1,0,0]
        elif i == "i3":
            return [1,1,0]           
        elif i == "i4":
            return [1,1,1]            
        elif i == "i5":
            return [0,1,1]
        else:
            return [0,0,1]
    elif n_bits == 4:
        if i == "i1":
            return [0,0,0,0]
        elif i == "i2":
            return [1,0,0,0]
        elif i == "i3":
            return [1,1,0,0]           
        elif i == "i4":
            return [1,1,1,0]            
        elif i == "i5":
            return [1,1,1,1]
        elif i == "i6":
            return [0,1,1,1]
        elif i == "i7":
            return [0,0,1,1]
        else:
            return [0,0,0,1]
    elif n_bits == 5:
        if i == "i1":
            return [0,0,0,0,0]
        elif i == "i2":
            return [1,0,0,0,0]
        elif i == "i3":
            return [1,1,0,0,0]           
        elif i == "i4":
            return [1,1,1,0,0]            
        elif i == "i5":
            return [1,1,1,1,0]
        elif i == "i6":
            return [1,1,1,1,1]
        elif i == "i7":
            return [0,1,1,1,1]
        elif i == "i8":
            return [0,0,1,1,1]
        elif i == "i9":
            return [0,0,0,1,1]
        else:
            return [0,0,0,0,1]


def from_addr_to_RS(jump_src, jump_dst, n_bits = 3): 
    i_src = from_addr_to_i(jump_src, n_bits)
    i_dst = from_addr_to_i(jump_dst, n_bits)    
    #print(i_src)
    #print(i_dst)
        
    R = [0]*n_bits
    S = [0]*n_bits

    for i in range(len(jump_src)):
        if jump_src[i] > jump_dst[i]:
                R[i] = i_src
        elif jump_src[i] < jump_dst[i]:
                S[i] = i_src
        
    
    return R,S

def from_i_to_RS(i_src, i_dst, n_bits = 3):
    
    jump_src = from_i_to_addr(i_src, n_bits)
    jump_dst = from_i_to_addr(i_dst, n_bits)
    #print(i_src, jump_src)
    #print(i_dst, jump_dst)
        
        
    R = [0]*n_bits
    S = [0]*n_bits

    for i in range(len(jump_src)):
        if jump_src[i] > jump_dst[i]:
                R[i] = i_src
        elif jump_src[i] < jump_dst[i]:
                S[i] = i_src
        
    
    return R,S




####################################
####################################
####################################    

def generate_model(program_name, output_name, n_bits, prog_alpha, prog_delta, prog_n, prog_Kd):
    code = []
    variables = ""
    pars = ""


    R = []
    S = []
    condition = []
    inhibition = []

    ## clock
    code.append("\tclk = get_clock(T)\n")



    ## program
    f_prog = open(program_name)

    addr = 1
    prog = defaultdict(str)
    prog_params = set()
    operands = set()


    for line in f_prog:
        for l in line.strip().split(";"):        
            try:
                l = l.split(",")
                ops = l[0].split()[1:] + l[1:]    
                instr = (l[0].split()[0]).lower().strip()   
            except:
                instr = l[0].strip()
        
            if instr == "if" or instr=="do-while":
                
                if instr == "if":
                    i_src = "i"+str(addr)
                    i_dst = "i"+str(addr+1)
                    inhibition.append(True)
                else: # do-while
                    i_src = "i"+str(addr + 1)
                    i_dst = "i"+str(addr)
                    inhibition.append(False)
                
                r,s  = from_i_to_RS(i_src, i_dst, n_bits = n_bits)
                R.append(r)
                S.append(s)
                
                condition.append(ops[0])
                operands.add(ops[0])
                instr = ops[1]    
                ops = ops[2:]
            elif instr == "while":
                # preverim, ce je izpolnjen pogoj, ce ni skocim naprej
                i_src = "i"+str(addr)
                i_dst = "i"+str(addr+1)
                inhibition.append(True)
                r,s  = from_i_to_RS(i_src, i_dst, n_bits = n_bits)
                R.append(r)
                S.append(s)
                condition.append(ops[0])
            
                # brezpogojno skocim nazaj, kjer bom spet preveril pogoj
                # lahko pa tudi pogojno skocim nazaj
                # kjer brez dodatnih stroskov ponovno preverjam enak pogoj
                i_src = "i"+str(addr + 1)
                i_dst = "i"+str(addr)
                inhibition.append(False)
                r,s  = from_i_to_RS(i_src, i_dst, n_bits = n_bits)
                R.append(r)
                S.append(s)
                condition.append(ops[0])
                
                operands.add(ops[0])
                instr = ops[1]    
                ops = ops[2:]
            
            if instr == 'nop':
                pass
            elif instr == 'generate1':
                o = ops[0] + '0'
                operands.add(o)
                operands.add(ops[0] + '1')
                
                alpha = "prog_alpha_" + o 
                Kd = "prog_Kd_" + o
                n = "prog_n_" + o
                prog_params |= {alpha, Kd, n}    
                prog[o] += "+"+alpha+"*activate_1(i"+str(addr)+",prog_Kd_"+o+",prog_n_"+o+")"
            
            elif instr == 'generate2':
                o = ops[0] + '1'
                operands.add(o)
                operands.add(ops[0] + '0')

                alpha = "prog_alpha_" + o 
                Kd = "prog_Kd_" + o
                n = "prog_n_" + o
                
                prog_params |= {alpha, Kd, n}
                
                prog[o] += "+"+alpha+"*activate_1(i"+str(addr)+",prog_Kd_"+o+",prog_n_"+o+") + 10"
            elif instr == 'generate3':
                o = ops[0] + '0'
                operands.add(o)

                alpha = "prog_alpha_" + o 
                Kd = "prog_Kd_" + o
                n = "prog_n_" + o
                prog_params |= {alpha, Kd, n}
                prog[o] += "+"+alpha+"*activate_1(i"+str(addr)+",prog_Kd_"+o+",prog_n_"+o+")"
                # 1st bit -----
                o = ops[0] + '1'
                operands.add(o)
                
                alpha = "prog_alpha_" + o 
                Kd = "prog_Kd_" + o
                n = "prog_n_" + o
                prog_params |= {alpha, Kd, n}
                prog[o] += "+"+alpha+"*activate_1(i"+str(addr)+",prog_Kd_"+o+",prog_n_"+o+") + 10"


            elif instr == 'add' or instr == 'sub':
                o = ops[0].strip()
                operands.add(o + '0')
                operands.add(o + '1')
                op1 = ops[1].strip()
                operands.add(op1 + '0')
                operands.add(op1 + '1')
                op2 = ops[2].strip()
                operands.add(op2 + '0')
                operands.add(op2 + '1')
                
                alpha = "prog_alpha_" + o + '0'
                Kd = "prog_Kd_" + o + '0'
                n = "prog_n_" + o + '0'
                
                alpha1 = "prog_alpha_" + o + '1' 
                Kd1 = "prog_Kd_" + o + '1'
                n1 = "prog_n_" + o + '1'

                prog_params |= {alpha, Kd, n, alpha1, Kd1, n1}
                if instr == 'add':
                    #prog[o] += alpha+"*activate_3(i"+str(addr)+","+op1+','+op2+",prog_Kd_"+o+",prog_n_"+o+")"
                    # 0-th bit
                    o = ops[0].strip() + '0'
                    #prog[o] += "+"+alpha+"*activate_2(i"+str(addr)+","+op1+"0,prog_Kd_"+o+",prog_n_"+o+")"
                    #prog[o] += "+"+alpha+"*activate_2(i"+str(addr)+","+op2+"0,prog_Kd_"+o+",prog_n_"+o+")"
                    #prog[o] += "+"+alpha+"*activate_XOR_3(i"+str(addr)+","+op1+"0,"+op2+"0,prog_Kd_"+o+",prog_n_"+o+")"
                    prog[o] += "+"+alpha+"*zeroth_bit(i"+str(addr)+","+op1+"0,"+op2+"0,prog_Kd_"+o+",prog_n_"+o+")"
                    #prog[o] += "+"+alpha+"*activate_XOR_2(i"+str(addr)+","+op2+"0,prog_Kd_"+o+",prog_n_"+o+")"
                    
                    # 1-st bit
                    o = ops[0].strip() + '1'
                    #prog[o] += "+"+alpha1+"*activate_2(i"+str(addr)+","+op1+"1,prog_Kd_"+o+",prog_n_"+o+")"
                    #prog[o] += "+"+alpha1+"*activate_2(i"+str(addr)+","+op2+"1,prog_Kd_"+o+",prog_n_"+o+")"
                    prog[o] += "+"+alpha1+"*first_bit(i"+str(addr)+","+op1+'0,'+op1+"1,"+op2+'0,'+op2+"1,prog_Kd_"+o+",prog_n_"+o+") + 10"
                    
                    print(addr, instr, '———>', prog[o])
                else:
                    prog[o] += "+"+alpha+"*hybrid_AAR(i"+str(addr)+","+op1+','+op2+",prog_Kd_"+o+",prog_n_"+o+")"
                    print(instr, '———>', prog[o])
            elif instr == "halt":
                i_src = "i"+str(addr + 1)
                i_dst = "i"+str(addr)
                inhibition.append(True)
                condition.append(0)
                r, s = from_i_to_RS(i_src, i_dst, n_bits = n_bits)
                R.append(r)
                S.append(s)

                break
            elif instr == "jump" or instr=="jumpif":
                
                i_src = "i"+str(addr)
                i_dst = "i"+str(addr+1)

                if instr == "jump":   
                    inhibition.append(True)
                    condition.append("0")
                else: # jumpif
                    inhibition.append(False)
                    condition.append(ops[0])
                    operands.add(ops[0])
                
                r,s  = from_i_to_RS(i_src, i_dst, n_bits = n_bits)
                R.append(r)
                S.append(s)
                
                
            
            
        addr += 1
    for o in operands: # dodam razgradnjo tudi za tiste, ki sicer nimajo enacbe
        prog[o] += "-prog_delta_"+o +"*"+o
        prog_params.add("prog_delta_" + o)

    #print(prog)
    #print(prog_params)


    for p in prog_params:
        value = eval("_".join(p.split("_")[:-1]))
        code.append("\t"+p+"="+str(value)+"\n")

    variables_prog = ""   
    for op in prog:
        code.append("\td"+op+"_dt="+prog[op]+"\n")
        
    variables_list = list(prog.keys())
    variables_list.sort()
    for op in variables_list:        
        variables_prog +=op+","


    ## jumps
    
    #jump_src = ""
    #jump_dst = ""

    #jump_src=[0,1,1]
    #jump_dst = [1,1,1]

    
    #R,S = from_addr_to_RS(jump_src, jump_dst, n_bits)

    #print(R,S)
    for i in range(len(condition)):
        if condition[i]:
            code.append("\tcond"+str(i)+"="+condition[i]+"\n")



    for j in range(n_bits):
        str_R = "\tRESET"+str(j)+"="
        str_S = "\tSET"+str(j)+"="
        str_R_inputs = ""
        str_S_inputs = ""        
        for i in range(len(condition)):
            cond = condition[i]
            inh = inhibition[i]
            r = R[i][j] 
            s = S[i][j]
            if not cond:              
                if r:
                    str_R_inputs += r + ","
                if s:
                    str_S_inputs += s + ","
            else:
                if inh:
                    #inhibition
                    if r:
                        str_R_inputs += "inhibition("+r+", cond"+str(i)+", KD_cond),"
                    if s:
                        str_S_inputs += "inhibition("+s+", cond"+str(i)+", KD_cond),"
                else:
                    #induction
                    if r:
                        str_R_inputs += "induction("+r+", cond"+str(i)+", KD_cond),"
                    if s:
                        str_S_inputs += "induction("+s+", cond"+str(i)+", KD_cond),"




        if str_R_inputs:
            str_R += "max(("+str_R_inputs+"))"
        else:
            str_R += "0"
        if str_S_inputs:
            str_S += "max(("+str_S_inputs+"))"
        else:
            str_S += "0"
        
        str_R += " if T > 1 else 100\n"
        str_S += " \n"
        
        code.append(str_R)
        code.append(str_S)





    """
    
    for i in range(len(condition)):
        cond = condition[i]
        inhibit = inhibition[i]
        r = R[i]
        s = S[i]


    if condition:
        code.append("\tcond="+condition+"\n")
    for i in range(n_bits):
        if not condition:
            code.append("\tRESET"+str(i)+"="+str(R[i])+" if T > 1 else 100\n")
            code.append("\tSET"+str(i)+"="+str(S[i])+"\n")
        else:
            if inhibition:
                code.append("\tRESET"+str(i)+"=inhibition("+str(R[i])+", cond, KD_cond) if T > 1 else 100\n")
                code.append("\tSET"+str(i)+"=inhibition("+str(S[i])+", cond, KD_cond)\n")
            else:
                code.append("\tRESET"+str(i)+"=induction("+str(R[i])+", cond, KD_cond) if T > 1 else 100\n")
                code.append("\tSET"+str(i)+"=induction("+str(S[i])+", cond, KD_cond)\n")
    """


    # flip_flops
    f = open("flip_flop")
    ff_vars = f.readline().strip()
    ff_pars = f.readline().strip()
    pars += ff_pars
    ff = f.read()
    f.close()

    for i in range(n_bits):
        if i > 0:
            code.append("\td"+str(i)+"=q"+str(i-1)+"\n")
        else:
            code.append("\td"+str(i)+"=not_q"+str(n_bits-1)+"\n")

        variables += ff_vars.replace("$", str(i))
        ff1 = ff.replace("$", str(i))
        code.append(ff1)

    # addressing
    f = open("addressing"+str(n_bits))
    variables += f.readline().strip()
    addr_pars = f.readline().strip()
    pars += addr_pars
    addr = f.read()
    code.append(addr)
    f.close()

    # jump conditions
    f = open("conditions")
    cond_pars = f.readline().strip()
    pars += cond_pars
    f.close()


    # finishing...
    variables += variables_prog

    #f_m = open("model"+str(n_bits)+".py", "w")
    f_m = open(output_name+".py", "w")
    #f_m.write("from flip_flops2 import *\n")
    f_m.write("from bioproc.hill_functions import *\n")

    f_m.write("def model(Y, T, params):\n")
    f_m.write("\t"+variables+"=Y\n")
    f_m.write("\t"+pars+"=params\n")
    for c in code:
            f_m.write(c)

    dvars = variables.split(",")[:-1]
    dvars = "["+",".join(map(lambda s:'d'+s.strip()+'_dt' ,dvars))+"]"

    f_m.write("\treturn "+dvars)
    f_m.close()

    #f_d = open("model"+str(n_bits)+"description.txt", "w")
    f_d = open(output_name+"description.txt", "w")
    f_d.write(variables_prog)
    f_d.close()


####################
####################
####################

if __name__ == '__main__':
    n_bits = 3
    program_name = "programs\\program_if.txt"

    prog_alpha = 10
    prog_delta = 0.01
    prog_n = 2
    prog_Kd = 10

    generate_model(program_name, "model3", n_bits, prog_alpha, prog_delta, prog_n, prog_Kd)

