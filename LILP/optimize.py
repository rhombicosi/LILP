import time
from utils.constants_paths import *
from utils.prepro_run import *
from lilp import *

def make_callback():
    best_obj = float("inf")
    best_obj_time = None

    def callback(model, where):
        nonlocal best_obj, best_obj_time

        if where == GRB.Callback.MIPSOL:
            obj = model.cbGet(GRB.Callback.MIPSOL_OBJ)

            if obj < best_obj:
                best_obj = obj
                best_obj_time = model.cbGet(GRB.Callback.RUNTIME)

    return callback, lambda: (best_obj, best_obj_time)

def optimize_lilp(rna: str, lp_file_name: str, model_name: str, stem: bool, hairpin: bool, internal: bool, bulge: bool, branch: bool, cbranch: bool, lp_dir: str, incumbent_dir: str, sol_dir: str, first = None, last = None, start = None, start_name = None, solstart_dir = None) -> None:
    
    model_start_time = time.time()
    rna_model = LILP(rna, model_name)
    rna_model.create_variables(stem, hairpin, internal, bulge, branch, cbranch, first, last)
    rna_model.create_constraints(stem, hairpin, internal, bulge, branch, cbranch, first, last)
    rna_model.create_objective(stem, hairpin, internal, bulge, branch, cbranch)

    # rna_model.model.addConstr(rna_model.model.getVarByName(f'P_7_47') == 1)
    # rna_model.model.addConstr(rna_model.model.getVarByName(f'P_24_40') == 1)
    # rna_model.model.addConstr(rna_model.model.getVarByName(f'X_41') == 1)
    # rna_model.model.addConstr(rna_model.model.getVarByName(f'X_42') == 1)
    # rna_model.model.addConstr(rna_model.model.getVarByName(f'X_43') == 1)
    # rna_model.model.addConstr(rna_model.model.getVarByName(f'X_44') == 1)
    # rna_model.model.addConstr(rna_model.model.getVarByName(f'X_45') == 1)
    # rna_model.model.addConstr(rna_model.model.getVarByName(f'X_46') == 1)
    # rna_model.model.addConstr(rna_model.model.getVarByName(f'Y_7_47_24_40') == 1)      
    # rna_model.model.addConstr(rna_model.model.getVarByName(f'CBRANCH_7_47_24_40') == 1)

    model_time = time.time() - model_start_time
    print(f'MODEL CONSTRUCTION TIME :: {model_time}')
    rna_model.model.write(f'{lp_dir}/{lp_file_name}-{model_name}.lp')

    rna_model.model.setParam("LogFile", f'{grb_log_dir}/{lp_file_name}-log-{model_name}')
    rna_model.model.setParam(GRB.Param.SolFiles, f'{incumbent_dir}/{lp_file_name}-incumbent-{model_name}') 
    
    # rna_model.model.setParam("Cuts", 3)
    # rna_model.model.setParam("BQPCuts", 2)
    # rna_model.model.setParam("PSDCuts", 2)
    # rna_model.model.setParam("ModKCuts", 2)    
    # rna_model.model.setParam("CliqueCuts", 1)
    # rna_model.model.setParam("RLTCuts", 1)
    # rna_model.model.setParam("ZeroHalfCuts", 1)
    # rna_model.model.setParam("RelaxLiftCuts", 2)
    # rna_model.model.setParam("MIPFocus", 3)
    # rna_model.model.setParam("Heuristics", 0)
    # rna_model.model.setParam('Cuts', 2)        # aggressiveness 0-3
    # rna_model.model.setParam('GomoryPasses', 5) # explicit Gomory cut rounds
    # rna_model.model.setParam('CoverCuts', 2)    # knapsack cover aggressiveness
    # rna_model.model.setParam('MIRCuts', 2)      # MIR cut aggressiveness
    # rna_model.model.setParam('CliqueCuts', 2)   # clique cut aggressiveness
    rna_model.model.setParam("TimeLimit", 2400)
    rna_model.model.setParam("MIPGap", 0.002)
    rna_model.model.setParam("Threads", 24)
    rna_model.model.setParam("NodefileStart", 0.5)  # start disk swapping earlier
    
    if start:
        rna_model.model.NumStart = 1
        rna_model.model.update()
        solvars = read_sol(f'{solstart_dir}/{lp_file_name}-{start_name}.sol')

        # start values
        for v in rna_model.model.getVars(): 
            if v.VarName in solvars.keys():
                v.Start = round(int(solvars[v.VarName]), 1)
        rna_model.model.update()

    opt_start_time = time.time()
    callback, get_results = make_callback()
    rna_model.model.optimize(callback)
    best_obj, best_obj_time = get_results()
    opt_time = time.time() - opt_start_time
    print(f'OPTIMIZATION TIME :: {opt_time}')

    rna_model.model.write(f'{sol_dir}/{lp_file_name}-{model_name}.sol')

    print(f'Last incumbent time: {best_obj_time}')
    print(f'Gap: {rna_model.model.MIPGap}')
    print(f'Obj: {rna_model.model.ObjVal:g}')

    if rna_model.model.ObjVal is not None:
        return rna_model.model.ObjVal, lp_file_name, opt_time, rna_model.model.MIPGap, best_obj_time
    else:
        print("Object value was not assigned due to an error.")

# seq_number = 7
# #46/39/34#19#29#18

# chain_file = seq_files[seq_number]
# chain_name_with_ext = os.path.basename(chain_file)
# chain_name_without_ext = os.path.splitext(chain_name_with_ext)[0]
# lp_file_name = chain_name_without_ext
# seq_data = parse_seq_file(chain_file)

# rna = seq_data['sequence']
# print(chain_file)
# print(rna)
# print(len(rna))

# start_name = 'lilp-H1start'
# stem = True
# hairpin = True
# internal = False
# bulge = False
# multi = False
# start = False
# optimize_lilp(rna, start_name, stem, hairpin, internal, bulge, multi, lpstart_dir, incumbent_start_dir, solstart_dir, start)
# model_name = 'lilp'
# stem = True
# hairpin = True
# internal = True
# bulge = True
# multi = True
# start = True
# optimize_lilp(rna, model_name, stem, hairpin, internal, bulge, multi, lp_dir, incumbent_dir, sol_dir)
# optimize_lilp(rna, model_name, stem, hairpin, internal, bulge, multi, lp_dir, incumbent_dir, sol_dir, start, start_name, solstart_dir)

# process solution(s) to dot-bracket
# filepath = f'{sol_dir}/{lp_file_name}-{model_name}.sol'
# pairs2brackets(filepath, rna)

# calculate_sol_energy(filepath, rna)

# script_dir = os.path.dirname(os.path.abspath(__file__))
# parent_dir = os.path.dirname(script_dir)

# for f in range(8):
#     filepath = f'{incumbent_dir}\lilp_{seq_number}_incumbent_{f}.sol'
#     pairs2brackets(filepath, rna)