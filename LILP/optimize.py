import time
from utils.constants_paths import *
from utils.prepro_run import *
from lilp import *

def cut_callback(model, where):
    if where == GRB.Callback.MIPNODE:
        # cuts found at node
        cuts = model.cbGet(GRB.Callback.MIP_CUTCNT)
        obj_best = model.cbGet(GRB.Callback.MIPNODE_OBJBST)
        # You can print or save this info
        print(f'Node cuts: {cuts}, Best Obj: {obj_best}')
    elif where == GRB.Callback.MIP:
        cuts = model.cbGet(GRB.Callback.MIP_CUTCNT)
        obj_best = model.cbGet(GRB.Callback.MIP_OBJBST)
        obj_bound = model.cbGet(GRB.Callback.MIP_OBJBND)
        print(f'Cuts used so far: {cuts}, At MIP callback: Best={obj_best}, Bound={obj_bound}')

def optimize_lilp(rna: str, lp_file_name: str, model_name: str, stem: bool, hairpin: bool, internal: bool, bulge: bool, branch: bool, lp_dir: str, incumbent_dir: str, sol_dir: str, first = None, last = None, start = None, start_name = None, solstart_dir = None) -> None:
    
    model_start_time = time.time()
    rna_model = LILP(rna, model_name)
    rna_model.create_variables(stem, hairpin, internal, bulge, branch, first, last)
    rna_model.create_constraints(stem, hairpin, internal, bulge, branch, first, last)
    rna_model.create_objective(stem, hairpin, internal, bulge, branch)

    # rna_model.create_cut(stem, hairpin, internal, bulge, 0, -1523)
      
    # rna_model.model.addConstr(rna_model.model.getVarByName(f'STEM_5_68_6_67') == 1)

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
    rna_model.model.setParam("TimeLimit", 3600)
    rna_model.model.setParam("MIPGap", 0.002)
    rna_model.model.setParam("Threads", 32)
    
    # sorted_bp = sorted(rna_model.base_pairs, key=lambda x: x.distance)
    # for bp in sorted_bp:
    #     bp.var.setAttr("BranchPriority", 10 * bp.distance)

    # if hairpin:
    #     for hl in rna_model.hairpin_loops:
    #         if hl.is_valid_size():
    #             hl.var.setAttr("BranchPriority", round(200/hl.size))
    #         else:
    #             hl.var.setAttr("BranchPriority", 0)
    
    # if stem:
    #     # sorted_stems = sorted(rna_model.stem_loops, key=lambda x: x.distance)   
    #     for sl in rna_model.stem_loops:
    #         sl.var.setAttr("BranchPriority", 10 * sl.distance)
  
    # if internal:
    #     for il in rna_model.internal_loops:
    #         il.var.setAttr("BranchPriority", round(2000/(il.energy+1)))    

    # if bulge:
    #     sorted_bulges = sorted(rna_model.bulge_loops, key=lambda x: x.size)
    #     for bl in sorted_bulges:
    #         bl.var.setAttr("BranchPriority", round(100/bl.size))

    # if multi:
    #     for ml in rna_model.multi_loops:
    #         if ml.size == 0:
    #             ml.var.setAttr("BranchPriority", round(100/(ml.size + 1)))
    #         else:
    #             ml.var.setAttr("BranchPriority", round(100/ml.size))

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
    rna_model.model.optimize()
    opt_time = time.time() - opt_start_time
    print(f'OPTIMIZATION TIME :: {opt_time}')

    rna_model.model.write(f'{sol_dir}/{lp_file_name}-{model_name}.sol')

    print(f'Obj: {rna_model.model.ObjVal:g}')

    if rna_model.model.ObjVal is not None:
        return rna_model.model.ObjVal, lp_file_name, opt_time
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