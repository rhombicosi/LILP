import pandas as pd
from utils.sol_converter import *
from utils.prepro_run import *
from utils.constants_paths import *
from optimize import *

def add_column(df, column_name, values):
    df[column_name] = values
    return df

ref_mfe_files = get_filenames(efn2_archive_dir, '.txt')
ref_mfe_names = [os.path.splitext(os.path.basename(f))[0] for f in ref_mfe_files]
results_df = pd.DataFrame(index=ref_mfe_names)

ref_MFEs = [get_energy_from_ct_file(f) for f in ref_mfe_files]

rna_mfe_names = get_filenames(ct_rnastruct_dir, '.ct')
rna_MFEs = [get_energy_from_ct_file(f) for f in rna_mfe_names]

vienna_mfe_names = get_filenames(dot_bracket_viennaRNA_dir, '_enrg.txt')
vienna_MFEs = [float(open(f).read().strip()) for f in vienna_mfe_names] 

unafold_mfe_names = get_filenames(unafold_fold_dir, '_enrg.txt')
unafold_MFEs = [float(open(f).read().strip()) for f in unafold_mfe_names]

add_column(results_df, 'MFE_ref', ref_MFEs)
add_column(results_df, 'MFE_rna', rna_MFEs)
add_column(results_df, 'MFE_vienna', vienna_MFEs)
add_column(results_df, 'MFE_unafold', unafold_MFEs)
print(results_df)

n1 = 0
n2 = 1 #len(seq_files)

for seq_no in range (n1, n2):

    chain_file = seq_files[seq_no]
    chain_name_with_ext = os.path.basename(chain_file)        
    chain_name_without_ext = os.path.splitext(chain_name_with_ext)[0]
    lp_file_name = chain_name_without_ext
    seq_data = parse_seq_file(chain_file)
    rna = seq_data['sequence'].upper()

    first = 1
    last = len(rna)

    print(chain_name_without_ext)
    print(rna)
    print(len(rna))

    ###### START SOL OPTIMIZATION #######    
    # start_name = 'lilp-branch-init'
    # stem = True
    # hairpin = True
    # internal = True
    # bulge = True
    # branch = False 
    # start = False
    
    # cut_MFE, lp_name, opt_time = optimize_lilp(rna, lp_file_name, start_name, stem, hairpin, internal, bulge, branch, lpstart_dir, incumbent_dir, solstart_dir, first, last)

    # f1_gen, fbeta_gen, MCC_gen, f1_rnastruct, fbeta_rnastruct, MCC_rnastruct, f1_vienna, fbeta_vienna, MCC_vienna, f1_unafold, fbeta_unafold, MCC_unafold = sol_analyse(seq_files, seq_no, solstart_dir, start_name, dot_bracket_dir, dot_bracket_archive_dir, dot_bracket_rnastructure_dir, dot_bracket_viennaRNA_dir, unafold_fold_dir, 0)

    ###### MAIN SOL OPTIMIZATION #######
    model_name = 'lilp-branch'#'lilp-multi'#'lilp-nomulti'#'lilp-nobranch'   
    stem = True
    hairpin = True
    internal = True
    bulge = True
    multi = False
    branch = True
    start = False 

    #### NO START VERSION ####
    gen_MFE, lp_name, opt_time = optimize_lilp(rna, lp_file_name, model_name, stem, hairpin, internal, bulge, branch, lp_dir, incumbent_dir, sol_dir, first, last)
    #### WITH START VERSION ####
    # gen_MFE, lp_name, opt_time = optimize_lilp(rna, lp_file_name, model_name, stem, hairpin, internal, bulge, branch, lp_dir, incumbent_dir, sol_dir, first, last, start, start_name, solstart_dir)
   
    f1_gen, fbeta_gen, MCC_gen, f1_rnastruct, fbeta_rnastruct, MCC_rnastruct, f1_vienna, fbeta_vienna, MCC_vienna, f1_unafold, fbeta_unafold, MCC_unafold = sol_analyse(seq_files, seq_no, sol_dir, model_name, dot_bracket_dir, dot_bracket_archive_dir, dot_bracket_rnastructure_dir, dot_bracket_viennaRNA_dir, unafold_fold_dir, start)

    write_results_to_file(lp_name, len(rna), opt_time, gen_MFE/100, ref_MFEs[seq_no], rna_MFEs[seq_no], vienna_MFEs[seq_no], unafold_MFEs[seq_no], round(f1_gen,2), round(f1_rnastruct,2), round(f1_vienna,2), round(f1_unafold,2), round(fbeta_gen,2), round(fbeta_rnastruct,2), round(fbeta_vienna,2), round(fbeta_unafold,2), round(MCC_gen,2), round(MCC_rnastruct,2), round(MCC_vienna,2), round(MCC_unafold,2), results_dir, f'LILP_{len_start}_{len_end}_{model_name}.txt')