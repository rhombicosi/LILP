import os
# from pathlib import Path
from lilp import *
# from utils.prepro_run import *
from utils.sol_converter import *
from start_solution import *

#len_start = 60
seq_number = 37

# cwd = Path.cwd()
# code_path = Path(__file__).parent.parent
# arch_rel_path = '../../ARCHIVE II/'
# archive_path = (code_path/arch_rel_path).resolve()
# # seq_len_dir = f'RNA_seq_{seq_len}'

# chain_dir = os.path.join(archive_path, seq_len_dir)
# seq_files = get_filenames(chain_dir, '.seq')

chain_file = seq_files[seq_number]
chain_name_with_ext = os.path.basename(chain_file)   
print(chain_name_with_ext)     
chain_name_without_ext = os.path.splitext(chain_name_with_ext)[0]
lp_file_name = chain_name_without_ext
seq_data = parse_seq_file(chain_file)

rna = seq_data['sequence'].upper()
print(rna)


model_name = 'cb-xyz'
start_name = 'lilp-branch-init'
# for f in range(14,15):
#     filepath = f'{incumbent_dir}/{lp_file_name}-incumbent-{model_name}_{f}.sol'
#     fold, pairs = pairs2brackets(filepath, rna)
#     print(fold)
#     calculate_sol_energy(filepath, rna)


# filepath = f'{sol_dir}/{lp_file_name}-{model_name}.sol'
# filepath = os.path.join(solstart_dir, f'start.sol')
# fold, pairs = pairs2brackets(filepath, rna)
# print(fold)
# calculate_sol_energy(filepath, rna)


# chain_file = seq_files[seq_number]
# chain_name_with_ext = os.path.basename(chain_file)        
# chain_name_without_ext = os.path.splitext(chain_name_with_ext)[0]
# lp_file_name = chain_name_without_ext

# seq_data = parse_seq_file(chain_file)
# this_RNA = seq_data['sequence']    

# filepath = os.path.join(sol_dir, f'{lp_file_name}-{model_name}.sol')



# archive referece
dotbracket_ref = f'{dot_bracket_archive_dir}/{lp_file_name}.txt'
sol_ref = f'{solstart_dir}/{lp_file_name}-ref.sol'

# rnafold solution
dotbracket_vienna = f'{dot_bracket_viennaRNA_dir}/{lp_file_name}_db.txt'
sol_vienna = f'{solstart_dir}/{lp_file_name}-vienna.sol'

# rnastructure solution
dotbracket_rnastruct = f'{dot_bracket_rnastructure_dir}/{lp_file_name}.txt'
sol_rnastruct = f'{solstart_dir}/{lp_file_name}-rnastruct.sol'

# unafold solution
dotbracket_una = f'{unafold_fold_dir}/{lp_file_name}_db.txt'
sol_una = f'{solstart_dir}/{lp_file_name}-una.sol'

print("\nREFERENCE\n")
generate_start_sol(dotbracket_ref, sol_ref)
fold, pairs = pairs2brackets(sol_ref, rna)
calculate_sol_energy(sol_ref, rna)

print("\nRNAFOLD\n")
generate_start_sol(dotbracket_vienna, sol_vienna)
fold, pairs = pairs2brackets(sol_vienna, rna)
calculate_sol_energy(sol_vienna, rna)

print("\nRNASTRUCTURE\n")
generate_start_sol(dotbracket_rnastruct, sol_rnastruct)
fold, pairs = pairs2brackets(sol_rnastruct, rna)
calculate_sol_energy(sol_rnastruct, rna)

print("\nUNAFOLD\n")
generate_start_sol(dotbracket_una, sol_una)
fold, pairs = pairs2brackets(sol_una, rna)
calculate_sol_energy(sol_una, rna)
