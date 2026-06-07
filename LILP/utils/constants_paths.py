from pathlib import Path

# archive ii path
# with .ct and .seq files
cwd = Path.cwd()
code_path = Path(__file__).parent.parent
arch_rel_path = '../../archiveII/'
archive_path = (code_path/arch_rel_path).resolve()

# folder name for files that should be excluded from test
badfiles_dir = 'noctfiles'  

# length for sequences to be tested
len_start = 0
len_end = 50

# folders to save .seq and .ct files with sequences of seq_len
seq_len_dir = f'RNA_seq_{len_start}_{len_end}'
ct_len_dir = f'RNA_ct_{len_start}_{len_end}'

# folders to save lp models and solutions
ilp_parent_dir = Path(__file__).parent
lp_rel_path = '../../'
lp_path = (ilp_parent_dir/lp_rel_path).resolve()
lp_folder_name = 'lp_branch'
sol_folder_name = 'sol_branch'
lpstart_folder_name = 'lp_start'
solstart_folder_name = 'sol_start'
incumbent_folder_name = 'incumbent'
incumbent_start_folder_name = 'incumbent_start'
dot_bracket_folder_name = f'dot_bracket_{len_start}_{len_end}'
dot_bracket_start_folder_name = f'dot_bracket_start_{len_start}_{len_end}'
dot_bracket_archive_folder_name = f'dot_bracket_archive_{len_start}_{len_end}'
dot_bracket_viennaRNA_folder_name = f'dot_bracket_viennaRNA_{len_start}_{len_end}'
rnastructure_folder_name = f'RNAstructure_fold_{len_start}_{len_end}'
unafold_folder_name = f'UNAfold_fold_{len_start}_{len_end}'
dot_bracket_rnastructure_folder_name = f'dot_bracket_RNAstructure_{len_start}_{len_end}'
efn2_archive_folder_name = f'efn2_archive_{len_start}_{len_end}'
grb_log = f'gurobi_log'
results_folder_name = f'results'
results_file = f'LILP_{len_start}_{len_end}_results.txt'
