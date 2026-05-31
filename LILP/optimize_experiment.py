import os
# from pathlib import Path
from lilp import *
# from utils.prepro_run import *
from utils.sol_converter import *

#len_start = 60
seq_number = 89

# cwd = Path.cwd()
# code_path = Path(__file__).parent.parent
# arch_rel_path = '../../ARCHIVE II/'
# archive_path = (code_path/arch_rel_path).resolve()
# seq_len_dir = f'RNA_seq_{seq_len}'

# chain_dir = os.path.join(archive_path, seq_len_dir)
# seq_files = get_filenames(chain_dir, '.seq')

chain_file = seq_files[seq_number]
chain_name_with_ext = os.path.basename(chain_file)        
chain_name_without_ext = os.path.splitext(chain_name_with_ext)[0]
lp_file_name = chain_name_without_ext
seq_data = parse_seq_file(chain_file)

rna = seq_data['sequence'].upper()
print(rna)

# bp1 = BasePair(3,54,rna)
# bp2 = BasePair(12,53,rna)

model_name = 'lilp-cbranch'
start_name = 'lilp-branch-init'
for f in range(8):
    filepath = f'{incumbent_dir}/{lp_file_name}-incumbent-{model_name}_{f}.sol'
    fold, pairs = pairs2brackets(filepath, rna)
    print(fold)
    # calculate_sol_energy(filepath, rna)


filepath = f'{sol_dir}/{lp_file_name}-{model_name}.sol'
fold, pairs = pairs2brackets(filepath, rna)
print(fold)
calculate_sol_energy(filepath, rna)
