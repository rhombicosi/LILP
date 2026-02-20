import os
# from pathlib import Path
from lilp import *
# from utils.prepro_run import *
from utils.sol_converter import *

#len_start = 60
seq_number = 72

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

model_name = 'lilp-multi'
start_name = 'cut-start'
for f in range(18):
    filepath = f'{incumbent_dir}/{lp_file_name}-incumbent-{model_name}_{f}.sol'
    fold, pairs = pairs2brackets(filepath, rna)
    print(fold)
    # calculate_sol_energy(filepath, rna)


# filepath = f'{solstart_dir}/{lp_file_name}-{start_name}.sol'
# fold, pairs,lngth = pairs2brackets(filepath, rna)
# print(fold)
# calculate_sol_energy(filepath, rna)

# lp_file_name ='tRNA_tdbR00000592-Ascaris_suum-6253-Gln-UG'
# 'tRNA_tdbR00000009-Escherichia_coli-562-Ala-VGC'
# rna = 'GGGGGCAUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCUGCGGUUCGAUCCCGCGCGCUCCCACCA'


# filepath = f'{sol_dir}/{lp_file_name}-{model_name}.sol'
# fold, pairs,lngth = pairs2brackets(filepath, rna)
# print(fold)
# calculate_sol_energy(filepath, rna)

# filepath = f'{incumbent_dir}\{lp_file_name}-incumbent-{model_name}_3.sol'
# pairs2brackets(filepath, rna)
# calculate_sol_energy(filepath, rna)

# filepath = f'{dot_bracket_archive_dir}/{lp_file_name}.txt'
# print(filepath)
# print(dot_from_txt(filepath))
# print(brackets2pairs(dot_from_txt(filepath)))


# rna_model = LILP(rna, model_name)
# rna_model.create_variables(1, 1, 1, 1, 1)

# sorted_hairpins = sorted(rna_model.hairpin_loops, key=lambda x: (x.size, x.energy))
# for hl in sorted_hairpins:
#     if hl.is_valid_size():
#         print(f'{hl.var} :: {hl.size} :: {hl.energy}')

# sorted_stems = sorted(rna_model.stem_loops, key=lambda x: x.distance)   
# for sl in sorted_stems:
#     print(f'{hl.var} :: {hl.size} :: {hl.energy}')

#     sl.var.setAttr("BranchPriority", 10 * sl.distance)



# counter = 0
# for b in rna_model.bulge_loops:   
#     if b.var.VarName == 'BULGE_3_54_12_53':
#         print(b.size)
#     # print(f"[{counter}]:: Variable name: {b.var.VarName}, Energy: {b.energy}, Size: {b.size} Gurobi var: {b.var}")
#     counter+=1

# script_dir = os.path.dirname(os.path.abspath(__file__))
# parent_dir = os.path.dirname(script_dir)
# filepath = os.path.join(parent_dir, f'lilp-start-{seq_no}.sol')

# pairs2brackets(filepath, rna)
# calculate_sol_energy(filepath, rna)
# for f in range(15):
#     # filepath = os.path.join(parent_dir, f'lilp-start-{seq_no}-incumbent_{f}.sol')
#     filepath = os.path.join(parent_dir, f'lilp-start-{seq_no}.sol')
#     pairs2brackets(filepath, rna)