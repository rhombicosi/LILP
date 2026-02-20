import re
import math
from stemloop import *
from hairpinloop import *
from internalloop import *
from bulgeloop import *
from multiloop import *
from utils.constants_paths import *
from utils.prepro_run import *
from utils.prepro_utils import *


def pairs2brackets(filepath, RNA): 
    lngth = len(RNA)
    pattern = r'_(\d+)_(\d+)' #"\((.*?)\)"
    fold = ["." for _ in range(lngth)] 
        
    with open(filepath) as fp:

        line = fp.readline()      
        cnt = 1
        pairs = []
        
        while line: 
            if " 1" in line and "P_" in line:
                # print("{}".format(line.strip()))
                match = re.search(pattern, line)
                i = int(match.group(1))
                j = int(match.group(2))

                pairs.append((i,j))
                fold[i - 1] = "("
                fold[j - 1] = ")"
                
            line = fp.readline()
            cnt += 1

    fold = ''.join([str(elem) for elem in fold])
    return (fold,pairs)

def brackets2pairs(dot_bracket):
    pair_stack = []
    base_pairs = []

    for i, symbol in enumerate(dot_bracket, start=1):
        if symbol in "({[":
            pair_stack.append(i)
        elif symbol in ")}]":
            if pair_stack:
                j = pair_stack.pop()
                base_pairs.append((j, i))

    sorted_bp = sorted(base_pairs, key=lambda x: x[0])

    return sorted_bp

def compare2folds(generated, reference):
    beta = 2
    TP = len(reference.intersection(generated)) # pairs that have been identified correctly
    FN = len(reference - generated) # pairs that have not been identified
    FP = len(generated - reference) # pairs that have been identified as part of the structure incorrectly

    if TP == 0 and (FP == 0 or FN == 0):
        f1 = 0
        fbeta = 0
        PPV = 0
        STY = 0
        MCC = 0
    else:    
        f1 = 2*TP/(2*TP + FP + FN)
        fbeta = (1+beta**2)*TP/((1+beta**2)*TP + (beta**2)*FN + FP)

        PPV = TP/(TP+FP)
        STY = TP/(TP+FN)
        MCC = math.sqrt(PPV*STY)

    return (f1,fbeta,MCC)

def calculate_sol_energy(filepath, rna):
    energy = 0
    with open(filepath, 'r') as file:
        lines = file.readlines()

    # Process lines
    for line in lines:
        if line.strip().endswith('1'):
            # Extract type strictly as a whole word before the first underscore
            type_match = re.match(r'^(STEM|HAIRPIN|INTERNAL|BULGE|MULTI)_', line)
            if type_match:
                element_type = type_match.group(1)
                
                # Extract all numbers after underscores
                numbers = re.findall(r'_(\d+)', line)
                indices = [int(num) for num in numbers]

                if element_type == 'STEM':
                    i1, j1, i2, j2 = indices
                    bp1 = BasePair(i1, j1, rna)
                    bp2 = BasePair(i2, j2, rna)
                    stemloop = StemLoop((bp1, bp2), rna)
                    print(f'{element_type} :: ({i1}, {j1}), ({i2}, {j2}) :: {stemloop.energy}')
                    energy += stemloop.energy

                if element_type == 'HAIRPIN':
                    i1, j1 = indices
                    bp1 = BasePair(i1, j1, rna)
                    hairpinloop = HairpinLoop([bp1], rna)
                    print(f'{element_type} :: ({i1}, {j1}) :: {hairpinloop.energy}')
                    energy += hairpinloop.energy

                if element_type == 'INTERNAL':
                    i1, j1, i2, j2 = indices
                    bp1 = BasePair(i1, j1, rna)
                    bp2 = BasePair(i2, j2, rna)
                    internalloop = InternalLoop((bp1, bp2), rna)
                    print(f'{element_type} :: ({i1}, {j1}), ({i2}, {j2}) :: {internalloop.energy}')
                    energy += internalloop.energy
                
                if element_type == 'BULGE':
                    i1, j1, i2, j2 = indices
                    bp1 = BasePair(i1, j1, rna)
                    bp2 = BasePair(i2, j2, rna)
                    bulgeloop = BulgeLoop((bp1, bp2), rna)
                    print(f'{element_type} :: ({i1}, {j1}), ({i2}, {j2}) :: {bulgeloop.energy}')
                    energy += bulgeloop.energy

                if element_type == 'MULTI':
                    i1, j1, i2, j2, i3, j3 = indices
                    bp1 = BasePair(i1, j1, rna)
                    bp2 = BasePair(i2, j2, rna)
                    bp3 = BasePair(i3, j3, rna)
                    multiloop = MultiLoop((bp1, bp2, bp3), rna)
                    print(f'{element_type} :: ({i1}, {j1}), ({i2}, {j2}) , ({i3}, {j3}) :: {multiloop.energy}')
                    energy += multiloop.energy
    print(energy)

def sol_analyse(seq_files, seq_number, sol_dir, model_name, dot_bracket_dir, dot_bracket_archive_dir, dot_bracket_rnastructure_dir, dot_bracket_viennaRNA_dir, dot_bracket_unafold_dir, start):

    chain_file = seq_files[seq_number]
    chain_name_with_ext = os.path.basename(chain_file)        
    chain_name_without_ext = os.path.splitext(chain_name_with_ext)[0]
    lp_file_name = chain_name_without_ext

    seq_data = parse_seq_file(chain_file)
    this_RNA = seq_data['sequence']    
    
    if start:
        print("START")
        filepath = os.path.join(solstart_dir, f'{lp_file_name}-start.sol')
    else:
        filepath = os.path.join(sol_dir, f'{lp_file_name}-{model_name}.sol')

    ic(filepath)

    # archive referece
    ref_bracket_path = f'{dot_bracket_archive_dir}/{lp_file_name}.txt'
    ref_brackets = dot_from_txt(ref_bracket_path)

    # rnastructure reference
    rnastruct_bracket_path = f'{dot_bracket_rnastructure_dir}/{lp_file_name}.txt'
    rnastruct_brackets = dot_from_txt(rnastruct_bracket_path)
    
    # viennaRNA reference
    viennaRNA_bracket_path = f'{dot_bracket_viennaRNA_dir}/{lp_file_name}_db.txt'
    with open(viennaRNA_bracket_path, 'r') as file:
        # Read all lines into a list
        lines = file.readlines()        
        veinnaRNA_brackets = str(lines[0]).strip()
        
    # UNAfold reference
    unafold_bracket_path = f'{dot_bracket_unafold_dir}/{lp_file_name}_db.txt'
    with open(unafold_bracket_path, 'r') as file:
        # Read all lines into a list
        lines = file.readlines()        
        unafold_brackets = str(lines[0]).strip()
    
    gen_brackets, gen_pairs = pairs2brackets(filepath, this_RNA)

    print(gen_brackets)
    print(gen_pairs)
    
    if start:
        file_bracket = f'{dot_bracket_dir}/{lp_file_name}-dotbrackets.txt'
    else: 
        file_bracket = f'{dot_bracket_dir}/{lp_file_name}-dotbrackets-{model_name}.txt'

    # with open(file_bracket, 'a') as file:
    with open(file_bracket, 'w') as file:
        file.write(gen_brackets)
        file.write("\n")

    # pairs generated by the model and pairs of a reference structure
    generated = set(gen_pairs)
    reference = set(brackets2pairs(ref_brackets))
    rnastruct = set(brackets2pairs(rnastruct_brackets))
    viennaRNA = set(brackets2pairs(veinnaRNA_brackets))
    unafold = set(brackets2pairs(unafold_brackets))

    print(generated)
    print(reference)
    print(rnastruct)
    print(viennaRNA)
    print(unafold)

    (f1_lilp,fbeta_lilp,MCC_lilp) = compare2folds(generated, reference)
    (f1_rnastruct,fbeta_rnastruct,MCC_rnastruct) = compare2folds(rnastruct, reference)
    (f1_viennaRNA,fbeta_viennaRNA,MCC_viennaRNA) = compare2folds(viennaRNA, reference)
    (f1_unafold,fbeta_unafold,MCC_unafold) = compare2folds(unafold, reference)

    print(f1_lilp)
    print(fbeta_lilp)
    print(MCC_lilp)

    print(f1_rnastruct)
    print(fbeta_rnastruct)
    print(MCC_rnastruct)
    
    print(f1_viennaRNA)
    print(fbeta_viennaRNA)
    print(MCC_viennaRNA)
    
    print(f1_unafold)
    print(fbeta_unafold)
    print(MCC_unafold)

    return f1_lilp, fbeta_lilp, MCC_lilp, f1_rnastruct, fbeta_rnastruct, MCC_rnastruct, f1_viennaRNA, fbeta_viennaRNA, MCC_viennaRNA, f1_unafold, fbeta_unafold, MCC_unafold
