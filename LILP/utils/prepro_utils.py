import os
from collections import defaultdict
from icecream import ic
import shutil
import subprocess
import re
import RNA
import tempfile
# uncomment to run code in lilp.py
# from utils.constants_paths import *
# uncomment to run prepro_run
from constants_paths import *

def sort_numeric_alpha(lst):
    def numeric_alpha_key(s):
        return [int(part) if part.isdigit() else part for part in re.findall(r'\d+|\D+', s)]
    return sorted(lst, key=numeric_alpha_key)

# gets all filenames of a given type in a given directory
def get_filenames(dir_path, f_type):
    f_list = []

    for file in os.listdir(dir_path):
        if file.endswith(f_type):
            f_list.append(os.path.join(dir_path, file))        
    return f_list

# parses a .ct file and returns a list of nucleotides with their pairings
def parse_ct_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    num_nucleotides = int(lines[0].split()[0])
    nucleotides = []
    
    for line in lines[1:num_nucleotides + 1]:
        parts = line.split()
        nucleotide_info = {
            'number': int(parts[0]),
            'nucleotide': parts[1],
            'prev': int(parts[2]),
            'next': int(parts[3]),
            'pair': int(parts[4])
        }
        nucleotides.append(nucleotide_info)
    
    return nucleotides

# parses a .seq file and returns a dictionary with the sequence identifier and sequence
def parse_seq_file(filename):
    sequence_data = {}

    with open(filename, 'r') as file:
        lines = file.readlines()
    
    sequence_data['identifier'] = lines[1]
    seq = ''.join(lines[2])
    sequence_data['sequence'] = seq[:-2]
    return sequence_data

# finds sequences that do not have 2D structure
def find_files_without_pairs(directory):
    files_dict = defaultdict(set)

    for filename in os.listdir(directory):
        if filename.endswith('.seq') or filename.endswith('.ct'):
            base_name, ext = os.path.splitext(filename)
            files_dict[base_name].add(ext)
    
    files_without_pairs = []

    for base_name, extensions in files_dict.items():
        if len(extensions) == 1:
            ext = list(extensions)[0]
            files_without_pairs.append(base_name + ext)
    
    return files_without_pairs

# moves files from one directory to the new directory
def create_directory_and_move_selected_files(src_dir, dest_dir, new_folder_name, selected_files):
    new_folder_path = os.path.join(dest_dir, new_folder_name)
    os.makedirs(new_folder_path, exist_ok=True)
    
    for filename in selected_files:
        src_file = os.path.join(src_dir, filename)
        
        if os.path.isfile(src_file):
            shutil.move(src_file, new_folder_path)
        else:
            ic(f"File not found: {filename}")

# copies files from one directory to the new directory
def create_directory_and_copy_selected_files(src_dir, dest_dir, new_folder_name, selected_files):
    new_folder_path = os.path.join(dest_dir, new_folder_name)
    os.makedirs(new_folder_path, exist_ok=True)
   
    for filename in selected_files:
        src_file = os.path.join(src_dir, filename)        
        
        if os.path.isfile(src_file):
            shutil.copy(src_file, new_folder_path)
        else:
            print(f"File not found: {filename}")

def ct_has_pairs(ct_path):
    """
    Check if a CT file has any base pairs.
    Returns True if at least one base pair exists, otherwise False.
    """
    with open(ct_path) as f:
        for line in f:
            parts = line.strip().split()

            # skip empty or non-data lines (header lines like 'through 511')
            if not parts or not parts[0].isdigit() or len(parts) < 5:
                continue

            try:
                j = int(parts[4])
                if j > 0:
                    return True  # found at least one pairing
            except ValueError:
                continue  # skip lines like 'through'

    return False  # no base pairs found

# selects sequences of a specific length less then some value
def get_seq_of_len(seq_list, seq_len, ct_list):
    seq_len_list = []
    ct_len_list = []
    seq_len_files = []
    ct_len_files = []
    
    seq_list = sort_numeric_alpha(seq_list)
    ct_list = sort_numeric_alpha(ct_list)

    for (seq,ct) in zip(seq_list, ct_list):
        seq_data = parse_seq_file(seq)
        ct_data = parse_ct_file(ct)        

        if len(seq_data['sequence']) <= seq_len:
            seq_len_files.append(seq)
            ct_len_files.append(ct)
            seq_len_list.append(seq_data)
            ct_len_list.append(ct_data)
    return seq_len_files,ct_len_files

# selects sequences of length between x nts and y nts
def get_seq_btwn_len(seq_list, len_start, len_end, ct_list):
    seq_len_list = []
    ct_len_list = []
    seq_len_files = []
    ct_len_files = []
    
    seq_list = sort_numeric_alpha(seq_list)
    ct_list = sort_numeric_alpha(ct_list)

    for (seq,ct) in zip(seq_list, ct_list):
        seq_data = parse_seq_file(seq)
        ct_data = parse_ct_file(ct)

        pairs = ct_has_pairs(ct)

        if not pairs:        
            print(f'NO PAIRS :: {ct}')
        
        if pairs and len(seq_data['sequence']) >= len_start and len(seq_data['sequence']) <= len_end:
            seq_len_files.append(seq)
            ct_len_files.append(ct)
            seq_len_list.append(seq_data)
            ct_len_list.append(ct_data)
    return seq_len_files,ct_len_files

# new_dir_name = 'lp'; rel_path_to_save = '../'
def create_dir(path_to_save, new_dir_name):
    # file_parent_dir = Path(__file__).parent
    # path_to_save = (file_parent_dir/rel_path_to_save).resolve()
    new_dir_path = os.path.join(path_to_save, new_dir_name)
    os.makedirs(new_dir_path, exist_ok=True)
    return new_dir_path

# converts .ct files to dot-bracket notation
def ct2dot(ct_files, first_file, last_file, dot_bracket_dir):
    for ct_number in range(first_file, last_file):
        ct = ct_files[ct_number]
        ct_name_with_ext = os.path.basename(ct)
        ct_name_without_ext = os.path.splitext(ct_name_with_ext)[0]

        result = subprocess.run(['ct2dot', ct, '1', f'{dot_bracket_dir}/{ct_name_without_ext + ".txt"}'], capture_output=True, text=True)
        print(result.stdout)
        


# generate .ct with RNA structure
def rnastruct_fold(seq_files, first_file, last_file, fold_dir):
    for seq_number in range(first_file, last_file):
        seq = seq_files[seq_number]
        seq_name_with_ext = os.path.basename(seq)
        seq_name_without_ext = os.path.splitext(seq_name_with_ext)[0]

        ct = os.path.join(fold_dir,f'{seq_name_without_ext + ".ct"}')

        result = subprocess.run(['Fold', seq, ct], capture_output=True, text=True)
        print(result.stdout)

# calculate energy based on .ct file data
def rnastruct_efn2(ct_files, first_file, last_file, efn2_dir):
    for ct_number in range(first_file, last_file):
        ct = ct_files[ct_number]
        ct_name_with_ext = os.path.basename(ct)
        ct_name_without_ext = os.path.splitext(ct_name_with_ext)[0]

        efn2 = os.path.join(efn2_dir,f'{ct_name_without_ext + ".txt"}')

        result = subprocess.run(['efn2', ct, efn2], capture_output=True, text=True)
        print(result.stdout)

# generate structs and enrgies with viennaRNA 
def viennaRNA_fold(seq_files, first_file, last_file, fold_dir):
    for seq_number in range(first_file, last_file):
        seq = seq_files[seq_number]
        seq_name_with_ext = os.path.basename(seq)
        seq_name_without_ext = os.path.splitext(seq_name_with_ext)[0]
        seq_data = parse_seq_file(seq)
        sequence = seq_data['sequence']

        structure = os.path.join(fold_dir,f'{seq_name_without_ext + "_db.txt"}')
        energy = os.path.join(fold_dir,f'{seq_name_without_ext + "_enrg.txt"}')
        
        #with open(seq_parsed, "r") as f:
            #sequence = f.read().strip()  # Remove any newline or whitespace 
            #print(sequence)       
        
        struct, enrg = RNA.fold(sequence)

        # Write structure to a file
        with open(structure, "w") as f:
            f.write(struct + "\n")

        # Write energy to a file
        with open(energy, "w") as f:
            f.write(f"{enrg}\n")
            
def parse_ct(ct_path: str):
    """Parse UNAfold .ct file into dot-bracket notation and extract free energy."""
    sequence = []
    pairs = {}
    energy = None

    with open(ct_path) as f:
        lines = [line.strip() for line in f if line.strip()]

    # Header line: contains energy
    header = lines[0]
    if "dG" in header:
        try:
            energy = float(header.split("dG =")[1].split()[0])
            print(energy)
        except Exception:
            energy = None

    # Parse the base pairing information
    for line in lines[1:]:
        cols = line.split()
        i = int(cols[0])
        base = cols[1]
        pair = int(cols[4])
        sequence.append(base)
        if pair != 0:
            pairs[i] = pair

    n = len(sequence)
    dot = ["." for _ in range(n)]

    for i, j in pairs.items():
        if i < j:
            dot[i - 1] = "("
            dot[j - 1] = ")"

    dot_bracket = "".join(dot)
    seq_str = "".join(sequence)
    return seq_str, dot_bracket, energy

# generate structs UNAFold
def UNAfold_fold(seq_files, first_file, last_file, fold_dir):
    for seq_number in range(first_file, last_file):
        seq_file = seq_files[seq_number]
        seq_name_with_ext = os.path.basename(seq_file)
        seq_name_without_ext = os.path.splitext(seq_name_with_ext)[0]    
        seq_data = parse_seq_file(seq_file)
        sequence = seq_data['sequence']
        
        with tempfile.TemporaryDirectory() as tmp:
            seqfile = os.path.join(tmp, "seq.txt")
            with open(seqfile, "w") as f:
                f.write(sequence + "\n")
            
            subprocess.run(["hybrid-ss-min", seqfile], cwd=tmp, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            
            ctfile = seqfile + ".ct"
            
            structure = os.path.join(fold_dir,f'{seq_name_without_ext + "_db.txt"}')
            energy = os.path.join(fold_dir,f'{seq_name_without_ext + "_enrg.txt"}')            
            
            seq_str, dot_bracket, enrg = parse_ct(ctfile)
            
             # Write structure to a file
            with open(structure, "w") as f:
                f.write(dot_bracket + "\n")
    
            # Write energy to a file
            with open(energy, "w") as f:
                f.write(f"{enrg}\n")

# retrieves dotbrackets from txt file
def dot_from_txt(f_txt):
    with open(f_txt, 'r') as file:
        lines = file.readlines()
        
        if len(lines) >= 3:
            third_line = str(lines[2]).strip()
            return(third_line)
        else:
            pass

# retrieves energy from ct file
def get_energy_from_ct_file(file_path):
    with open(file_path, 'r') as file:
        first_line = file.readline().strip()
        
        if "ENERGY" or "Energy" in first_line:
            parts = first_line.split()
            for i, part in enumerate(parts):
                if part == "ENERGY" or part == "Energy":
                    energy_value = float(parts[i + 2])
                    return energy_value
        else:
            print("ENERGY value not found in the first line.")
            return None

# reads variable values of start solution     
def read_sol(f_name):
    solvars = {}
    with open(f_name,'r') as f:
            lines = f.readlines()
            for line in lines[2:]:
                parts = line.strip().split()
                if len(parts) >= 2:  # Ensure at least two parts exist
                    key, value = parts[0], parts[1]
                    solvars[key] = value
    return solvars

# store optimization results
def write_results_to_file(sequence_name, rna_len, time, mfe_gen, mfe_ref, mfe_rna, mfe_vienna, mfe_unafold, f1_gen, f1_rna, f1_vienna, f1_unafold, fb_gen, fb_rna, fb_vienna, fb_unafold, mcc_gen, mcc_rna, mcc_vienna, mcc_unafold, results_dir, results_filename):
    
    headers = ["RNAseqname", "numofnts", "Time(s)", "MFEILP", "MFEARCHIVE", "MFERNAstr", "MFERNAFold", "MFEUNAfold", "F1ILP", "F1RNAstr", "F1RNAFold", "F1UNAFold", "FbILP", "FbRNAstr", "FbRNAFold", "FbUNAFold", "INFILP", "INFRNAstr", "INFRNAFold", "INFUNAFold"]
    
    filename = f'{results_dir}/{results_filename}'
    print(filename)
    
    # Check if the file exists and is not empty
    file_exists = os.path.exists(filename) and os.path.getsize(filename) > 0
    
    # Format the floating point numbers to 2 decimal places and ensure all values are strings
    values = [
        sequence_name, 
        rna_len,
        f"{time:.2f}",
        f"{mfe_gen:.2f}", 
        f"{mfe_ref:.2f}", 
        f"{mfe_rna:.2f}", 
        f"{mfe_vienna:.2f}",
        f"{mfe_unafold:.2f}",
        f"{f1_gen:.2f}", 
        f"{f1_rna:.2f}",
        f"{f1_vienna:.2f}",
        f"{f1_unafold:.2f}",
        f"{fb_gen:.2f}", 
        f"{fb_rna:.2f}",
        f"{fb_vienna:.2f}",
        f"{fb_unafold:.2f}",
        f"{mcc_gen:.2f}",
        f"{mcc_rna:.2f}",
        f"{mcc_vienna:.2f}",
        f"{mcc_unafold:.2f}"
    ]

    headers_line = "{:<45}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\n".format(*headers)
    
    line = "{:<45}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\n".format(*values)
    
    with open(filename, 'a') as file:
        if not file_exists:
            file.write(headers_line)
        
        file.write(line)
