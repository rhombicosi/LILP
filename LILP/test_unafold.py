import subprocess, tempfile, os
from pathlib import Path
from utils.prepro_utils import *
from utils.constants_paths import *

def run_unafold(seq: str):
    with tempfile.TemporaryDirectory() as tmp:
        seqfile = os.path.join(tmp, "seq.txt")
        with open(seqfile, "w") as f:
            f.write(seq + "\n")

        subprocess.run(["hybrid-ss-min", seqfile], cwd=tmp, check=True)
        ctfile = seqfile + ".ct"

        with open(ctfile) as f:
            return f.read()

# print(run_unafold("GCGCGAUUAGCUAGCGC"))

# ct = run_unafold("GCGCGAUUAGCUAGCGC")

def ct_to_dotbracket(ct_text: str) -> str:
    lines = [l.strip() for l in ct_text.splitlines() if l.strip()]
    pairs = {}
    for l in lines:
        if l[0].isdigit():
            cols = l.split()
            i = int(cols[0])
            j = int(cols[4])
            pairs[i] = j
    n = len(pairs)
    dot = []
    for i in range(1, n + 1):
        j = pairs[i]
        if j == 0:
            dot.append('.')
        elif j > i:
            dot.append('(')
        else:
            dot.append(')')
    return ''.join(dot)
    
# print(ct_to_dotbracket(ct))

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
        print("dG")
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


def run_unafold_from_file(seq_file_input: str, fold_dir, dot_bracket_dir):

    seq = seq_file_input
    seq_name_with_ext = os.path.basename(seq)
    seq_name_without_ext = os.path.splitext(seq_name_with_ext)[0]    
    seq_data = parse_seq_file(seq)
    sequence = seq_data['sequence']
    
    with tempfile.TemporaryDirectory() as tmp:
        seqfile = os.path.join(tmp, "seq.txt")
        with open(seqfile, "w") as f:
            f.write(sequence + "\n")
        
        subprocess.run(["hybrid-ss-min", seqfile], cwd=tmp, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        ctfile = seqfile + ".ct"
        
        # change working directory to path to RNAstructure exe files
        # os.chdir(rnastruct_path)
    
        result = subprocess.run(['./ct2dot', ctfile,'1', f'{dot_bracket_dir}/{seq_name_without_ext + ".txt"}'], capture_output=True, text=True)
        
        structure = os.path.join(fold_dir,f'{seq_name_without_ext + "_db.txt"}')
        energy = os.path.join(fold_dir,f'{seq_name_without_ext + "_enrg.txt"}')
        
        
        seq_str, dot_bracket, enrg = parse_ct(ctfile)
        
         # Write structure to a file
        with open(structure, "w") as f:
            f.write(dot_bracket + "\n")

        # Write energy to a file
        with open(energy, "w") as f:
            f.write(f"{enrg}\n")
            
        print(dot_bracket)
        print(enrg)
        
        with open(ctfile) as f:
            return f.read()

# Example usage

seq_file_input = "/path/to/srp_Acin.spec._CR543861.seq"
fold_dir = "/path/to/UNAfold_fold_0_50/"
dot_bracket_dir = "/path/to/UNAfold_fold_0_50/"


ct_content = run_unafold_from_file(seq_file_input, fold_dir, dot_bracket_dir)
print(ct_content)
