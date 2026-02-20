import os
from pathlib import Path
import time
from lilp import *

seq_len = 60
seq_no = 0

cwd = Path.cwd()
code_path = Path(__file__).parent.parent
arch_rel_path = '../../ARCHIVE II/'
archive_path = (code_path/arch_rel_path).resolve()
seq_len_dir = f'RNA_seq_{seq_len}'

chain_dir = os.path.join(archive_path, seq_len_dir)
seq_files = get_filenames(chain_dir, '.seq')

chain_file = seq_files[seq_no]
chain_name_with_ext = os.path.basename(chain_file)        
chain_name_without_ext = os.path.splitext(chain_name_with_ext)[0]
lp_file_name = chain_name_without_ext
seq_data = parse_seq_file(chain_file)

rna = seq_data['sequence']
print(rna)

rna_model = LILP(rna)
rna_model.create_variables(1, 1, 1, 1, 1)
rna_model.create_constraints(1, 1, 1, 1, 1)
rna_model.create_objective(1, 1, 1, 1, 1)

# # Read the model
# model = gp.read(sys.argv[1])

# Set a time limit for the model
rna_model.model.Params.TimeLimit = 300

# Set a time limit for the whole tuning run
rna_model.model.Params.TuneTimeLimit = 1800

# Set the TuneResults parameter to 2
#
# The first parameter setting is the result for the first solved
# setting. The second entry the parameter setting of the best parameter
# setting.
rna_model.model.Params.TuneResults = 2

# Tune the model
rna_model.model.tune()

if rna_model.model.TuneResultCount >= 2:
    # Load the best tuned parameters into the model
    #
    # Note, the first parameter setting is associated to the first solved
    # setting and the second parameter setting to best tune result.
    rna_model.model.getTuneResult(1)

    # Write tuned parameters to a file
    rna_model.model.write("tune.prm")

    # Solve the model using the tuned parameters
    rna_model.model.optimize()