import os
# uncomment to run code in lilp.py
# from utils.prepro_utils import *
# from utils.constants_paths import *
# uncomment to run prepro_run
from prepro_utils import *
from constants_paths import *

# create folders to write .lp model and optimization .sol results
lp_dir = create_dir(lp_path, lp_folder_name)
sol_dir = create_dir(lp_path, sol_folder_name)
lpstart_dir = create_dir(lp_path, lpstart_folder_name)
solstart_dir = create_dir(lp_path, solstart_folder_name)
incumbent_dir = create_dir(lp_path, incumbent_folder_name)
incumbent_start_dir = create_dir(lp_path, incumbent_start_folder_name)
dot_bracket_dir = create_dir(lp_path, dot_bracket_folder_name)
dot_bracket_start_dir = create_dir(lp_path, dot_bracket_start_folder_name)
dot_bracket_archive_dir = create_dir(archive_path, dot_bracket_archive_folder_name)
dot_bracket_viennaRNA_dir = create_dir(archive_path, dot_bracket_viennaRNA_folder_name)
rnastructure_fold_dir =  create_dir(archive_path, rnastructure_folder_name)
unafold_fold_dir =  create_dir(archive_path, unafold_folder_name)
dot_bracket_rnastructure_dir = create_dir(archive_path, dot_bracket_rnastructure_folder_name)
efn2_archive_dir = create_dir(archive_path, efn2_archive_folder_name)
grb_log_dir = create_dir(lp_path, grb_log)
results_dir = create_dir(lp_path, results_folder_name)

# get .seq and .ct files after clearing redundant files
seq_list = get_filenames(archive_path, ".seq")
ct_list = get_filenames(archive_path, ".ct")

# create path to archive .ct files 
ct_archive_dir = os.path.join(archive_path, ct_len_dir)

# create path to rnastructure .ct files 
ct_rnastruct_dir = os.path.join(archive_path, rnastructure_fold_dir)

# create path to .seq files 
seq_dir = os.path.join(archive_path, seq_len_dir)

if __name__ == "__main__":

    # clean up ARCHIVE dir by getting rid off .ct/.seq files that do not .seq/.ct counterpart

    # find files that do not corresponding .ct or .seq paired file
    directory_path = archive_path 
    unpaired_files = find_files_without_pairs(directory_path)

    # create 'noctfiles' folder and move .seq files without .ct to it 
    source_directory = archive_path  
    destination_directory = archive_path 
    selected_files = unpaired_files 

    create_directory_and_move_selected_files(source_directory, destination_directory, badfiles_dir, selected_files)

    # seq_len_files,ct_len_files = get_seq_of_len(seq_list,seq_len,ct_list)
    seq_len_files,ct_len_files = get_seq_btwn_len(seq_list, len_start, len_end, ct_list)

    create_directory_and_copy_selected_files(source_directory, destination_directory, seq_len_dir, seq_len_files)
    create_directory_and_copy_selected_files(source_directory, destination_directory, ct_len_dir, ct_len_files)    

    # prepare dot bracket notations and 
    # get energies for ARCHIVE references;
    # generate reference structures with RNAstructure and
    # prepare dot bracket notations and 
    # get energies for RNAstructure references

    # change working directory to path to RNAstructure exe files
    # os.chdir(rnastruct_path)

    # get all paths to archive .ct files 
    ct_archive_files = get_filenames(ct_archive_dir, '.ct')

    # obtain dot-bracket notations for all archive .ct references
    ct2dot(ct_archive_files, 0, len(ct_archive_files), dot_bracket_archive_dir)

    # obtain enrgies for all archive .ct files 
    rnastruct_efn2(ct_archive_files, 0, len(ct_archive_files), efn2_archive_dir)

    # get all paths to .seq files    
    seq_files = get_filenames(seq_dir, '.seq')

    # generate secondary structures with RNAstructure algorithm, base settings
    rnastruct_fold(seq_files, 0, len(seq_files), rnastructure_fold_dir)   

    # get all paths to rnastructure .ct files
    ct_rnastruct_files = get_filenames(ct_rnastruct_dir, '.ct') 

    # obtain dot-bracket notations for all rnastructure .ct references
    ct2dot(ct_rnastruct_files, 0, len(ct_rnastruct_files), dot_bracket_rnastructure_dir)
    
    # obtain structures and enrgies with viennaRNA
    viennaRNA_fold(seq_files, 0, len(seq_files), dot_bracket_viennaRNA_dir)
    
    # change working directory to path to unafold files
    # os.chdir(unafold_path)
    
    # unafold
    UNAfold_fold(seq_files, 0, len(seq_files), unafold_fold_dir)

# path to sequences from archive ii that are currently under test
# which should be available for further calculations
chain_dir = os.path.join(archive_path, seq_len_dir)
seq_files = get_filenames(chain_dir, '.seq')

# print(f' # of sequences of lengths <= {seq_len} :: {len(seq_files)}')
print(f' # of sequences between {len_start} and {len_end} :: {len(seq_files)}')