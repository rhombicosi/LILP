# Integer Programming Framework (LILP) for RNA 2D structure prediction.

Model performance is verified on Archive II dataset of sequences and structures, and compared to three DP methods RNAstructure, RNAfold and UNAfold.

# Installation guide for WSL VSCode

install Python 3.12 on WSL
install WSL extension for VScode
open working directory in WSL:UBUNTU in VSCode
install Gurobi with academic license

## ENVIRONMENT SETUP

`python3 -m venv path/to/venv/`

`source path/to/venv/bin/activate`

`pip install --upgrade pip`

`pip install -r requirements.txt`

Set correct python interpreter: Ctrl+Shift+P -> Python: Select Interpreter -> Enter interpreter path: path/to/venv/bin/python3.12

## PREREQUISITES

extract archiveii.tar to the code parent directory 

### RNAstructure 
install RNAstructure from tarball

add the path

`echo 'export PATH=$HOME/path/to/RNAstructure/exe:$PATH' >> ~/.bashrc`

### ViennaRNA (RNAfold)
sudo apt update
sudo apt install -y build-essential libgsl-dev python3-dev swig python3-pip pkg-config
wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_7_x/ViennaRNA-2.7.2.tar.gz
tar -xvzf ViennaRNA-2.7.2.tar.gz
cd ViennaRNA-2.7.2
./configure --prefix=$HOME/path/to/ViennaRNA
make
make check
make install

echo 'export PATH="$HOME/path/to/ViennaRNA/bin:$PATH"' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH="$HOME/path/to/ViennaRNA/lib:$LD_LIBRARY_PATH"' >> ~/.bashrc

### UNAfold

install unafold from tarball

export PATH="$HOME/path/to/unafold/bin:$PATH"

## CODE RUN
1. In constants_paths.py set the range of lengths of the sequences to test. for example:
    len_start = 0
    len_end = 50

2. Run prepro_run.py to initialize all folder names, paths and to prepare files for testing and references in correct format.

3. In store_results.py set the range of sequences to test. For example to run model on one first sequence in the list set:
    n1 = 0
    n2 = 1
    
    To run code on all sequences set:

    n1 = 0
    n2 = len(seq_files)

4. Run store_results.py to solve the LILP model and store the results in *ilp_results.txt*

## DIRECTORIES ORGANIZING 

**constant_paths.py**

names and paths to directories for storing reference structures, DP generated structures, LILP optimization solutions and logs:

*arch_rel_path* -- path to archiveii directory

*badfiles_dir* -- folder to store .seq files that do not have corresponding .ct structures, not useful for testing

*seq_len_dir*  -- directory to store sequences of size less or equal to {seq_len}

*ct_len_dir*  -- directory to store structures of sequences of size less or equal to {seq_len}

*lp_folder_name* -- folder to store models as .lp files

*sol_folder_name* --  folder to store solutions .sol files

*dot_bracket_folder_name* -- folder to store solutions in dot-bracket format

*dot_bracket_archive_folder_name* -- folder to store archive reference structures in dot-bracket format

*rnastructure_folder_name* -- folder to store structures .ct files ganerated by RNAstructure algorithm

*dot_bracket_rnastructure_folder_name* -- folder to store RNAstructure structures in dot-bracket format

*efn2_archive_folder_name* -- folder with .txt files with energies calculated with RNAstructure software for archive reference structures

*rnastruct_path* -- path to "RNAstructure6.5/exe" execution files

## DATA PREPROCESSING 

**prepro_run.py**

executes several commands to prepare reference structures to test loop-decomposition milp model 

**prepro_utils.py**

functions to preprocess .seq and .ct files, create folders and writing results to files
