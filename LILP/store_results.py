import pandas as pd
from LILP.utils.sol_converter import *
from LILP.utils.prepro_run import *
from LILP.utils.constants_paths import *
from optimize import *

def add_column(df, column_name, values):
    df[column_name] = values
    return df

ref_mfe_files = get_filenames(efn2_archive_dir, '.txt')
ref_mfe_names = [os.path.splitext(os.path.basename(f))[0] for f in ref_mfe_files]
results_df = pd.DataFrame(index=ref_mfe_names)

ref_MFEs = [get_energy_from_ct_file(f) for f in ref_mfe_files]

rna_mfe_names = get_filenames(ct_rnastruct_dir, '.ct')
rna_MFEs = [get_energy_from_ct_file(f) for f in rna_mfe_names]

vienna_mfe_names = get_filenames(dot_bracket_viennaRNA_dir, '_enrg.txt')
vienna_MFEs = [float(open(f).read().strip()) for f in vienna_mfe_names] 

unafold_mfe_names = get_filenames(unafold_fold_dir, '_enrg.txt')
unafold_MFEs = [float(open(f).read().strip()) for f in unafold_mfe_names]

add_column(results_df, 'MFE_ref', ref_MFEs)
add_column(results_df, 'MFE_rna', rna_MFEs)
add_column(results_df, 'MFE_vienna', vienna_MFEs)
add_column(results_df, 'MFE_unafold', unafold_MFEs)
print(results_df)

def completed_sequences(path):
    # Names of sequences that already have a full row in the results file.
    # A row is written only after a sequence is finished, so a sequence interrupted
    # by a server restart is missing here and will be computed again from scratch.
    if not os.path.exists(path):
        return set()
    with open(path) as fp:
        content = fp.read()
    if content and not content.endswith('\n'):
        # process killed while writing the last row: drop the partial row,
        # otherwise the next row would be glued to it
        content = content[:content.rfind('\n') + 1]
        with open(path, 'w') as fp:
            fp.write(content)
    rows = [line.split() for line in content.splitlines()[1:]]
    return {row[0] for row in rows if len(row) == RESULT_COLUMNS}


def sample_names(path, size):
    # benchmark/sample.tsv: every sequence of every stratum in a fixed random order
    # (column 'order'); the sample is the first `size` of each stratum, so a larger
    # size keeps all sequences of a smaller one
    sample = pd.read_csv(path, sep='\t')
    return set(sample[sample['order'] < size]['name'])


RESULT_COLUMNS = 22
model_name = os.environ.get('LILP_MODEL_NAME', 'cb-xyz') #'lilp-bcoaxial' #'lilp-branch'#'lilp-hairpin'
sample = (sample_names(os.environ['LILP_SAMPLE'], int(os.environ.get('LILP_SAMPLE_SIZE', 100)))
          if os.environ.get('LILP_SAMPLE') else None)
results_filename = f'LILP_{len_start}_{len_end}_{model_name}.txt'

# warm start (LILP/warmstart): LILP_START = none | lilp_init | file; with a start, the cost of
# producing it is taken from the budget LILP_TIME_LIMIT, statistics go to results/WS_<range>_<model>.tsv
from warmstart import run as warm_start
start_source = warm_start.source_name()
budget = int(os.environ.get('LILP_TIME_LIMIT', 7200))
ws_filename = f'WS_{len_start}_{len_end}_{model_name}.tsv'
done = completed_sequences(f'{results_dir}/{results_filename}')

n1 = int(os.environ.get('LILP_N1', 0))
n2 = int(os.environ.get('LILP_N2', len(seq_files)))
todo = [os.path.splitext(os.path.basename(f))[0] for f in seq_files[n1:n2]]
if sample is not None:
    todo = [name for name in todo if name in sample]
print(f'{model_name}: {len(todo)} sequences to compute, already done: {len(done & set(todo))}')

for seq_no in range (n1, n2):

    chain_file = seq_files[seq_no]
    chain_name_with_ext = os.path.basename(chain_file)
    chain_name_without_ext = os.path.splitext(chain_name_with_ext)[0]
    if chain_name_without_ext in done or chain_name_without_ext not in todo:
        continue
    lp_file_name = chain_name_without_ext
    seq_data = parse_seq_file(chain_file)
    rna = seq_data['sequence'].upper()

    first = 1
    last = len(rna)

    print(chain_name_without_ext)
    print(rna)
    print(len(rna))

    ###### START SOL OPTIMIZATION #######    
    start_name = 'lilp-bcoaxial-init'
    stem = True
    hairpin = True
    internal = True
    bulge = False
    branch = False 
    cbranch = False
    start = False
    
    # cut_MFE, lp_name, opt_time, gap, best_obj_time = optimize_lilp(rna, lp_file_name, start_name, stem, hairpin, internal, bulge, branch, cbranch, lpstart_dir, incumbent_dir, solstart_dir, first, last)

    # PPV_gen, STY_gen, MCC_gen, PPV_rnastruct, STY_rnastruct, MCC_rnastruct, PPV_vienna, STY_vienna, MCC_vienna, PPV_unafold, STY_unafold, MCC_unafold = sol_analyse(seq_files, seq_no, solstart_dir, start_name, dot_bracket_start_dir, dot_bracket_archive_dir, dot_bracket_rnastructure_dir, dot_bracket_viennaRNA_dir, unafold_fold_dir)

    ###### MAIN SOL OPTIMIZATION #######
    stem = True
    hairpin = True
    internal = True
    bulge = True
    branch = True
    cbranch = True    
    # multi = False
    start = False 

    #### NO START VERSION (or warm start from LILP_START) ####
    start_values, time_limit, ws, info = None, None, None, None
    if start_source != 'none':
        ws = warm_start.prepare(start_source, rna, chain_name_without_ext, model_name, budget, grb_log_dir)
        start_values, time_limit, info = ws.pop('start_values'), ws['time_limit'], {}
        print(f"WARM START {start_source}: {ws['start_pairs']} pairs, {ws['start_seconds']} s, time limit {time_limit:.1f} s")
    gen_MFE, lp_name, opt_time, gap, best_obj_time = optimize_lilp(rna, lp_file_name, model_name, stem, hairpin, internal, bulge, branch, cbranch, lp_dir, incumbent_dir, sol_dir, first, last, start_values=start_values, time_limit=time_limit, info=info)
    #### WITH START VERSION ####
    # gen_MFE, lp_name, opt_time, gap, best_obj_time = optimize_lilp(rna, lp_file_name, model_name, stem, hairpin, internal, bulge, branch, cbranch, lp_dir, incumbent_dir, sol_dir, first, last, start, start_name, solstart_dir)
   
    nan = float('nan')
    if best_obj_time is None:
        best_obj_time = nan
    if gen_MFE is None:
        # no solution within the time limit: still write a row (with NaN), so the
        # sequence counts as done and is not retried after every restart
        gen_MFE, gap = nan, nan
        PPV_gen, STY_gen, MCC_gen, PPV_rnastruct, STY_rnastruct, MCC_rnastruct, PPV_vienna, STY_vienna, MCC_vienna, PPV_unafold, STY_unafold, MCC_unafold = [nan] * 12
    else:
        PPV_gen, STY_gen, MCC_gen, PPV_rnastruct, STY_rnastruct, MCC_rnastruct, PPV_vienna, STY_vienna, MCC_vienna, PPV_unafold, STY_unafold, MCC_unafold = sol_analyse(seq_files, seq_no, sol_dir, model_name, dot_bracket_dir, dot_bracket_archive_dir, dot_bracket_rnastructure_dir, dot_bracket_viennaRNA_dir, unafold_fold_dir)

    if ws is not None:
        warm_start.write_row(f'{results_dir}/{ws_filename}', ws, info)
    write_results_to_file(lp_name, len(rna), opt_time, best_obj_time, gap, gen_MFE/100, ref_MFEs[seq_no], rna_MFEs[seq_no], vienna_MFEs[seq_no], unafold_MFEs[seq_no], round(PPV_gen,2), round(PPV_rnastruct,2), round(PPV_vienna,2), round(PPV_unafold,2), round(STY_gen,2), round(STY_rnastruct,2), round(STY_vienna,2), round(STY_unafold,2), round(MCC_gen,2), round(MCC_rnastruct,2), round(MCC_vienna,2), round(MCC_unafold,2), results_dir, results_filename)