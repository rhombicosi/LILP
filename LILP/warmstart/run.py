"""Glue between store_results.py and the warm start: prepare the start of one sequence,
and write one row of warm-start statistics (results/WS_<range>_<model>.tsv)."""
import math
import os
import re

from warmstart import sources
from warmstart.structure import parse, repair, start_values, to_dot_bracket

MIP_GAP = 0.002                      # the same as the full model in optimize.py
SOURCES = ('none', 'lilp_init', 'file')
OBJECTIVE = re.compile(r'objective\s+([-+0-9.eE]+)')

COLUMNS = ['name', 'length', 'source', 'start_seconds', 'time_limit', 'pairs_in', 'pseudoknots',
           'removed_noncanonical', 'removed_short_loop', 'start_pairs', 'start_vars_set',
           'start_vars_missing', 'start_obj_first', 'start_obj', 'start_log', 'first_sol_time', 'first_sol_obj',
           'n_sols', 'status', 'runtime', 'final_obj', 'init_build', 'init_runtime', 'init_status',
           'init_obj', 'init_gap', 'structure_in', 'structure_start']


def source_name():
    source = os.environ.get('LILP_START', 'none')
    if source not in SOURCES:
        raise ValueError(f'LILP_START={source!r}, expected one of {SOURCES}')
    return source


def prepare(source, rna, name, model_name, budget, log_dir):
    """Structure from the source -> repaired partial MIP start + time limit of the full model."""
    threads = int(os.environ.get('LILP_THREADS', 8))
    if source == 'file':
        structure, seconds, extra = sources.from_file(os.environ['LILP_START_FILE'], name)
    elif source == 'lilp_init':
        structure, seconds, extra = sources.lilp_init(
            rna, name, float(os.environ['LILP_T_INIT']), threads, MIP_GAP,
            f'{log_dir}/{name}-log-{model_name}-init')
    else:
        raise ValueError(source)
    if len(structure) != len(rna):
        raise ValueError(f'{name}: structure of length {len(structure)} for a sequence of length {len(rna)}')
    pairs, pseudoknots = parse(structure)
    kept, removed = repair(pairs, rna)
    return {
        'name': name, 'length': len(rna), 'source': source,
        'start_seconds': round(seconds, 3),
        'time_limit': max(budget - seconds, 1.0),
        'pairs_in': len(pairs), 'pseudoknots': pseudoknots,
        'removed_noncanonical': removed['noncanonical'], 'removed_short_loop': removed['short_loop'],
        'start_pairs': len(kept),
        'structure_in': structure, 'structure_start': to_dot_bracket(kept, len(rna)),
        'start_values': start_values(kept) if kept else None,
        **extra,
    }


def start_objective(messages, first=False):
    """Objective of the MIP start reported by Gurobi, NaN if the start was not used.

    Gurobi completes the partial start ('User MIP start produced solution with objective ...',
    first=True: the first completion) and may improve it while processing the start; the
    search then begins from 'Loaded user MIP start with objective ...' (first=False).
    """
    produced, loaded = [], []
    for line in messages:
        low = line.lower()
        match = OBJECTIVE.search(line)
        if not match or 'did not' in low:
            continue
        if 'produced solution' in low:
            produced.append(float(match.group(1)))
        elif 'loaded' in low:
            loaded.append(float(match.group(1)))
    if first:
        return produced[0] if produced else (loaded[0] if loaded else math.nan)
    return loaded[-1] if loaded else (produced[-1] if produced else math.nan)


def write_row(path, ws, info):
    """One row per sequence; a sequence computed again after a restart gets a second row
    (the analysis keeps the last one)."""
    messages = info.get('start_messages', [])
    row = dict(ws)
    row.update({
        'start_vars_set': info.get('start_vars_set', 0),
        'start_vars_missing': info.get('start_vars_missing', 0),
        'start_obj_first': start_objective(messages, first=True),
        'start_obj': start_objective(messages),
        'start_log': ' | '.join(m.strip() for m in messages).replace('\t', ' '),
        'first_sol_time': info.get('first_sol_time', math.nan),
        'first_sol_obj': info.get('first_sol_obj', math.nan),
        'n_sols': info.get('n_sols', 0),
        'status': info.get('status', ''),
        'runtime': info.get('runtime', math.nan),
        'final_obj': info.get('final_obj', math.nan),
    })
    new = not os.path.exists(path)
    with open(path, 'a') as fp:
        if new:
            fp.write('\t'.join(COLUMNS) + '\n')
        fp.write('\t'.join(str(row.get(c, '')) for c in COLUMNS) + '\n')
