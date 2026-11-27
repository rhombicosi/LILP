"""Sources of a start structure. Each returns (dot-bracket structure, cost in seconds, extra info)."""
import csv
import math
import re
import time

from lilp import LILP
from warmstart.structure import to_dot_bracket

# the simplified model = 2-degree LIP: stem, hairpin, internal and bulge loops, no multiloop
# variables. A loop type without variables is not forbidden but free (energy not counted,
# size limits not checked): 'lilp-bcoaxial-init' in store_results.py also leaves out the bulges,
# which made 4 of 8 test starts positive in the full model (1 of 8 with bulges); multiloops
# stay free here, so a start with a multiloop the full model does not allow is rejected.
INIT_FLAGS = dict(stem=True, hairpin=True, internal=True, bulge=True, branch=False, cbranch=False)
PAIR_VAR = re.compile(r'^P_(\d+)_(\d+)$')

_files = {}


def from_file(path, name):
    """Structure predicted beforehand by any program, e.g. benchmark/mxfold2.tsv.

    Columns: name, structure, seconds (cost of the prediction, taken from the budget).
    """
    if path not in _files:
        with open(path, newline='') as fp:
            _files[path] = {row['name']: row for row in csv.DictReader(fp, delimiter='\t')}
    row = _files[path].get(name)
    if row is None:
        raise KeyError(f'{name} is not in {path}')
    return row['structure'], float(row['seconds']), {}


def lilp_init(rna, name, time_limit, threads, mip_gap, log_file):
    """Build and solve the simplified LILP model; its structure is the start of the full model."""
    start = time.perf_counter()
    m = LILP(rna, f'{name}-init')
    m.create_variables(*INIT_FLAGS.values(), 1, len(rna))
    m.create_constraints(*INIT_FLAGS.values(), 1, len(rna))
    m.create_objective(*INIT_FLAGS.values())
    build = time.perf_counter() - start
    model = m.model
    model.setParam('LogToConsole', 0)
    model.setParam('LogFile', log_file)
    model.setParam('TimeLimit', time_limit)
    model.setParam('MIPGap', mip_gap)
    model.setParam('Threads', threads)
    model.optimize()
    pairs = []
    if model.SolCount > 0:
        for v in model.getVars():
            match = PAIR_VAR.match(v.VarName)
            if match and v.X > 0.5:
                pairs.append((int(match.group(1)), int(match.group(2))))
    info = {
        'init_build': round(build, 2),
        'init_runtime': round(model.Runtime, 2),
        'init_status': model.Status,
        'init_obj': model.ObjVal if model.SolCount > 0 else math.nan,
        'init_gap': model.MIPGap if model.SolCount > 0 else math.nan,
    }
    model.dispose()
    return to_dot_bracket(sorted(pairs), len(rna)), time.perf_counter() - start, info
