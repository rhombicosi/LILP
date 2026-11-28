"""Solves one sequence with LILP and prints the result as JSON (used by tests in a fresh process,
because the parameter set is chosen at import time from LILP_PARAMS).

    python test/solve_example.py GGAGCAGAAAGCAUCC
"""
import json
import sys
from collections import defaultdict

import gurobipy as gp
from lilp import LILP

seq = sys.argv[1]
flags = (True,) * 6
rna = LILP(seq, 'regression')
rna.create_variables(*flags, 1, len(seq))
rna.create_constraints(*flags, 1, len(seq))
rna.create_objective(*flags)
rna.model.setParam('OutputFlag', 0)
rna.model.setParam('Threads', 4)
rna.model.optimize()

objective = rna.model.getObjective()
energy_by_type = defaultdict(int)
for k in range(objective.size()):
    var = objective.getVar(k)
    if var.X > 0.5:
        energy_by_type[var.VarName.split('_')[0]] += round(objective.getCoeff(k))

structure = ['.'] * len(seq)
for var in rna.model.getVars():
    if var.VarName.startswith('P_') and var.X > 0.5:
        i, j = map(int, var.VarName.split('_')[1:3])
        structure[i - 1], structure[j - 1] = '(', ')'

print(json.dumps({'energy': round(rna.model.ObjVal), 'structure': ''.join(structure),
                  'energy_by_type': dict(sorted(energy_by_type.items()))}))
