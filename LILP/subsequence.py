from typing import List
from basepair import *
from lilp_config import *

class Subsequence:

    def __init__(self, start: int, end: int):
        self.start = start
        self.end = end
        self.var = None

    def add_variable(self, model: gp.Model):
        name = f'SSQ_{self.start}_{self.end}'
        self.var = model.addVar(vtype=GRB.BINARY, name=name)
        return self.var

    def create_unpaired_subsequence_constaint(self, model: gp.Model, base_pairs: List["BasePair"]):
        inequality = gp.LinExpr(0)
        matches = BasePair._find_base_pair_in_ssq(base_pairs, self.start, self.end)

        if matches:
            for bp in matches:
                inequality.add(gp.LinExpr([1.0], [bp.var]))
            
        inequality.add(gp.LinExpr([1.0], [self.var]))
        model.addConstr(inequality >= 1, f'USSQ_{self.start}_{self.end}')
        