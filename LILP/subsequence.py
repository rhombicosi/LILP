from typing import List
from basepair import *
from lilp_config import *

class Nt2Subsequence:

    def __init__(self, nt: int, start: int, end: int):
        self.nt = nt
        self.start = start
        self.end = end
        self.var = None

    def add_variable(self, model: gp.Model):

        name = f'NTSSQ_{self.nt}_{self.start}_{self.end}'
        self.var = model.addVar(vtype=GRB.BINARY, name=name)

    # def _find_base_pair_in_ssq(base_pairs: List["BasePair"], nt: int, start: int, end: int) -> List["BasePair"]:
    #     return [bp for bp in base_pairs if (bp.i > start and bp.i < end and bp.j == nt) or (bp.i == nt and  bp.j > start and bp.j < end)]

    def create_unpaired_nt_constaint(self, model: gp.Model, base_pairs: List["BasePair"]):
        inequality = gp.LinExpr(0)
        matches = BasePair._find_nt_base_pair_in_ssq(base_pairs, self.nt, self.start, self.end)

        if matches:
            for bp in matches:
                inequality.add(gp.LinExpr([1.0], [bp.var]))
            
        inequality.add(gp.LinExpr([1.0], [self.var]))
        model.addConstr(inequality >= 1, f'UNTSSQ_{self.nt}_{self.start}_{self.end}')
        