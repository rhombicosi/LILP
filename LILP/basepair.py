from typing import List
import gurobipy as gp
from gurobipy import GRB
from lilp_config import *
from tnn_parameters.stem_parameters import *

class BasePair:

    def __init__(self, i: int, j: int, RNA: str):
        self.i = i
        self.j = j
        self.rna = RNA
        self.nt1 = RNA[i-1]
        self.nt2 = RNA[j-1]
        self.var = None
        self.pair_penalty_energy = self.calculate_pair_penalty_energy()
        self.distance = self.j - self.i - 1

    def distance(self) -> int:
        return self.j - self.i

    def is_valid(self) -> bool:
        return f'{self.nt1}{self.nt2}' in VALID_PAIRS and BasePair.distance(self) > MIN_D
    
    def calculate_pair_penalty_energy(self) -> int:
        if self.nt1 + self.nt2 == 'AU' or self.nt1 + self.nt2 == 'UA' or self.nt1 + self.nt2 == 'GU' or self.nt1 + self.nt2 == 'UG':
            G = wcf_AU_end_penalty
        else:
            G = 0
        return round(G)
            
    def add_variable(self, model: gp.Model, label: str) -> gp.Var:
        if self.is_valid():
            name = f'{label}_{self.i}_{self.j}'
            self.var = model.addVar(vtype=GRB.BINARY, name=name)
        return self.var
    
    def _find_base_pairs_with_index(base_pairs: List["BasePair"], index: int) -> List["BasePair"]:
        return [bp for bp in base_pairs if bp.i == index or bp.j == index]
    
    def _find_base_pair_with_indices(base_pairs: List["BasePair"], index1: int, index2: int) -> List["BasePair"]:
        return [bp for bp in base_pairs if bp.i == index1 and bp.j == index2]
    
    def _find_base_pairs_matches(base_pairs: List["BasePair"], i: int, j: int) -> "BasePair":
        return next((bp for bp in base_pairs if bp.i == i and bp.j == j), None)
    
    def create_single_pair_constraint(model: gp.Model, base_pairs: List["BasePair"], i: int) -> None:        
        inequality = gp.LinExpr(0)
        matches = BasePair._find_base_pairs_with_index(base_pairs, i)

        if matches:
            for bp in matches:
                inequality.add(gp.LinExpr([1.0], [bp.var]))
            model.addConstr(inequality <= 1, f'SP-{i}')

    def create_no_crossing_constraint(model: gp.Model, bp1: "BasePair", bp2: "BasePair") -> None:        
        if bp2.i > bp1.i and bp2.i < bp1.j and bp2.j > bp1.j:
            inequality = gp.LinExpr(0)
            inequality.add(gp.LinExpr([1.0, 1.0],[bp1.var, bp2. var]))
            model.addConstr(inequality <= 1, f'NC-{bp1.i}-{bp1.j}-{bp2.i}-{bp2.j}')
