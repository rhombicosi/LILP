from typing import List, Tuple, Optional
from enum import Enum
from basepair import *
from lilp_config import *

class Loop:

    def __init__(self, base_pairs: List[BasePair], RNA: str):
        self.RNA = RNA
        self.base_pairs = base_pairs
        self.degree = len(base_pairs)
        self.type = self._infer_loop_type()
        self.size = self._compute_loop_size()
        self.var = None

    def _infer_loop_type(self) -> LoopType:
        if self.degree == 1:
            return LoopType.HAIRPIN
        elif self.degree == 2:
            bp1 = self.base_pairs[0]
            bp2 = self.base_pairs[1]

            i1, j1 = bp1.i, bp1.j
            i2, j2 = bp2.i, bp2.j

            if i2 == i1 + 1 and j2 == j1 - 1:
                return LoopType.STEM
            elif (i2 == i1 + 1 and j2 < j1 - 1) or (j2 == j1 - 1 and i2 > i1 + 1):
                return LoopType.BULGE
            elif i2 > i1 + 1 and j2 < j1 - 1:
                return LoopType.INTERNAL
        elif self.degree > 2:
            return LoopType.MULTI
        else:
            raise ValueError("Invalid number of base pairs to determine loop type.")
        
    def _compute_loop_size(self) -> int:
        pairs = self.base_pairs
        d = self.degree

        if d != len(pairs):
            raise ValueError("Mismatch between declared degree and number of base pairs")

        if d == 1:
            bp = pairs[0]
            size = bp.j - bp.i - 1
        else:
            i1, j1 = pairs[0].i, pairs[0].j
            i2 = pairs[1].i
            jd = pairs[-1].j

            size = i2 - i1 + (j1 - jd) - d

            for k in range(1, d - 1):
                ik, jk = pairs[k + 1].i, pairs[k].j
                size += ik - jk

        return size    

    def is_valid_size(self) -> bool:
        pairs = self.base_pairs
        d = self.degree
        is_valid = True
        
        if d == 1:
            is_valid = self.size <= MAX_LOOP_SIZES[self.type]
        else:
            is_valid = pairs[1].i - pairs[0].i - 1 <= MAX_LOOP_SIZES[self.type] and pairs[0].j - pairs[d - 1].j - 1 <= MAX_LOOP_SIZES[self.type]
            for k in range(1, d - 1):
                ik, jk = pairs[k + 1].i, pairs[k].j
                is_d_valid = is_valid and ik - jk - 1 <= MAX_LOOP_SIZES[self.type]
                is_valid = is_d_valid

        return is_valid
    
    def add_variable(self, model: gp.Model):       
            
        indices = [str(idx) for bp in self.base_pairs for idx in (bp.i, bp.j)]
        index_str = "_".join(indices)
        
        name = f"{self.type.name}_{index_str}"
        self.var = model.addVar(vtype=GRB.BINARY, name=name)
        
        return self.var    
