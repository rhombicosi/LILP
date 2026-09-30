from dloop import *
from tnn_parameters.initiation_parameters import *
from tnn_parameters.stem_parameters import *

class StemLoop(Loop):

    def __init__(self, base_pairs, RNA):
        super().__init__(base_pairs, RNA)  
        self.first_pair = base_pairs[0]  
        self.last_pair = base_pairs[1]
        self.energy = self.calculate_energy()
        self.distance = self.first_pair.j - self.last_pair.i - 1
            
    def _find_stem_neighbors(self, stem_loops: List[Loop], order: str) -> Loop:
        if order == "next":
            curr_loop_lp = self.last_pair
            return next((sl for sl in stem_loops if sl.first_pair.i == curr_loop_lp.i and sl.first_pair.j == curr_loop_lp.j), None)
        elif order == "previous":
            curr_loop_fp = self.first_pair
            return next((sl for sl in stem_loops if sl.last_pair.i == curr_loop_fp.i and sl.last_pair.j == curr_loop_fp.j), None)

    def calculate_energy(self) -> int:
        G = wcf_df.loc[self.first_pair.nt1 + self.first_pair.nt2, self.last_pair.nt1 + self.last_pair.nt2]
        return round(G)

    def create_stem_ifthen_constraint(self, model: gp.Model) -> None:
        inequality = gp.LinExpr([1, 1, -1], [self.first_pair.var, self.last_pair.var, self.var])
        model.addConstr(inequality <= 1, f'SLIT-{self.base_pairs[0].i}-{self.first_pair.j}')

    def create_stem_onlyif_constraint(self, model: gp.Model) -> None:
        inequality = gp.LinExpr([2, -1, -1], [self.var, self.first_pair.var, self.last_pair.var])
        model.addConstr(inequality <= 0, f'SLOI-{self.first_pair.i}-{self.first_pair.j}')
    
    def create_stem_constraints(self, model: gp.Model) -> None:
        # else:        
            inequality = gp.LinExpr([2, -1, -1], [self.var, self.first_pair.var, self.last_pair.var])
            model.addConstr(inequality <= 0, f'SLOI-{self.first_pair.i}-{self.first_pair.j}')
        # if self.energy > 0:
            inequality = gp.LinExpr([1, 1, -1], [self.first_pair.var, self.last_pair.var, self.var])
            model.addConstr(inequality <= 1, f'SLIT-{self.base_pairs[0].i}-{self.first_pair.j}')
        

    def create_first_pair_constraints(self, model: gp.Model, stem_loops: List["StemLoop"], first_pairs: List[BasePair]) -> None:
        
        sl_prev = self._find_stem_neighbors(stem_loops, "previous")   
        first_pair = BasePair._find_base_pairs_matches(first_pairs, self.first_pair.i, self.first_pair.j)    

        if sl_prev:
            inequality = gp.LinExpr([2, -1, 1], [first_pair.var, self.var, sl_prev.var])
            model.addConstr(inequality <= 1, f'FPIT-{first_pair.i}-{first_pair.j}')
            inequality = gp.LinExpr([1, -1, -1], [self.var, sl_prev.var, first_pair.var])
            model.addConstr(inequality <= 0, f'FPOI-{first_pair.i}-{first_pair.j}')
        else:
            inequality = gp.LinExpr([2, -1], [first_pair.var, self.var])
            model.addConstr(inequality <= 1, f'FPIT-{first_pair.i}-{first_pair.j}')
            inequality = gp.LinExpr([1, -1], [self.var, first_pair.var])
            model.addConstr(inequality <= 0, f'FPOI-{first_pair.i}-{first_pair.j}')

    def create_last_pair_constraints(self, model: gp.Model, stem_loops: List["StemLoop"], last_pairs: List[BasePair]) -> None:
                   
        sl_next = self._find_stem_neighbors(stem_loops, "next")
        last_pair = BasePair._find_base_pairs_matches(last_pairs, self.last_pair.i, self.last_pair.j) 

        if sl_next:
            inequality = gp.LinExpr([2, -1, 1], [last_pair.var, self.var, sl_next.var])
            model.addConstr(inequality <= 1, f'LPIT-{last_pair.i}-{last_pair.j}')
            inequality = gp.LinExpr([1, -1, -1], [self.var, sl_next.var, self.last_pair.var])
            model.addConstr(inequality <= 0, f'LPOI-{last_pair.i}-{last_pair.j}')                    
        else:
            inequality = gp.LinExpr([1, -1], [self.var, last_pair.var])
            model.addConstr(inequality <= 0, f'LPIT-{last_pair.i}-{last_pair.j}')
            inequality = gp.LinExpr([2, -1], [last_pair.var, self.var])
            model.addConstr(inequality <= 1, f'LPOI-{last_pair.i}-{last_pair.j}')   

# rna = "GCCGCGAACCCCGCCAGGCCCGGAAGGGAGCAACGGUAGUGGUGGAU"
# bp1 = BasePair(20,27,rna)
# bp2 = BasePair(21,26,rna)
# i_loop = StemLoop((bp1, bp2), rna)

# print(i_loop.energy)
# print(i_loop.size)
     
    