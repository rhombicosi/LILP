from dloop import *
from tnn_parameters.initiation_parameters import *
from tnn_parameters.hairpin_parameters import *
from tnn_parameters.terminal_mismatch import *
from lilp_config import M

class HairpinLoop(Loop):

    def __init__(self, base_pairs, RNA):
        super().__init__(base_pairs, RNA)
        self.mistmatch_nt1 = RNA[base_pairs[0].i]
        self.mistmatch_nt2 = RNA[base_pairs[0].j - 2]
        self.energy = self.calculate_energy()

    def calculate_energy(self) -> int:
        if self.size == 3:
            G = initiation_df.loc[self.size, "hairpin"]
        elif self.is_valid_size():
            G = initiation_df.loc[self.size, "hairpin"] + mismatch_df.loc[self.base_pairs[0].nt1 + self.base_pairs[0].nt2, self.mistmatch_nt1][self.mistmatch_nt2]

            if self.base_pairs[0].nt1 + self.base_pairs[0].nt2 == 'GU' and self.RNA[self.base_pairs[0].i - 3] + self.RNA[self.base_pairs[0].j - 2] == 'GG':
                G += spec_GU_clos
            
            if self.mistmatch_nt1 + self.mistmatch_nt2 == 'GA' or self.mistmatch_nt1 + self.mistmatch_nt2 == 'UU':
                G += hp_mismatch['GA']

            if self.mistmatch_nt1 + self.mistmatch_nt2 == 'GG':
                G += hp_mismatch['GG']
        else:
            G = M
        return round(G)

    def create_hairpin_size_constraint(self, model: gp.Model) -> None:
        if not self.is_valid_size():
            inequality = gp.LinExpr([1], [self.var])
            model.addConstr(inequality == 0, f'HS-{self.base_pairs[0].i}-{self.base_pairs[0].j}')
            
    def create_hairpin_ifthen_constraint(self, model: gp.Model) -> None:
        inequality = gp.LinExpr(0)
        
        for u in range(self.base_pairs[0].i + 1, self.base_pairs[0].j):
            nucleotide = model.getVarByName(f'X_{u}')
            inequality.add(gp.LinExpr([1], [nucleotide]))
        
        inequality.add(gp.LinExpr([1, -1],[self.base_pairs[0].var, self.var]))            
        model.addConstr(inequality <= self.size, f'HIT-{self.base_pairs[0].i}-{self.base_pairs[0].j}')

    def create_hairpin_onlyif_constraint(self, model: gp.Model, base_pairs: List[BasePair]) -> None:        
        for u in range(self.base_pairs[0].i + 1, self.base_pairs[0].j):
            inequality = gp.LinExpr([2], [self.var])
            matches = BasePair._find_base_pairs_with_index(base_pairs, u)

            for bp in matches:
                inequality.add(gp.LinExpr([1], [bp.var]))
            
            inequality.add(gp.LinExpr([-1], [self.base_pairs[0].var]))
            model.addConstr(inequality <= 1, f'HOI-{self.base_pairs[0].i}-{self.base_pairs[0].j}-{u}')

    def create_hairpin_max_number_constraint(model: gp.Model, hairpin_loops: List["HairpinLoop"]) -> None:
        inequality = gp.LinExpr(0)

        for hl in hairpin_loops:
            inequality.add(gp.LinExpr([1], [hl.var]))
        model.addConstr(inequality <= MAX_NUM_OF_LOOPS[hl.type], f'HMN')
        model.update()    


# rna = "GCCGCGAACCCCGCCAGGCCCGGAAGGGAGCAACGGUAGUGGUGGAU"
# bp1 = BasePair(21,26,rna)
# h_loop = HairpinLoop([bp1], rna)

# print(h_loop.energy)
# print(h_loop.size)