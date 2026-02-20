from dloop import *
from tnn_parameters.initiation_parameters import *
from tnn_parameters.internal_parameters import *
from lilp_config import *

class InternalLoop(Loop):

    def __init__(self, base_pairs, RNA):
        super().__init__(base_pairs, RNA)
        self.bp1 = base_pairs[0]
        self.bp2 = base_pairs[1]
        self.subtype = self._infer_subtype()
        self.energy = self.calculate_energy()

    def _infer_subtype(self) -> InternalType:        
        if self.bp2.i - self.bp1.i == 2 and self.bp1.j - self.bp2.j == 2:
            return InternalType.INT11
        elif self.bp2.i - self.bp1.i == 2 and self.bp1.j - self.bp2.j == 3:
            return InternalType.INT12
        elif self.bp2.i - self.bp1.i == 3 and self.bp1.j - self.bp2.j == 2:
            return InternalType.INT21
        elif self.bp2.i - self.bp1.i == 2 or self.bp1.j - self.bp2.j == 2:
            return InternalType.INT1N        
        elif self.bp2.i - self.bp1.i == 3 and self.bp1.j - self.bp2.j == 3:
            return InternalType.INT22
        elif self.bp2.i - self.bp1.i == 3 and self.bp1.j - self.bp2.j == 4:
            return InternalType.INT23
        elif self.bp2.i - self.bp1.i == 4 and self.bp1.j - self.bp2.j == 3:
            return InternalType.INT32
        else:
            return InternalType.INTGEN

    def calculate_energy(self) -> int:
        mismatch1_nt1 = self.RNA[self.base_pairs[0].i]
        mismatch1_nt2 = self.RNA[self.base_pairs[0].j - 2]
        mismatch2_nt1 = self.RNA[self.base_pairs[1].i - 2]
        mismatch2_nt2 = self.RNA[self.base_pairs[1].j]

        common_term = initiation_df.loc[self.size, "internal"] + asymmetry * abs(self.bp2.i - self.bp1.i - self.bp1.j + self.bp2.j) if self.is_valid_size() else M        
        AU_closure_1 = self.bp1.nt1 + self.bp1.nt2 == 'AU' #or self.bp1.nt1 + self.bp1.nt2 == 'UA'
        GU_closure_1 = self.bp1.nt1 + self.bp1.nt2 == 'GU' #or self.bp1.nt1 + self.bp1.nt2 == 'UG'
        # AU_closure_2 = self.bp2.nt1 + self.bp2.nt2 == 'UA' # self.bp2.nt1 + self.bp2.nt2 == 'AU' or self.bp2.nt1 + self.bp2.nt2 == 'UA'  
        # GU_closure_2 = self.bp2.nt1 + self.bp2.nt2 == 'UG' #self.bp2.nt1 + self.bp2.nt2 == 'GU' or self.bp2.nt1 + self.bp2.nt2 == 'UG'

        # double_penalty = (AU_closure_1 or GU_closure_1) and (AU_closure_2 or GU_closure_2)
        penalty = (AU_closure_1 or GU_closure_1) #or (AU_closure_2 or GU_closure_2)

        if self.subtype == InternalType.INT11:
            G = int11_df[self.bp1.nt1 + self.bp1.nt2, mismatch1_nt1][self.bp2.nt1 + self.bp2.nt2, mismatch1_nt2]
        elif self.subtype == InternalType.INT12:
            G = int12_df.loc[self.bp1.nt1 + self.bp1.nt2, mismatch1_nt1][mismatch2_nt2, self.bp2.nt1 + self.bp2.nt2, mismatch1_nt2]
        elif self.subtype == InternalType.INT21:
            G = int12_df.loc[self.bp2.nt2 + self.bp2.nt1, mismatch2_nt2][mismatch1_nt1, self.bp1.nt2 + self.bp1.nt1, mismatch2_nt1]
        elif self.subtype == InternalType.INT1N:
            if penalty:
                G = common_term + AU_end_penalty
            # if double_penalty:
            #     G = common_term + 2*AU_end_penalty
            # elif penalty:
            #     G = common_term + AU_end_penalty
            else:
                G = common_term 
        elif self.subtype == InternalType.INT22:
            G = int22_df.loc[self.bp1.nt1 + self.bp1.nt2, mismatch1_nt1 + mismatch1_nt2][self.bp2.nt1 + self.bp2.nt2, mismatch2_nt1 + mismatch2_nt2]
        elif self.subtype == InternalType.INT23 or self.subtype == InternalType.INT32:
            int23_term = common_term + int23_df.loc[self.bp1.nt1 + self.bp1.nt2][mismatch1_nt1 + mismatch1_nt2] + int23_df.loc[self.bp2.nt2 + self.bp2.nt1][mismatch2_nt2 + mismatch2_nt1]
            if penalty:
                G = common_term + AU_end_penalty
            # if double_penalty:
            #     G = int23_term + 2*AU_end_penalty
            # elif penalty:
            #     G = int23_term + AU_end_penalty
            else:
                G = int23_term
        elif self.subtype == InternalType.INTGEN and self.is_valid_size():
            intgen_term = common_term + intnn_df.loc[self.bp1.nt1 + self.bp1.nt2][mismatch1_nt1 + mismatch1_nt2] + intnn_df.loc[self.bp2.nt1 + self.bp2.nt2][mismatch2_nt1 + mismatch2_nt2]
            if penalty:
                G = common_term + AU_end_penalty
            # if double_penalty:
            #     G = intgen_term + 2*AU_end_penalty
            # elif penalty:
            #     G = intgen_term + AU_end_penalty
            else:
                G = intgen_term
        else:
            G = M
        return round(G)

    def create_internal_size_constraint(self, model: gp.Model) -> None:        
        if not self.is_valid_size():
            inequality = gp.LinExpr([1], [self.var])
            model.addConstr(inequality == 0, f'IS-{self.base_pairs[0].i}-{self.base_pairs[0].j}-{self.base_pairs[1].i}-{self.base_pairs[1].j}')

    def create_internal_ifthen_constraint(self, model: gp.Model) -> None:
        bp1 = self.base_pairs[0]
        bp2 = self.base_pairs[1]

        inequality = gp.LinExpr(0)                       

        for u in range(bp1.i + 1, bp2.i):
            nucleotide = model.getVarByName(f'X_{u}')
            inequality.add(gp.LinExpr([1],[nucleotide]))

        for u in range(bp2.j + 1, bp1.j):
            nucleotide = model.getVarByName(f'X_{u}')
            inequality.add(gp.LinExpr([1],[nucleotide]))
            
        inequality.add(gp.LinExpr([1, 1, -1],[bp1.var, bp2.var, self.var]))
        model.addConstr(inequality <= self.size + 1, f'IIT-{bp1.i}-{bp1.j}-{bp2.i}-{bp2.j}')

    def create_internal_onlyif_constraint(self, model: gp.Model, base_pairs: List[BasePair]) -> None:
        bp1 = self.base_pairs[0]
        bp2 = self.base_pairs[1]

        for u in range(bp1.i + 1, bp2.i):
            inequality = gp.LinExpr(0)
            inequality.add(gp.LinExpr([3], [self.var]))

            matches = BasePair._find_base_pairs_with_index(base_pairs, u)
            
            for bp in matches:
                inequality.add(gp.LinExpr([1], [bp.var]))
            
            inequality.add(gp.LinExpr([-1, -1],[bp1.var, bp2.var]))
            model.addConstr(inequality <= 1, f'IOI-{bp1.i}-{bp1.j}-{bp2.i}-{bp2.j}-{u}')

        for u in range(bp2.j + 1, bp1.j):
            inequality = gp.LinExpr(0)
            inequality.add(gp.LinExpr([3], [self.var]))

            matches = BasePair._find_base_pairs_with_index(base_pairs, u)

            for bp in matches:
                inequality.add(gp.LinExpr([1], [bp.var]))

            inequality.add(gp.LinExpr([-1, -1], [bp1.var, bp2.var]))
            model.addConstr(inequality <= 1, f'IOI-{bp1.i}-{bp1.j}-{bp2.i}-{bp2.j}-{u}')

    def create_internal_max_number_constraint(model: gp.Model, internal_loops: List["InternalLoop"]) -> None:
        inequality = gp.LinExpr(0)

        for il in internal_loops:
            inequality.add(gp.LinExpr([1], [il.var]))
        model.addConstr(inequality == MAX_NUM_OF_LOOPS[il.type], f'IMN')
        model.update()

# rna = "GCCGCGAACCCCGCCAGGCCCGGAAGGGAGCAACGGUAGUGGUGGAU"
# bp1 = BasePair(14,30,rna)
# bp2 = BasePair(19,28,rna)
# i_loop = InternalLoop((bp1, bp2), rna)

# print(i_loop.energy)
# print(i_loop.subtype)
# print(i_loop.size)

