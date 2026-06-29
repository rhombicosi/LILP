from basepair import *
from lilp_config import *
from tnn_parameters.internal_parameters import *

class InternalBranch():

    def __init__(self, base_pairs: List[BasePair], RNA: str):
        self.RNA = RNA
        self.base_pairs = base_pairs
        self.bp1 = base_pairs[0]
        self.bp2 = base_pairs[1]
        self.type = LoopType.BRANCH
        self.distance = self.compute_distance()
        self.energy = self.calculate_energy()
        self.var = None

    def compute_distance(self):
        return self.bp2.i - self.bp1.j - 1
    
    def is_valid_size(self) -> bool:
        return self.distance <= MAX_LOOP_SIZES[self.type]
    
    def calculate_energy(self) -> int:
        if self.is_valid_size():
            G = SCALE * B
        else:
            G = M
        return round(G)
        # bp1 = self.bp1
        # bp2 = self.bp2
        # mismatch1_nt1 = self.RNA[self.base_pairs[0].i - 2]
        # mismatch1_nt2 = self.RNA[self.base_pairs[0].j]
        
        # if self.is_valid_size():
        #     if bp2.i - bp1.j == 1:
        #         G = SCALE * B + wcf_df.loc[bp1.nt1 + bp1.nt2, bp2.nt1 + bp2.nt2]
        #     elif bp2.i - bp1.j == 2:
        #         G = SCALE * B + wcf_df.loc[bp1.nt1 + bp1.nt2, bp2.nt1 + bp2.nt2] + intnn_df.loc[bp1.nt1 + bp1.nt2][mismatch1_nt1 + mismatch1_nt2]
        #     else:
        #         G = SCALE * B
        # else:
        #     G = M
        # return round(G)
    
    def _find_branches_with_pair(branches : List["InternalBranch"], pair: BasePair) -> List["InternalBranch"]:
        return [b for b in branches if (b.bp1.i == pair.i and b.bp1.j == pair.j) or (b.bp2.i == pair.i and b.bp2.j == pair.j)]
    
    def _find_branches_end_pair(branches : List["InternalBranch"], pair: BasePair) -> List["InternalBranch"]:
        return [b for b in branches if (b.bp2.i == pair.i and b.bp2.j == pair.j)]

    def create_branch_distance_constraint(self, model: gp.Model) -> None:
        n = len(self.RNA)      
        if self.distance > MAX_LOOP_SIZES[self.type] and self.bp1.i > BRANCH_START and self.bp2.j < n - BRANCH_START:
            inequality = gp.LinExpr([1], [self.var])
            model.addConstr(inequality == 0, f'BD-{self.base_pairs[0].i}-{self.base_pairs[0].j}-{self.base_pairs[1].i}-{self.base_pairs[1].j}')

    def create_branch_ifthen_constraint(self, model: gp.Model) -> None:
        bp1 = self.bp1
        bp2 = self.bp2
        inequality = gp.LinExpr(0)                       

        for u in range(bp1.j + 1, bp2.i):
            nucleotide = model.getVarByName(f'X_{u}')
            inequality.add(gp.LinExpr([1],[nucleotide]))
            
        inequality.add(gp.LinExpr([1, 1, -1],[bp1.var, bp2.var, self.var]))
        model.addConstr(inequality <= self.distance + 1, f'BIT-{bp1.i}-{bp1.j}-{bp2.i}-{bp2.j}')

    def create_branch_onlyif_constraint(self, model: gp.Model, base_pairs: List[BasePair]) -> None:
        bp1 = self.bp1
        bp2 = self.bp2

        if self.distance > 0:
            for u in range(bp1.j + 1, bp2.i):
                inequality = gp.LinExpr(0)
                inequality.add(gp.LinExpr([3], [self.var]))

                matches = BasePair._find_base_pairs_with_index(base_pairs, u)
                
                if matches:
                    for bp in matches:
                        inequality.add(gp.LinExpr([1], [bp.var]))
                    
                inequality.add(gp.LinExpr([-1, -1],[bp1.var, bp2.var]))
                model.addConstr(inequality <= 1, f'BOI-{bp1.i}-{bp1.j}-{bp2.i}-{bp2.j}-{u}') 
        else:
            inequality = gp.LinExpr(0)
            inequality.add(gp.LinExpr([3], [self.var]))
            inequality.add(gp.LinExpr([-1, -1],[bp1.var, bp2.var]))
            model.addConstr(inequality <= 1, f'BOI-{bp1.i}-{bp1.j}-{bp2.i}-{bp2.j}')        

    def create_branch_max_number_constraint(model: gp.Model, branches: List["InternalBranch"]) -> None:
        inequality = gp.LinExpr(0)

        for b in branches:
            inequality.add(gp.LinExpr([1], [b.var]))
        model.addConstr(inequality <= MAX_NUM_OF_LOOPS[b.type], f'BMN')
        model.update()

    def add_variable(self, model: gp.Model):
        indices = [str(idx) for bp in self.base_pairs for idx in (bp.i, bp.j)]
        index_str = "_".join(indices)
        
        name = f"{self.type.name}_{index_str}"
        self.var = model.addVar(vtype=GRB.BINARY, name=name)        
        return self.var 

class BranchPair(BasePair):
    def __init__(self, i, j, RNA):
        super().__init__(i, j, RNA)
        self.type = LoopType.BRANCHPAIR   
        self.energy = self.calculate_energy()

    def compute_distance(self):
        return self.i - self.j - 1
    
    def is_valid_size(self) -> bool:
        n = len(self.rna)
        return self.distance >= MIN_LOOP_SIZES[self.type] and self.i > BRANCH_START and self.j < n - BRANCH_START
    
    def calculate_energy(self) -> int:
        if self.is_valid_size():
            G = SCALE * B
        else:
            G = M
        return round(G)
    
    def _find_branchpairs_in_subsequence(branchpairs: List["BranchPair"], i: int, j: int) -> List["BranchPair"]:
        return [brp for brp in branchpairs if (brp.i > i and brp.j < j)]
    
    def create_branch_pair_ifthen_constraints(self, model: gp.Model, branches: List[InternalBranch]) -> None:
        bp_branches = InternalBranch._find_branches_with_pair(branches, self)
        if len(bp_branches) != 0:
            for bpb in bp_branches:
                model.addConstr(self.var >= bpb.var, f'BPIT-{self.i}-{self.j}-{bpb.bp1.i}-{bpb.bp1.j}-{bpb.bp2.i}-{bpb.bp2.j}')
        else:
            model.addConstr(self.var == 0)

    def create_branch_pair_onlyif_constraints(self, model: gp.Model, branches: List[InternalBranch]) -> None:
        bp_branches = InternalBranch._find_branches_with_pair(branches, self)
        if len(bp_branches) != 0:            
            model.addConstr(self.var <= gp.quicksum(bpb.var for bpb in bp_branches), f'BPOI-{self.i}-{self.j}')                
        else:
            model.addConstr(self.var == 0)

class ClosingBranch():
    def __init__(self, base_pairs: List[BasePair], RNA: str):
        self.RNA = RNA
        self.base_pairs = base_pairs
        self.bp1 = base_pairs[0]
        self.bp2 = base_pairs[1]        
        self.type = LoopType.CBRANCH
        self.distance = self.compute_distance()
        self.energy = self.calculate_energy()
        self.var = None

    def compute_distance(self):
        bp1 = self.bp1
        bp2 = self.bp2
        return bp1.j - bp2.j - 1
    
    def is_valid_size(self) -> bool:
        n = len(self.RNA)
        return self.distance <= MAX_LOOP_SIZES[self.type] and self.bp1.j - self.bp1.i > MIN_LOOP_SIZES[self.type] and self.bp1.i > BRANCH_START and self.bp1.j < n - BRANCH_START
    
    def add_variable(self, model: gp.Model):
        bp1 = self.bp1
        bp2 = self.bp2
        indices = [bp1.i, bp1.j, bp2.i, bp2.j]
        indices = list(map(str, indices))
        index_str = "_".join(indices)
        
        name = f"CBRANCH_{index_str}"
        self.var = model.addVar(vtype=GRB.BINARY, name=name)
        return self.var
    
    def calculate_energy(self) -> int:
        if self.is_valid_size():
            G = SCALE * A
        else:
            G = M
        return round(G)
    
        # bp1 = self.bp1
        # bp2 = self.bp2
        # mismatch1_nt1 = self.RNA[self.base_pairs[0].i]
        # mismatch1_nt2 = self.RNA[self.base_pairs[0].j - 2]

        # if self.is_valid_size():
        #     if bp1.j - bp2.j == 1:
        #         G = G = SCALE * A + wcf_df.loc[bp1.nt1 + bp1.nt2, bp2.nt1 + bp2.nt2]
        #     elif bp1.j - bp2.j == 2:
        #         G = G = SCALE * A + wcf_df.loc[bp1.nt1 + bp1.nt2, bp2.nt1 + bp2.nt2] + intnn_df.loc[bp1.nt1 + bp1.nt2][mismatch1_nt1 + mismatch1_nt2]
        #     else:
        #         G = SCALE * A
        # else:
        #     G = M
        # return round(G)
    
    def create_branch_distance_constraint(self, model: gp.Model) -> None: 
        bp1 = self.bp1
        bp2 = self.bp2       
        if not self.is_valid_size():
            inequality = gp.LinExpr([1], [self.var])
            model.addConstr(inequality == 0, f'CBD-{bp1.i}-{bp1.j}-{bp2.i}-{bp2.j}')

    # def create_auxiliary_constraints(self, model: gp.Model, branch_pairs: List["BranchPair"] ) -> None:
    #     bp1 = self.bp1
    #     bp2 = self.bp2 
        
    #     aux = model.getVarByName(f'Y_{bp1.i}_{bp1.j}_{bp2.i}_{bp2.j}')
    #     matches = BranchPair._find_branchpairs_in_subsequence(branch_pairs, bp1.i, bp2.i)
    #     if matches:
    #         for m in matches:
    #             model.addConstr(aux >= m.var, f'AUX-{bp1.i}-{bp1.j}-{bp2.i}-{bp2.j}-{m.i}-{m.j}')
    
    def create_closing_branch_ifthen_constraint(self, model: gp.Model) -> None:
        bp1 = self.bp1
        bp2 = self.bp2 
        inequality = gp.LinExpr(0)

        bpvar = model.getVarByName(f'BP_{bp2.i}_{bp2.j}')

        for u in range(bp2.j + 1, bp1.j):
            nucleotide = model.getVarByName(f'X_{u}')
            inequality.add(gp.LinExpr([1],[nucleotide]))

        inequality.add(gp.LinExpr([1, 1, 1, -1],[bp1.var, bp2.var, bpvar, self.var]))
        model.addConstr(inequality <= self.distance + 2, f'CBIT-{bp1.i}-{bp1.j}-{bp2.i}-{bp2.j}')

    # def create_closing_branch_ifthen_constraint(self, model: gp.Model) -> None:
    #     bp1 = self.bp1
    #     bp2 = self.bp2 
    #     inequality = gp.LinExpr(0)

    #     bp = model.getVarByName(f'BP_{bp2.i}_{bp2.j}')

    #     for u in range(bp2.j + 1, bp1.j):
    #         nucleotide = model.getVarByName(f'X_{u}')
    #         inequality.add(gp.LinExpr([1],[nucleotide]))

    #     inequality.add(gp.LinExpr([1, 1, 1, -1],[bp1.var, bp2.var, bp, self.var]))
    #     model.addConstr(inequality <= self.distance + 2, f'CBIT-{bp1.i}-{bp1.j}-{bp2.i}-{bp2.j}')

    

# class ClosingBranch():
#     def __init__(self, base_pair: BasePair, i: int, j: int, RNA: str):
#         self.RNA = RNA
#         self.base_pair = base_pair
#         self.i = i
#         self.j = j
#         self.distance = self.compute_distance()
#         self.energy = self.calculate_energy()
#         # self.type = LoopType.BRANCH
#         self.var = None

#     def compute_distance(self):
#         return self.i - self.base_pair.i + self.base_pair.j - self.j
    
#     def add_variable(self, model: gp.Model):
#         indices = [self.base_pair.i, self.i, self.j, self.base_pair.j]
#         indices = list(map(str, indices))
#         index_str = "_".join(indices)
        
#         name = f"CBRANCH_{index_str}"
#         self.var = model.addVar(vtype=GRB.BINARY, name=name)        
#         return self.var
    
#     def calculate_energy(self) -> int:
#         G = SCALE * B
#         return round(G)
    
#     def create_branch_distance_constraint(self, model: gp.Model) -> None:        
#         if self.distance > MAX_LOOP_SIZES[self.type]:
#             inequality = gp.LinExpr([1], [self.var])
#             model.addConstr(inequality == 0, f'CBD-{self.base_pair.i}-{self.i}-{self.j}-{self.base_pair.j}')
    
#     def create_closing_branch_ifthen_constraint(self, branch_pairs: List["BranchPair"], model: gp.Model) -> None:
#         bp = self.base_pair
#         inequality = gp.LinExpr(0)

#         matches_left = BranchPair._find_branchpairs_startwith_nt(branch_pairs, self.i)
#         matches_right = BranchPair._find_branchpairs_endwith_nt(branch_pairs, self.j)
        
#         inequality.add(gp.LinExpr([1], [bp.var]))

#         if matches_left:
#             for ml in matches_left:
#                 inequality.add(gp.LinExpr([1.0], [ml.var]))

#         if matches_right:
#             for mr in matches_right:
#                 inequality.add(gp.LinExpr([1.0], [mr.var]))

#         for u in range(bp.i + 1, self.i):
#             nucleotide = model.getVarByName(f'X_{u}')
#             inequality.add(gp.LinExpr([1],[nucleotide]))

#         for u in range(self.j + 1, bp.j):
#             nucleotide = model.getVarByName(f'X_{u}')
#             inequality.add(gp.LinExpr([1],[nucleotide]))

#         inequality.add(gp.LinExpr([-1],[self.var]))
#         model.addConstr(inequality <= self.distance, f'CBIT-{bp.i}-{self.i}-{self.j}-{bp.j}')

# rna = "GCCGCGAACCCCGCCAGGCCCGGAAGGGAGCAACGGUAGUGGUGGAU"
# bp = BasePair(10,39,rna)
# i = 14
# j = 38

# cb = ClosingBranch(bp,i,j,rna)
# model = gp.Model("probe")
# cb.add_variable(model)
# model.update()

# print(cb.distance)
# print(cb.var)
