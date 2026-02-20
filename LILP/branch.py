from basepair import *
from lilp_config import *

class Branch():

    def __init__(self, base_pairs: List[BasePair], RNA: str):
        self.RNA = RNA
        self.base_pairs = base_pairs
        self.bp1 = base_pairs[0]
        self.bp2 = base_pairs[1] 
        self.distance = self.compute_distance()
        # self.energy = self.calculate_energy()
        self.type = LoopType.BRANCH
        self.var = None

    def compute_distance(self):
        return self.bp2.i - self.bp1.j - 1
    
    # def is_valid_size(self) -> bool:
    #     return self.distance >= MAX_LOOP_SIZES[self.type]
    
    # def calculate_energy(self) -> int:
    #     if self.is_valid_size():
    #         G = SCALE * BB
    #     else:
    #         G = M
    #     return round(G)
    
    def _find_branches_with_pair(branches : List["Branch"], pair: BasePair) -> List["Branch"]:
        return [b for b in branches if (b.bp1.i == pair.i and b.bp1.j == pair.j) or (b.bp2.i == pair.i and b.bp2.j == pair.j)]

    def create_branch_distance_constraint(self, model: gp.Model) -> None:        
        if self.distance > BRANCH_D:
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

        for u in range(bp1.j + 1, bp2.i):
            inequality = gp.LinExpr(0)
            inequality.add(gp.LinExpr([3], [self.var]))

            matches = BasePair._find_base_pairs_with_index(base_pairs, u)
            
            for bp in matches:
                inequality.add(gp.LinExpr([1], [bp.var]))
            
            inequality.add(gp.LinExpr([-1, -1],[bp1.var, bp2.var]))
            model.addConstr(inequality <= 1, f'BOI-{bp1.i}-{bp1.j}-{bp2.i}-{bp2.j}-{u}')        

    def create_branch_max_number_constraint(model: gp.Model, branches: List["Branch"]) -> None:
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
        self.type = LoopType.BRANCH   
        self.energy = self.calculate_energy()

    def compute_distance(self):
        return self.i - self.j - 1
    
    def is_valid_size(self) -> bool:
        return self.distance >= MAX_LOOP_SIZES[self.type]
    
    def calculate_energy(self) -> int:
        if self.is_valid_size():
            G = SCALE * BB
        else:
            G = M
        return round(G)

    # def add_variable(self, model: gp.Model):
    #     indices = [self.i, self.j]
    #     index_str = "_".join(indices)
        
    #     name = f"BP_{index_str}"
    #     self.var = model.addVar(vtype=GRB.BINARY, name=name)        
    #     return self.var 


# rna = "GCCGCGAACCCCGCCAGGCCCGGAAGGGAGCAACGGUAGUGGUGGAU"
# bp1 = BasePair(14,20,rna)
# bp2 = BasePair(25,38,rna)
# branch = Branch((bp1, bp2), rna)

# print(branch.distance)
# print(branch.type)
# print(branch.var)
