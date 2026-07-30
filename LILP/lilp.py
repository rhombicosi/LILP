import gurobipy as gp
from gurobipy import GRB
from utils.prepro_utils import *
from utils.sol_converter import *
from lilp_config import *
from basepair import *
from dloop import *
from stemloop import *
from hairpinloop import *
from internalloop import *
from bulgeloop import *
from multiloop import *
from branch import *
from subsequence import *

class LILP:
    def __init__(self, rna_seq: str, name: str):
        self.rna_seq = rna_seq
        self.model = gp.Model(name)
        self.nucleotides : List[gp.Var] = []
        self.base_pairs : List[BasePair] = []
        self.not_pairs : List[gp.Var] = []
        self.subsequences: List[Subsequence] = []
        self.first_pairs : List[BasePair] = []
        self.last_pairs : List[BasePair] = []
        self.hairpin_loops : List[HairpinLoop] = []
        self.stem_loops : List[StemLoop] = []
        self.internal_loops : List[InternalLoop] = []
        self.bulge_loops : List[BulgeLoop] = []
        self.multi_loops : List[MultiLoop] = []
        self.branches: List[InternalBranch] = []
        self.branch_pairs: List[BranchPair] = []
        self.cbranches: List[ClosingBranch] = []
        self.objective : gp.LinExpr   

    def create_nucleotides(self, first, last) -> None:
        for i in range(first, last + 1):
            # var = self.model.addVar(vtype=GRB.BINARY, name=f'X_{i}')
            var = self.model.addVar(vtype=GRB.BINARY, name=f'Y_{i}')
            var = self.model.addVar(vtype=GRB.BINARY, name=f'Z_{i}') 
            self.nucleotides.append(var)
        self.model.update()

    def create_base_pairs(self, first, last) -> None:
        for i in range(first, last): 
            for j in range(i + MIN_D + 1, last + 1):
                bp = BasePair(i, j, self.rna_seq)
                var = bp.add_variable(self.model, 'P')
                if var:
                    self.base_pairs.append(bp)
        self.model.update()

    def create_not_pairs(self) -> None:        
        for bp in self.base_pairs:
            # var = self.model.addVar(vtype=GRB.BINARY, name=f'X_{bp.i}_{bp.j}')
            var = self.model.addVar(vtype=GRB.BINARY, name=f'Y_{bp.i}_{bp.j}') 
            var = self.model.addVar(vtype=GRB.BINARY, name=f'Z_{bp.i}_{bp.j}')  
            self.not_pairs.append(var)
        self.model.update()

    def create_subsequences(self) -> None:
        ssq = set()
        for hl in self.hairpin_loops:
            ssq.add((hl.base_pairs[0].i, hl.base_pairs[0].j))

        # for il in self.internal_loops:
        #     ssq.add((il.base_pairs[0].i, il.base_pairs[1].i))
        #     ssq.add((il.base_pairs[1].j, il.base_pairs[0].j))

        for bl in self.bulge_loops:
            bp1 = bl.base_pairs[0]
            bp2 = bl.base_pairs[1]

            if bp2.i == bp1.i + 1:
                ssq.add((bp2.j, bp1.j))
            elif bp2.j == bp1.j - 1:
                ssq.add((bp1.i, bp2.i))

        # for br in self.branches:
        #     bp1 = br.base_pairs[0]
        #     bp2 = br.base_pairs[1]

        #     ssq.add((br.bp1.j,br.bp2.i))

        for s,e in ssq:
            sq = Subsequence(s, e)
            sq.add_variable(self.model)
            self.subsequences.append(sq)
        self.model.update()

    def create_first_pairs(self) -> None:
        for bp in self.base_pairs:
            fp = BasePair(bp.i, bp.j, self.rna_seq)
            self.first_pairs.append(fp)
            fp.add_variable(self.model, 'F') 
        self.model.update()

    def create_last_pairs(self) -> None:
        for bp in self.base_pairs:
            lp = BasePair(bp.i, bp.j, self.rna_seq)
            self.last_pairs.append(lp)
            lp.add_variable(self.model,'L')
        self.model.update()

    def create_hairpin_loops(self) -> None:
        for bp in self.base_pairs:
            hairpin = HairpinLoop([bp], self.rna_seq)
            # if hairpin.size < U_MAX:
            hairpin.add_variable(self.model)
            self.hairpin_loops.append(hairpin)
        self.model.update()

    def create_stem_loops(self) -> None:
        for bp1 in self.base_pairs:
            bp2 = BasePair._find_base_pairs_matches(self.base_pairs, bp1.i + 1, bp1.j - 1)
            if bp2:
                stem = StemLoop([bp1, bp2], self.rna_seq)
                stem.add_variable(self.model)
                self.stem_loops.append(stem)
        self.model.update()

    def create_internal_loops(self) -> None:
        for bp1 in self.base_pairs:
            for bp2 in self.base_pairs:
                if bp2.i > bp1.i + 1 and bp2.j < bp1.j - 1:
                    internal = InternalLoop([bp1, bp2], self.rna_seq)
                    # if internal.bp2.i - internal.bp1.i - 1 < U_MAX and internal.bp1.j - internal.bp2.j - 1 < U_MAX:
                    internal.add_variable(self.model)
                    self.internal_loops.append(internal)
        self.model.update()

    def create_bulge_loops(self) -> None:
        for bp1 in self.base_pairs:
            for bp2 in self.base_pairs:
                if (bp2.i == bp1.i + 1 and bp2.j < bp1.j - 1) or (bp2.i > bp1.i + 1 and bp2.j == bp1.j - 1):
                    bulge = BulgeLoop([bp1, bp2], self.rna_seq)
                    # if bulge.size < U_MAX:
                    bulge.add_variable(self.model)
                    self.bulge_loops.append(bulge)
        self.model.update()

    def create_multi_loops(self) -> None:
        n = len(self.rna_seq)
        for bp1 in self.base_pairs:
            for bp2 in self.base_pairs:
                for bp3 in self.base_pairs:
                    if bp2.i > bp1.i and bp3.i > bp2.j and bp1.j > bp3.j:
                        if bp1.distance > MULTI_MIN_D and bp2.distance > MULTI_MIN_D and bp3.distance > MULTI_MIN_D:
                            if bp1.i >= MULTI_BOUND and bp1.j <= n - MULTI_BOUND:
                                multi = MultiLoop([bp1, bp2, bp3], self.rna_seq)
                                multi.add_variable(self.model)
                                self.multi_loops.append(multi)
        self.model.update()

    def create_branches(self) -> None:
        n = len(self.rna_seq)
        for bp1 in self.base_pairs:
            for bp2 in self.base_pairs:
                if bp2.i > bp1.j:
                    # if bp1.j - bp1.i <= n/2 and bp2.j - bp2.i <= n/2:
                    branch = InternalBranch([bp1, bp2], self.rna_seq)
                    # if branch.distance < U_MAX:
                    branch.add_variable(self.model)
                    self.branches.append(branch)
        self.model.update()

    def create_branch_pairs(self):
        for bp in self.base_pairs:
            branch_pair = BranchPair(bp.i, bp.j, self.rna_seq)
            branch_pair.add_variable(self.model, 'BP')
            self.branch_pairs.append(branch_pair)
        self.model.update()

    def create_closing_branches(self) -> None:
        for bp1 in self.base_pairs:
            for bp2 in self.base_pairs:
                if bp1.i < bp2.i and bp2.j < bp1.j and (bp2.i - bp1.i - 1 >= MIN_D + 2):
                    branch = ClosingBranch([bp1, bp2], self.rna_seq)
                    # if branch.distance < U_MAX:
                    branch.add_variable(self.model)
                    self.cbranches.append(branch)
        self.model.update()

    # def create_multiloop(self) -> None:
    #     bp1 = BasePair._find_base_pair_with_indices(self.base_pairs,5,49)[0]
    #     bp2 = BasePair._find_base_pair_with_indices(self.base_pairs,10,22)[0]
    #     bp3 = BasePair._find_base_pair_with_indices(self.base_pairs,24,40)[0]
    #     multi = MultiLoop([bp1,bp2,bp3], self.rna_seq)
    #     multi.add_variable(self.model)
    #     self.multi_loops.append(multi)
    #     self.model.update()

    ######################## ADD CONSTRAINTS ###########################

    def add_unpaired_nucleotides_constraints(self, first, last) -> None:
        for i in range(first, last + 1):
            # inequality1 = gp.LinExpr(0)
            inequality2 = gp.LinExpr(0)
            inequality3 = gp.LinExpr(0)
            matches = BasePair._find_base_pairs_with_index(self.base_pairs, i)

            if matches:
                for bp in matches:
                    # inequality1.add(gp.LinExpr([1], [bp.var]))
                    inequality2.add(gp.LinExpr([1], [bp.var]))
                    inequality3.add(gp.LinExpr([1], [bp.var]))
                
                # nucleotide = self.model.getVarByName(f'X_{i}')
                # inequality1.add(gp.LinExpr([1], [nucleotide]))
                # self.model.addConstr(inequality1 == 1, f'UNX_{i}')

                nucleotide = self.model.getVarByName(f'Y_{i}')
                inequality2.add(gp.LinExpr([1], [nucleotide]))
                self.model.addConstr(inequality2 == 1, f'UNY_{i}')

                nucleotide = self.model.getVarByName(f'Z_{i}')
                inequality3.add(gp.LinExpr([1], [nucleotide]))
                self.model.addConstr(inequality3 == 1, f'UNZ_{i}')
        self.model.update()

    def add_max_unpaired_region_constraints(self) -> None:
        n = len(self.rna_seq)
        for i in range(1, n - U_MAX + 2):
            inequality = gp.LinExpr(0)
            for j in range(i, i + U_MAX):
                nucleotide = self.model.getVarByName(f'X_{j}')
                inequality.add(gp.LinExpr([1], [nucleotide]))
            self.model.addConstr(inequality <= U_MAX - 1, f'UNR_{i}')
        self.model.update()
    
    def add_single_pair_constraints(self, first, last) -> None:
        for i in range(first, last + 1):
            BasePair.create_single_pair_constraint(self.model, self.base_pairs, i)
        self.model.update()

    def add_no_crossing_constraints(self) -> None: 
        for bp1 in self.base_pairs: 
            for bp2 in self.base_pairs:
                BasePair.create_no_crossing_constraint(self.model, bp1, bp2)
        self.model.update()

    def add_not_pairs_constraints(self) -> None: 
        for bp in self.base_pairs:
            # inequality = gp.LinExpr(0)
            # nbp = self.model.getVarByName(f'X_{bp.i}_{bp.j}')
            # inequality.add(gp.LinExpr([1.0, 1.0], [nbp, bp.var]))
            # self.model.addConstr(inequality == 1, f'NXBP_{bp.i}_{bp.j}')

            inequality = gp.LinExpr(0)
            nbp = self.model.getVarByName(f'Y_{bp.i}_{bp.j}')
            inequality.add(gp.LinExpr([1.0, 1.0], [nbp, bp.var]))
            self.model.addConstr(inequality == 1, f'NYBP_{bp.i}_{bp.j}')

            inequality = gp.LinExpr(0)
            nbp = self.model.getVarByName(f'Z_{bp.i}_{bp.j}')
            inequality.add(gp.LinExpr([1.0, 1.0], [nbp, bp.var]))
            self.model.addConstr(inequality == 1, f'NZBP_{bp.i}_{bp.j}')
        self.model.update()  

    def add_subsequences_constraints(self) -> None:
        for ssq in self.subsequences:
            ssq.create_unpaired_subsequence_constaint(self.model, self.base_pairs)
        self.model.update()

    def add_stem_ifthen_constraints(self) -> None:
        for sl in self.stem_loops:
            sl.create_stem_ifthen_constraint(self.model)
        self.model.update()

    def add_stem_onlyif_constraints(self) -> None:
        for sl in self.stem_loops:
            sl.create_stem_onlyif_constraint(self.model)
        self.model.update()

    def add_stem_constraints(self) -> None:
        for sl in self.stem_loops:
            if sl.energy > 0:
                sl.create_stem_ifthen_constraint(self.model)
            else:
                sl.create_stem_onlyif_constraint(self.model)
            # sl.create_stem_constraints(self.model)
        self.model.update()
        
    def add_first_pair_constraints(self) -> None:
        for sl in self.stem_loops:
            sl.create_first_pair_constraints(self.model, self.stem_loops, self.first_pairs)
        self.model.update()

    def add_last_pair_constraints(self) -> None:
        for sl in self.stem_loops:    
            sl.create_last_pair_constraints(self.model, self.stem_loops, self.last_pairs)
        self.model.update()

    def add_hairpin_size_constraints(self) -> None:
        for hl in self.hairpin_loops:
            hl.create_hairpin_size_constraint(self.model)
        self.model.update()

    def add_hairpin_constraints(self) -> None:
        for hl in self.hairpin_loops:
            if hl.energy > 0:
                hl.create_hairpin_ifthen_constraint(self.model, self.base_pairs)
            else:
                hl.create_hairpin_onlyif_constraint(self.model, self.base_pairs)

    def add_hairpin_max_number_constraint(self) -> None:
        HairpinLoop.create_hairpin_max_number_constraint(self.model, self.hairpin_loops)

    def add_internal_size_constraints(self) -> None:
        for il in self.internal_loops:
            il.create_internal_size_constraint(self.model)
        self.model.update()  
    
    def add_internal_constraints(self)-> None:        
        for il in self.internal_loops:                
            if il.energy > 0:
                il.create_internal_ifthen_constraint(self.model)
            else:
                il.create_internal_onlyif_constraint(self.model, self.base_pairs)

    def add_internal_max_number_constraint(self) -> None:
        InternalLoop.create_internal_max_number_constraint(self.model, self.internal_loops)

    def add_bulge_size_constraints(self) -> None:
        for il in self.bulge_loops:
            il.create_bulge_size_constraint(self.model)
        self.model.update()

    def add_bulge_constraints(self) -> None:
        for bl in self.bulge_loops:
            if bl.energy > 0:
                bl.create_bulge_ifthen_constraint(self.model)
            else:
                bl.create_bulge_onlyif_constraint(self.model, self.base_pairs)

    def add_bulge_max_number_constraint(self) -> None:
        BulgeLoop.create_bulge_max_number_constraint(self.model, self.bulge_loops)

    def add_multi_size_constraints(self) -> None:
        for ml in self.multi_loops:
            ml.create_multi_size_constraint(self.model)
        self.model.update()

    def add_multi_energy_constraints(self, energy) -> None:
        for ml in self.multi_loops:
            ml.create_multi_energy_constraint(self.model, energy)
        self.model.update()
    
    def add_multi_constraints(self) -> None:
        for ml in self.multi_loops:
            if ml.energy > 0:
                ml.create_multi_ifthen_constraint(self.model)
            else:
                ml.create_multi_onlyif_constraint(self.model, self.base_pairs)
        self.model.update()

    def add_multi_max_number_constraint(self) -> None:
        MultiLoop.create_multi_max_number_constraint(self.model, self.multi_loops)

    def add_branch_distance_constraints(self) -> None:
        for b in self.branches:
            b.create_branch_distance_constraint(self.model)
        self.model.update()

    def add_branch_constraints(self) -> None:
        for b in self.branches:
            if b.energy > 0:
                b.create_branch_ifthen_constraint(self.model)
            else:
                b.create_branch_onlyif_constraint(self.model, self.base_pairs)
        self.model.update()

    def add_branch_pair_distance_constraints(self) -> None:
        for b in self.branch_pairs:
            b.create_branchpair_distance_constraint(self.model)
        self.model.update()

    def add_branch_pairs_constraints(self) -> None:
        for bp in self.branch_pairs:
            # if bp.energy > 0:
                bp.create_branch_pair_ifthen_constraints(self.model, self.branches)
            # else:
            #     bp.create_branch_pair_onlyif_constraints(self.model, self.branches)
        self.model.update()

    def add_branch_base_pairs_constraints(self) -> None:
        for bp in self.base_pairs:
            branch_pair = self.model.getVarByName(f'BP_{bp.i}_{bp.j}')            
            self.model.addConstr(branch_pair <= bp.var , f'BBP-{bp.i}-{bp.j}')

    def add_cbranch_constraints(self) -> None:
        for cb in self.cbranches:
            if cb.energy > 0:
                cb.create_closing_branch_ifthen_constraint(self.model)
            else:
                cb.create_closing_branch_onlyif_constraint(self.model, self.base_pairs)
        self.model.update()

    def add_cbranch_distance_constraints(self) -> None:
        for cb in self.cbranches:
            cb.create_branch_distance_constraint(self.model)
        self.model.update()

    def add_energy_constraint(self, stem, hairpin, internal, bulge, branch, cbranch) -> gp.LinExpr:
        inequality = gp.LinExpr()
        if stem:
            inequality.add(self.create_stem_term())
        if hairpin:
            inequality.add(self.create_hairpin_term())
        if internal:
            inequality.add(self.create_internal_term())
        if bulge:
            inequality.add(self.create_bulge_term())
        if branch:
            # inequality.add(self.create_branch_term())
            inequality.add(self.create_branchpair_term())
        if cbranch:
            inequality.add(self.create_cbranch_term())        
        self.model.addConstr(inequality >= MFE, f'MFE')
        self.model.update()

    def add_lowest_energy_constraint(self, stem, hairpin, internal, bulge, branch, cbranch) -> gp.LinExpr:
        inequality = gp.LinExpr()
        stem_energies = []
        for sl in self.stem_loops:
            stem_energies.append(sl.energy)

        lowest_bound = min(stem_energies) * (len(self.rna_seq)/2 - 3)
        if stem:
            inequality.add(self.create_stem_term())
        if hairpin:
            inequality.add(self.create_hairpin_term())
        if internal:
            inequality.add(self.create_internal_term())
        if bulge:
            inequality.add(self.create_bulge_term())
        if branch:
            # inequality.add(self.create_branch_term())
            inequality.add(self.create_branchpair_term())
        if cbranch:
            inequality.add(self.create_cbranch_term())        
        self.model.addConstr(inequality >= lowest_bound, f'LOWEST_BOUND')
        self.model.update()

    ###################### CREATE OBJECTIVE TERMS ################################

    def create_stem_term(self) -> gp.LinExpr:
        objective = gp.LinExpr([sl.energy for sl in self.stem_loops], [sl.var for sl in self.stem_loops])
        return objective
    
    # def create_hairpin_term(self) -> gp.LinExpr:
    #     objective = gp.LinExpr([hl.energy for hl in self.hairpin_loops], [hl.var for hl in self.hairpin_loops])
    #     return objective
    
    def create_hairpin_term(self) -> gp.LinExpr:
        pairs = [(hl.energy, hl.var) for hl in self.hairpin_loops if hl.energy < M]
        return gp.LinExpr([e for e, _ in pairs], [v for _, v in pairs])
    
    # def create_internal_term(self) -> gp.LinExpr:
    #     objective = gp.LinExpr([il.energy for il in self.internal_loops], [il.var for il in self.internal_loops])
    #     return objective
    
    def create_internal_term(self) -> gp.LinExpr:
        pairs = [(il.energy, il.var) for il in self.internal_loops if il.energy < M]
        return gp.LinExpr([e for e, _ in pairs], [v for _, v in pairs])    
    
    # def create_bulge_term(self) -> gp.LinExpr:
    #     objective = gp.LinExpr([bl.energy for bl in self.bulge_loops], [bl.var for bl in self.bulge_loops])
    #     return objective
    
    def create_bulge_term(self) -> gp.LinExpr:
        pairs = [(bl.energy, bl.var) for bl in self.bulge_loops if bl.energy < M]
        return gp.LinExpr([e for e, _ in pairs], [v for _, v in pairs])
    
    def create_multi_term(self) -> gp.LinExpr:
        objective = gp.LinExpr([ml.energy for ml in self.multi_loops], [ml.var for ml in self.multi_loops])
        return objective
    
    # def create_branchpair_term(self) -> gp.LinExpr:
    #     objective = gp.LinExpr([bp.energy for bp in self.branch_pairs], [bp.var for bp in self.branch_pairs])
    #     return objective
    
    def create_branchpair_term(self) -> gp.LinExpr:
        pairs = [(bp.energy, bp.var) for bp in self.branch_pairs if bp.energy < M]
        return gp.LinExpr([e for e, _ in pairs], [v for _, v in pairs])
    
    # def create_branch_term(self) -> gp.LinExpr:
    #     objective = gp.LinExpr([bp.energy for bp in self.branches], [bp.var for bp in self.branches])
    #     return objective
    
    def create_branch_term(self) -> gp.LinExpr:
        pairs = [(bp.energy, bp.var) for bp in self.branches if bp.energy < M]
        return gp.LinExpr([e for e, _ in pairs], [v for _, v in pairs])
    
    # def create_cbranch_term(self) -> gp.LinExpr:
    #     objective = gp.LinExpr([cbp.energy for cbp in self.cbranches], [cbp.var for cbp in self.cbranches])
    #     return objective
    
    def create_cbranch_term(self) -> gp.LinExpr:
        pairs = [(cbp.energy, cbp.var) for cbp in self.cbranches if cbp.energy < M]
        return gp.LinExpr([e for e, _ in pairs], [v for _, v in pairs])

    ################ MODEL CONSCTRUCTION ###################################     
    
    def create_objective(self, stem, hairpin, internal, bulge, branch, cbranch) -> gp.LinExpr:
        objective = gp.LinExpr()
        if stem:
            objective.add(self.create_stem_term())
        if hairpin:
            objective.add(self.create_hairpin_term())
        if internal:
            objective.add(self.create_internal_term())
        if bulge:
            objective.add(self.create_bulge_term())
        if branch:
            # objective.add(self.create_branch_term())
            objective.add(self.create_branchpair_term())
        if cbranch:
            objective.add(self.create_cbranch_term())
        self.model.setObjective(objective, GRB.MINIMIZE)
        self.model.update()

    def create_variables(self, stem, hairpin, internal, bulge, branch, cbranch, first, last):
        self.create_base_pairs(first, last)                
        # self.create_not_pairs()                       
        self.create_nucleotides(first, last)
        if hairpin:
            self.create_hairpin_loops()       
        if stem:            
            self.create_stem_loops() 
        if internal:  
            self.create_internal_loops()
        if bulge:
            self.create_bulge_loops()          
        if branch:
            self.create_branches()
            self.create_branch_pairs()      
        if cbranch:
            self.create_closing_branches()         
        self.create_subsequences()

        print(f"PAIRS: {len(self.base_pairs)}, NUCLEOTIDES: {len(self.nucleotides)}, SSQ: {len(self.subsequences)}, HAIRPIN: {len(self.hairpin_loops)}, STEMS: {len(self.stem_loops)}, INTERNALS: {len(self.internal_loops)}, BULGES: {len(self.bulge_loops)}, BRANCHES: {len(self.branches)}, BRANCH PAIRS: {len(self.branch_pairs)}, CLOSING_BRANCHES: {len(self.cbranches)}")    
    
    def create_constraints(self, stem, hairpin, internal, bulge, branch, cbranch, first, last):
        self.add_single_pair_constraints(first, last)
        self.add_no_crossing_constraints()
        if stem:
            # self.add_stem_ifthen_constraints()
            # self.add_stem_onlyif_constraints()
            self.add_stem_constraints()
        self.add_unpaired_nucleotides_constraints(first, last)
        # self.add_max_unpaired_region_constraints()
        # self.add_not_pairs_constraints()
        self.add_subsequences_constraints()
        if hairpin:
            self.add_hairpin_size_constraints()
            self.add_hairpin_constraints()
            # self.add_hairpin_ifthen_constraints()
            # self.add_hairpin_onlyif_constraints()            
            # self.add_hairpin_max_number_constraint()
        if internal:
            self.add_internal_size_constraints()
            self.add_internal_constraints()
            # self.add_internal_onlyif_constraints()
            # self.add_internal_ifthen_constraints()            
            # self.add_internal_max_number_constraint()
            # self.model.addConstr(self.model.getVarByName(f'INTERNAL_5_39_11_36') == 1)
        if bulge:
            self.add_bulge_size_constraints()
            self.add_bulge_constraints()
            # self.add_bulge_ifthen_constraints()
            # self.add_bulge_onlyif_constraints()
            # self.add_bulge_max_number_constraint()
        if branch: 
            self.add_branch_distance_constraints()
            self.add_branch_constraints()
            self.add_branch_pair_distance_constraints()
            self.add_branch_pairs_constraints()           
            self.add_branch_base_pairs_constraints()
        if cbranch:
            self.add_cbranch_distance_constraints()
            self.add_cbranch_constraints()
        # self.add_energy_constraint(stem, hairpin, internal, bulge, branch, cbranch)
        # self.add_lowest_energy_constraint(stem, hairpin, internal, bulge, branch, cbranch)

# seq_number = 1
# chain_file = seq_files[seq_number]
# chain_name_with_ext = os.path.basename(chain_file)        
# chain_name_without_ext = os.path.splitext(chain_name_with_ext)[0]
# lp_file_name = chain_name_without_ext
# seq_data = parse_seq_file(chain_file)

# rna = seq_data['sequence'].upper()
# print(rna)
# print(len(rna))
# model_name = 'test'
# rna_model = LILP(rna, model_name)

# first = 1
# last = len(rna)

# rna_model.create_variables(0, 0, 0, 0, 0, 0, first, last)
# rna_model.create_constraints(0, 0, 0, 0, 0, 0, first, last)
# rna_model.create_objective(0, 0, 0, 0, 0, 1)

# print(len(rna_model.subsequences))
# for ssq in rna_model.subsequences:
#     print(ssq.var)


# bl_energies = []
# for bl in rna_model.bulge_loops:
#     bl_energies.append(bl.energy)

# ml_energies = []
# for ml in rna_model.multi_loops:
#     if ml.is_valid_size():
#         if ml.energy < 825:
#             ml_energies.append(ml.energy)


# il_energies = []
# for il in rna_model.internal_loops:
#     if il.energy < -85:
#         print(f'{il.bp1.i}_{il.bp1.j}_{il.bp2.i}_{il.bp2.j}')


# print(sorted(set(il_energies)))
# for nt in rna_model.nucleotides:
#     print(nt.VarName)
# rna_model.add_unpaired_nucleotides_constraints(first, last)