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

class LILP:
    def __init__(self, rna_seq: str, name: str):
        self.rna_seq = rna_seq
        self.model = gp.Model(name)
        self.nucleotides : List[gp.Var] = []
        self.base_pairs : List[BasePair] = []
        self.first_pairs : List[BasePair] = []
        self.last_pairs : List[BasePair] = []
        self.hairpin_loops : List[HairpinLoop] = []
        self.stem_loops : List[StemLoop] = []
        self.internal_loops : List[InternalLoop] = []
        self.bulge_loops : List[BulgeLoop] = []
        self.multi_loops : List[MultiLoop] = []
        self.branches: List[Branch] = []
        self.branch_pairs: List[BranchPair] = []
        self.objective : gp.LinExpr   

    def create_nucleotides(self, first, last) -> None:
        for i in range(first, last + 1):
            var = self.model.addVar(vtype=GRB.BINARY, name=f'X_{i}')   
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
                    internal.add_variable(self.model)
                    self.internal_loops.append(internal)
        self.model.update()

    def create_bulge_loops(self) -> None:
        for bp1 in self.base_pairs:
            for bp2 in self.base_pairs:
                if (bp2.i == bp1.i + 1 and bp2.j < bp1.j - 1) or (bp2.i > bp1.i + 1 and bp2.j == bp1.j - 1):
                    bulge = BulgeLoop([bp1, bp2], self.rna_seq)
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
                    # if bp1.i >= MULTI_BOUND and bp2.j <= n - MULTI_BOUND:
                        branch = Branch([bp1, bp2], self.rna_seq)
                        branch.add_variable(self.model)
                        self.branches.append(branch)
        self.model.update()

    def create_branch_pairs(self):
        for bp in self.base_pairs:
            branch_pair = BranchPair(bp.i, bp.j, self.rna_seq)
            branch_pair.add_variable(self.model, 'BP')
            self.branch_pairs.append(branch_pair)
        self.model.update()

    def create_multiloop(self) -> None:
        bp1 = BasePair._find_base_pair_with_indices(self.base_pairs,5,49)[0]
        bp2 = BasePair._find_base_pair_with_indices(self.base_pairs,10,22)[0]
        bp3 = BasePair._find_base_pair_with_indices(self.base_pairs,24,40)[0]
        multi = MultiLoop([bp1,bp2,bp3], self.rna_seq)
        multi.add_variable(self.model)
        self.multi_loops.append(multi)
        self.model.update()

    def add_unpaired_nucleotides_constraints(self, first, last) -> None:
        for i in range(first, last + 1):
            inequality = gp.LinExpr(0)
            matches = BasePair._find_base_pairs_with_index(self.base_pairs, i)

            if matches:
                for bp in matches:
                    inequality.add(gp.LinExpr([1.0], [bp.var]))
                
                nucleotide = self.model.getVarByName(f'X_{i}')
                inequality.add(gp.LinExpr([1.0], [nucleotide]))
                self.model.addConstr(inequality == 1, f'UN-{i}')
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

    # def add_hairpin_ifthen_constraints(self) -> None:
    #     for hl in self.hairpin_loops:
    #         hl.create_hairpin_ifthen_constraint(self.model)
    #     self.model.update()

    # def add_hairpin_onlyif_constraints(self) -> None:
    #     for hl in self.hairpin_loops:
    #         hl.create_hairpin_onlyif_constraint(self.model, self.base_pairs)
    #     self.model.update()

    def add_hairpin_constraints(self) -> None:
        for hl in self.hairpin_loops: 
            if hl.energy > 0:
                hl.create_hairpin_ifthen_constraint(self.model)
            else:
                hl.create_hairpin_onlyif_constraint(self.model, self.base_pairs)

    def add_hairpin_max_number_constraint(self) -> None:
        HairpinLoop.create_hairpin_max_number_constraint(self.model, self.hairpin_loops)

    def add_internal_size_constraints(self) -> None:
        for il in self.internal_loops:
            il.create_internal_size_constraint(self.model)
        self.model.update()

    # def add_internal_ifthen_constraints(self) -> None:
    #     for il in self.internal_loops:
    #         il.create_internal_ifthen_constraint(self.model)
    #     self.model.update()

    # def add_internal_onlyif_constraints(self) -> None:
    #     for il in self.internal_loops:
    #         il.create_internal_onlyif_constraint(self.model, self.base_pairs)
    #     self.model.update()    
    
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

    # def add_bulge_ifthen_constraints(self) -> None:
    #     for bl in self.bulge_loops:
    #         bl.create_bulge_ifthen_constraint(self.model)
    #     self.model.update()

    # def add_bulge_onlyif_constraints(self) -> None:
    #     for bl in self.bulge_loops:
    #         bl.create_bulge_onlyif_constraint(self.model, self.base_pairs)
    #     self.model.update()

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

    # def add_multi_ifthen_constraints(self) -> None:
    #     for ml in self.multi_loops:
    #         ml.create_multi_ifthen_constraint(self.model)
    #     self.model.update()

    # def add_multi_onlyif_constraints(self) -> None:
    #     for ml in self.multi_loops:
    #         ml.create_multi_onlyif_constraint(self.model, self.base_pairs)
    #     self.model.update()
    
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
            b.create_branch_ifthen_constraint(self.model)
        self.model.update()    

    def add_branch_pairs_constraints(self) -> None:
        for bp in self.branch_pairs:
            bp_branches = Branch._find_branches_with_pair(self.branches, bp)
            if len(bp_branches) != 0:
                for bpb in bp_branches:
                    self.model.addConstr(bp.var >= bpb.var, f'BP-{bp.i}-{bp.j}-{bpb.bp1.i}-{bpb.bp1.j}-{bpb.bp2.i}-{bpb.bp2.j}')
                
                self.model.addConstr(bp.var <= gp.quicksum(bpb.var for bpb in bp_branches), f'BPS-{bp.i}-{bp.j}')
            else:
                self.model.addConstr(bp.var == 0)

    def create_stem_term(self) -> gp.LinExpr:
        objective = gp.LinExpr([sl.energy for sl in self.stem_loops], [sl.var for sl in self.stem_loops])
        return objective
    
    # def cut_stem_term(self, first, last) -> gp.LinExpr:
    #     stems = []
    #     for sl in self.stem_loops:
    #         if sl.first_pair.i >= first and sl.first_pair.i <= last and sl.first_pair.j >= first and sl.first_pair.j <= last and sl.last_pair.i >= first and sl.last_pair.i <= last and sl.last_pair.j >= first and sl.last_pair.j <= last:
    #             stems.append(sl)
    #     objective = gp.LinExpr([sl.energy for sl in stems], [sl.var for sl in stems])
    #     return objective       
    
    def create_hairpin_term(self) -> gp.LinExpr:
        objective = gp.LinExpr([hl.energy for hl in self.hairpin_loops], [hl.var for hl in self.hairpin_loops])
        return objective
    
    # def cut_hairpin_term(self, first, last) -> gp.LinExpr:
    #     hairpins = []
    #     for hl in self.hairpin_loops:
    #         if hl.base_pairs[0].i >= first and hl.base_pairs[0].i <= last and hl.base_pairs[0].j >= first and hl.base_pairs[0].j <= last:
    #             hairpins.append(hl)
    #     objective = gp.LinExpr([hl.energy for hl in hairpins], [hl.var for hl in hairpins])
    #     return objective
    
    def create_internal_term(self) -> gp.LinExpr:
        objective = gp.LinExpr([il.energy for il in self.internal_loops], [il.var for il in self.internal_loops])
        return objective
    
    # def cut_internal_term(self, first, last) -> gp.LinExpr:
    #     internals = []
    #     for il in self.internal_loops:
    #         if il.bp1.i >= first and il.bp1.i <= last and il.bp1.j >= first and il.bp1.j <= last and il.bp2.i >= first and il.bp2.i <= last and il.bp2.j >= first and il.bp2.j <= last:
    #             internals.append(il)
    #     objective = gp.LinExpr([il.energy for il in internals], [il.var for il in internals])
    #     return objective
    
    def create_bulge_term(self) -> gp.LinExpr:
        objective = gp.LinExpr([bl.energy for bl in self.bulge_loops], [bl.var for bl in self.bulge_loops])
        return objective
    
    # def cut_bulge_term(self, first, last) -> gp.LinExpr:
    #     bulges = []
    #     for bl in self.bulge_loops:
    #         if bl.base_pairs[0].i >= first and bl.base_pairs[0].i <= last and bl.base_pairs[0].j >= first and bl.base_pairs[0].j <= last and bl.base_pairs[1].i >= first and bl.base_pairs[1].i <= last and bl.base_pairs[1].j >= first and bl.base_pairs[1].j <= last:
    #             bulges.append(bl)
    #     objective = gp.LinExpr([bl.energy for bl in bulges], [bl.var for bl in bulges])
    #     return objective
    
    def create_multi_term(self) -> gp.LinExpr:
        objective = gp.LinExpr([ml.energy for ml in self.multi_loops], [ml.var for ml in self.multi_loops])
        return objective
    
    def create_branch_term(self) -> gp.LinExpr:
        objective = gp.LinExpr([bp.energy for bp in self.branch_pairs], [bp.var for bp in self.branch_pairs])
        return objective
    
    def create_cut(self, stem, hairpin, internal, bulge, multi, energy) -> gp.LinExpr:
        cut = gp.LinExpr()
        if stem:
            cut.add(self.create_stem_term())
        if hairpin:
            cut.add(self.create_hairpin_term())
        if internal:
            cut.add(self.create_internal_term())
        if bulge:
            cut.add(self.create_bulge_term())
        if multi:
            cut.add(self.create_multi_term())
        self.model.addConstr(cut <= energy, f'CUT')   

    ################ MODEL CONSCTRUCTION ###################################     
    
    def create_objective(self, stem, hairpin, internal, bulge, branch) -> gp.LinExpr:
        objective = gp.LinExpr()
        if stem:
            objective.add(self.create_stem_term())
        if hairpin:
            objective.add(self.create_hairpin_term())
        if internal:
            objective.add(self.create_internal_term())
        if bulge:
            objective.add(self.create_bulge_term())
        #if multi:
           #objective.add(self.create_multi_term())
        if branch:
            objective.add(self.create_branch_term())
        self.model.setObjective(objective, GRB.MINIMIZE)
        self.model.update()

    def create_variables(self, stem, hairpin, internal, bulge, branch, first, last):
        self.create_base_pairs(first, last)                       
        self.create_nucleotides(first, last)    
        if hairpin:
            self.create_hairpin_loops()       
        if stem:            
            self.create_stem_loops()          
        if branch:
            self.create_branches()
            self.create_branch_pairs()       
        if internal:  
            self.create_internal_loops()  
        if bulge:
            self.create_bulge_loops() 
        #if multi:
            #self.create_multi_loops()

        print(f"PAIRS: {len(self.base_pairs)}, NUCLEOTIDES: {len(self.nucleotides)}, HAIRPIN: {len(self.hairpin_loops)}, STEMPS: {len(self.stem_loops)}, BRANCHES: {len(self.branches)}, BRANCH PAIRS: {len(self.branch_pairs)} INTERNALS: {len(self.internal_loops)}, BULGES: {len(self.bulge_loops)}")    
    
    def create_constraints(self, stem, hairpin, internal, bulge, branch, first, last):
        self.add_single_pair_constraints(first, last)
        self.add_no_crossing_constraints()
        if stem:
            # self.add_stem_ifthen_constraints()
            # self.add_stem_onlyif_constraints()
            self.add_stem_constraints()
        self.add_unpaired_nucleotides_constraints(first, last)
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
            self.add_branch_pairs_constraints()
        #if multi:
            #self.add_multi_size_constraints()
            #self.add_multi_constraints()
            ## self.add_multi_ifthen_constraints()
            ## self.add_multi_onlyif_constraints()
            ## self.add_multi_max_number_constraint()
            ## self.add_multi_energy_constraints(825)

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

# rna_model.create_variables(0, 0, 0, 0, 1, first, last)
# rna_model.create_constraints(0, 0, 0, 0, 1, first, last)
# rna_model.create_objective(0, 0, 0, 0, 1)

# print(len(rna_model.branches))
# print(len(rna_model.branch_pairs))

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