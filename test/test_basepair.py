import pytest
import LILP.basepair as bp_module
from LILP.basepair import BasePair
from LILP.utils.prepro_run import *
from LILP.tnn_parameters.stem_parameters import *
import gurobipy as gp
from gurobipy import GRB

# sequence used in test is srp_Acin.spec._CR543861.txt 47nts

seq_number = 0

chain_file = seq_files[seq_number]
chain_name_with_ext = os.path.basename(chain_file)        
chain_name_without_ext = os.path.splitext(chain_name_with_ext)[0]
lp_file_name = chain_name_without_ext

seq_data = parse_seq_file(chain_file)
RNA_seq = seq_data['sequence']  

# TEST 1: Initialization
def test_basepair_initialization():
    bp = BasePair(1, 5, RNA_seq)

    assert bp.i == 1
    assert bp.j == 5
    assert bp.rna == RNA_seq
    assert bp.nt1 == "G"
    assert bp.nt2 == "C"
    assert bp.var is None
    assert bp.distance == 3
    assert bp.pair_penalty_energy == 0


# TEST 2: Distance
@pytest.mark.parametrize(
    "i, j, expected",
    [
        (1, 2, 0),
        (1, 10, 8),
        (5, 20, 14),
        (10, 30, 19),
        (20, 47, 26),
    ]
)
def test_calculate_distance(i, j, expected):
    bp = BasePair(i, j, RNA_seq)

    assert bp.calculate_distance() == expected

# TEST 3: Penalty
@pytest.mark.parametrize(
    "i, j, expected",
    [
        (1, 5, 0),  # GC
        (2, 6, 0),  # CG
        (7, 37, round(wcf_AU_end_penalty)),  # AU
        (37, 46, round(wcf_AU_end_penalty)), # UA
        (1, 37, round(wcf_AU_end_penalty)),  # GU
        (37, 41, round(wcf_AU_end_penalty)), # UG
    ]
)
def test_calculate_pair_penalty_energy(i, j, expected):
    bp = BasePair(i, j, RNA_seq)

    assert bp.distance >= 3
    assert bp.calculate_pair_penalty_energy() == expected

# TEST 3: Add variable for valid base pair
def test_add_variable_valid():
    bp = BasePair(1, 5, RNA_seq)

    model = gp.Model("test_model")
    model.Params.OutputFlag = 0

    var = bp.add_variable(model, "BP")
    model.update()

    assert var is not None
    assert var.VType == GRB.BINARY
    assert var.VarName == "BP_1_5"
    assert bp.var is var

    model.dispose()

# TEST 4: Do NOT add variable for invalid base pair
def test_add_variable_invalid_pair(monkeypatch):
    monkeypatch.setattr(
        bp_module, "VALID_PAIRS",
        {"AU", "UA", "GC", "CG", "GU", "UG"}
    )
    monkeypatch.setattr(bp_module, "MIN_D", 3)

    bp = BasePair(12, 19, RNA_seq)

    model = gp.Model("test_model")
    model.Params.OutputFlag = 0

    var = bp.add_variable(model, "BP")

    assert var is None
    assert bp.var is None

    model.dispose()

# TEST 5: Do NOT add variable for pair of nts that are separated by less than 3 nts
def test_add_variable_invalid_distance(monkeypatch):
    monkeypatch.setattr(
        bp_module, "VALID_PAIRS",
        {"AU", "UA", "GC", "CG", "GU", "UG"}
    )
    monkeypatch.setattr(bp_module, "MIN_D", 3)

    bp = BasePair(16, 19, RNA_seq)

    model = gp.Model("test_model")
    model.Params.OutputFlag = 0

    var = bp.add_variable(model, "BP")

    assert var is None
    assert bp.var is None

    model.dispose()



