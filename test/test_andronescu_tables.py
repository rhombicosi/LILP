import pytest
from LILP.andronescu_parameters.stem_parameters import wcf_df, wcf_AU_end_penalty
from LILP.andronescu_parameters.initiation_parameters import initiation_df
from LILP.andronescu_parameters.terminal_mismatch import mismatch_df
from LILP.andronescu_parameters.internal_parameters import (int11_df, int12_df, int22_df, intnn_df,
                                                            AU_end_penalty, asymmetry)
from LILP.lilp_config import M

# Expected loop energies below were computed independently with ViennaRNA 2.7.2
# (eval_int_loop / eval_hp_loop) loaded with the same BL* file, converted by RNAparconv.


@pytest.mark.parametrize(
    ("first", "second", "expected"),
    [
        pytest.param("AU", "AU", -69, id="5'AA/3'UU"),
        pytest.param("AU", "UA", -84, id="5'AU/3'UA"),
        pytest.param("GC", "GC", -207, id="5'GG/3'CC"),
        pytest.param("GC", "CG", -205, id="5'GC/3'CG"),
        pytest.param("CG", "GC", -133, id="5'CG/3'GC"),
        pytest.param("GU", "UG", -22, id="5'GU/3'UG"),
        pytest.param("GC", "GU", -91, id="5'GG/3'CU"),
    ],
)
def test_stack_orientation_matches_blstar_file(first: str, second: str, expected: int) -> None:
    assert wcf_df.loc[first, second] == expected


def test_stack_table_is_symmetric() -> None:
    # the same stack read from the other strand: 5'AC/3'BD == 5'DB/3'CA
    for first in wcf_df.index:
        for second in wcf_df.columns:
            assert wcf_df.loc[first, second] == wcf_df.loc[second[::-1], first[::-1]]


@pytest.mark.parametrize(
    ("size", "loop", "expected"),
    [
        pytest.param(3, "hairpin", 364, id="hairpin-3"),
        pytest.param(9, "hairpin", 272, id="hairpin-9"),
        pytest.param(10, "hairpin", 283, id="hairpin-10-extrapolated"),
        pytest.param(1, "bulge", 281, id="bulge-1"),
        pytest.param(6, "bulge", 272, id="bulge-6"),
        pytest.param(4, "internal", 83, id="internal-4"),
        pytest.param(30, "hairpin", 401, id="hairpin-30"),
    ],
)
def test_initiation_values(size: int, loop: str, expected: int) -> None:
    assert initiation_df.loc[size, loop] == expected


def test_forbidden_loop_sizes_are_m() -> None:
    assert initiation_df.loc[2, "hairpin"] == M
    assert initiation_df.loc[3, "internal"] == M


def test_hairpin_matches_vienna() -> None:
    # GCGAAACGC, hairpin closed by (2,8) = C-G, mismatch G..C, size 5
    assert initiation_df.loc[5, "hairpin"] + mismatch_df.loc[("CG", "G"), "C"] == 273


def test_int11_matches_vienna() -> None:
    # GACAAAAGUC: (1,10) G-C, (3,8) C-G, unpaired A2 and U9
    assert int11_df.loc[("GC", "A"), ("CG", "U")] == 35


def test_int12_matches_vienna() -> None:
    # CUGAAAACAAG: (1,11) C-G, (3,8) G-C, unpaired U2 and A9 A10
    assert int12_df.loc[("CG", "U"), ("A", "GC", "A")] == 171


def test_int22_matches_vienna() -> None:
    # AGAUAAAAAGCU: (1,12) A-U, (4,9) U-A, unpaired G2 A3 and G10 C11
    assert int22_df.loc[("AU", "GC"), ("UA", "AG")] == 242


def test_internal_loop_penalties() -> None:
    assert AU_end_penalty == 63
    assert asymmetry == 50
    assert wcf_AU_end_penalty == 56
    assert intnn_df.loc["AU", "AG"] == -50
    assert intnn_df.loc["CG", "UU"] == -45
    assert intnn_df.loc["GC", "CC"] == 0
