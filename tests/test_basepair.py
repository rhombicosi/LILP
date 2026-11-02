import pytest
from LILP.basepair import BasePair


@pytest.mark.parametrize(
    ("rna", "i", "j", "expected"),
    [
        pytest.param("AAAAU", 1, 5, True, id="pair-AU"),
        pytest.param("UAAAA", 1, 5, True, id="pair-UA"),
        pytest.param("CAAAG", 1, 5, True, id="pair-CG"),
        pytest.param("GAAAC", 1, 5, True, id="pair-GC"),
        pytest.param("GAAAU", 1, 5, True, id="pair-GU"),
        pytest.param("UAAAG", 1, 5, True, id="pair-UG"),
        pytest.param("AAAAA", 1, 5, False, id="non-canonical-AA"),
        pytest.param("AAAAG", 1, 5, False, id="non-canonical-AG"),
        pytest.param("CAAAU", 1, 5, False, id="non-canonical-CU"),
        pytest.param("GAAC", 1, 4, False, id="loop-of-2-too-short"),
        pytest.param("GAAAC", 1, 5, True, id="loop-of-3-shortest-allowed"),
        pytest.param("AGAAACA", 2, 6, True, id="indices-are-1-based"),
    ],
)
def test_is_valid_checks_nucleotides_and_loop_size(
    rna: str,
    i: int,
    j: int,
    expected: bool,
) -> None:
    assert BasePair(i, j, rna).is_valid() is expected
