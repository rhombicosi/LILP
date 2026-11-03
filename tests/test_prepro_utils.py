from pathlib import Path
import pytest
from LILP.utils.prepro_utils import ct_has_pairs


@pytest.mark.parametrize(
    ("content", "expected"),
    [
        pytest.param("", False, id="empty-file"),
        pytest.param("\n \n\t\n", False, id="blank-lines"),
        pytest.param("1 A 0 2 0 1\n", False, id="single-unpaired-base"),
        pytest.param("1 A 0 2 2 1\n", True, id="single-paired-base"),
        pytest.param("1 A 0 2 -2 1\n", False, id="negative-partner"),
        pytest.param("4 example sequence\n", False, id="header-only"),
        pytest.param(
            "2 example sequence\n1 A 0 2 0 1\n2 U 1 0 0 2\n",
            False,
            id="valid-ct-without-pairs",
        ),
        pytest.param(
            "2 example sequence\n1 A 0 2 2 1\n2 U 1 0 1 2\n",
            True,
            id="valid-ct-with-pair",
        ),
        pytest.param("through 511\n", False, id="through-header"),
        pytest.param("# 1 A 0 2 3 1\n", False, id="comment-line"),
        pytest.param("1 A 0 2\n", False, id="too-few-fields"),
        pytest.param("base A 0 2 3 1\n", False, id="non-digit-index"),
        pytest.param("1.0 A 0 2 3 1\n", False, id="decimal-index"),
        pytest.param("-1 A 0 2 3 1\n", False, id="negative-index"),
        pytest.param("+1 A 0 2 3 1\n", False, id="explicit-positive-index"),
        pytest.param("1 A 0 2 paired 1\n", False, id="non-integer-partner"),
        pytest.param("1 A 0 2 2.0 1\n", False, id="decimal-partner"),
        pytest.param(
            "1 A 0 2 invalid 1\n2 U 1 0 1 2\n",
            True,
            id="invalid-partner-before-valid-pair",
        ),
        pytest.param(
            "1 A 0 2 0 1\n2 U 1 3 0 2\n3 G 2 0 1 3\n",
            True,
            id="pair-on-later-line",
        ),
        pytest.param("1 A 0 2 3 1 ignored fields\n", True, id="extra-fields"),
        pytest.param("  1\tA  0\t2  3\t1  \n", True, id="mixed-whitespace"),
        pytest.param("1 A 0 2 0002 1\n", True, id="zero-padded-partner"),
        pytest.param("1 A 0 2 +2 1\n", True, id="explicit-positive-partner"),
        pytest.param("1 A 0 2 -0 1\n", False, id="negative-zero-partner"),
        pytest.param("001 A 0 2 2 1\n", True, id="zero-padded-index"),
        pytest.param("1 anything here is 7\n", True, id="only-required-columns-matter"),
        pytest.param("1 A 0 2 2e0 1\n", False, id="scientific-notation-partner"),
    ],
)
def test_ct_has_pairs_detects_positive_partner_indices(
    tmp_path: Path,
    content: str,
    expected: bool,
) -> None:
    ct_file = tmp_path / "structure.ct"
    ct_file.write_text(content, encoding="utf-8")

    assert ct_has_pairs(ct_file) is expected


def test_ct_has_pairs_propagates_missing_file_error(tmp_path: Path) -> None:
    missing_file = tmp_path / "missing.ct"

    with pytest.raises(FileNotFoundError):
        ct_has_pairs(missing_file)
