import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SOLVER = ROOT / 'test' / 'solve_example.py'


def solve(seq: str, params: str | None) -> dict:
    # fresh process: the parameter set is chosen at import time from LILP_PARAMS
    env = {**os.environ, 'PYTHONPATH': f'{ROOT}:{ROOT / "LILP"}'}
    env.pop('LILP_PARAMS', None)
    if params is not None:
        env['LILP_PARAMS'] = params
    out = subprocess.run([sys.executable, str(SOLVER), seq], env=env, cwd=ROOT,
                         capture_output=True, text=True, timeout=300)
    assert out.returncode == 0, out.stderr[-2000:]
    return json.loads(next(line for line in out.stdout.splitlines() if line.startswith('{')))


# Regression: results of the TNN model computed BEFORE the parameter switch was added
# (commit d28b96a + Andronescu tables). If any of these changes, the refactoring
# changed the TNN model.
TNN_BASELINE = [
    ('GGAGCAGAAAGCAUCC', -312, '(((((.....)).)))', {'BULGE': 171, 'HAIRPIN': 420, 'STEM': -903}),
    ('GGCUAGCAGAAAGCAAUGCC', -370, '(((..((.....))...)))', {'HAIRPIN': 420, 'INTERNAL': 220, 'STEM': -1010}),
    ('GGACUCGAUUCGCGAUAGUCC', -536, '(((((((.....)))..))))', {'BULGE': 280, 'HAIRPIN': 440, 'STEM': -1256}),
    ('GGAUCAGCUUCGGCAAGAUCC', -558, '(((((.((....))..)))))', {'HAIRPIN': 460, 'INTERNAL': 230, 'STEM': -1248}),
]


@pytest.mark.parametrize('params', [None, 'tnn'], ids=['default', 'explicit'])
@pytest.mark.parametrize(('seq', 'energy', 'structure', 'by_type'), TNN_BASELINE,
                         ids=[case[0] for case in TNN_BASELINE])
def test_tnn_regression(seq, energy, structure, by_type, params):
    result = solve(seq, params)
    assert result == {'energy': energy, 'structure': structure, 'energy_by_type': by_type}


def test_andronescu_matches_vienna_blstar():
    # RNAfold (dangles 0) with the same BL* file: MFE structure (((((.....)).))), -1.31 kcal/mol;
    # LILP differs by 1 in the bulge-1 energy (rounding/unit issue I8)
    result = solve('GGAGCAGAAAGCAUCC', 'andronescu')
    assert result['structure'] == '(((((.....)).)))'
    assert result['energy'] == -132


def test_loop_with_unpairable_nucleotide_is_paid():
    # issue I13: A5, A15, A16 cannot pair (no U); the 1x3 internal loop that contains
    # them used to be left out of the objective (energy -826 instead of -596)
    result = solve('GGCGAGCAGAAAGCAAGCGCC', 'tnn')
    assert result['energy_by_type'].get('INTERNAL') == 230
    assert result['energy'] == -596


def test_unknown_parameter_set_fails():
    env = {**os.environ, 'PYTHONPATH': f'{ROOT}:{ROOT / "LILP"}', 'LILP_PARAMS': 'turner'}
    out = subprocess.run([sys.executable, '-c', 'import LILP.parameters'], env=env, cwd=ROOT,
                         capture_output=True, text=True)
    assert out.returncode != 0 and 'LILP_PARAMS' in out.stderr
