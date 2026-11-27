import math

import pytest

from warmstart import run
from warmstart.structure import base_pair_distance, parse, repair, start_values, to_dot_bracket

SEQ = 'GGGAAACCCAGGGAUUCCCU'          # 20 nt


def test_parse_nested_pairs() -> None:
    pairs, pseudoknots = parse('(((...)))..((....))')
    assert pairs == [(1, 9), (2, 8), (3, 7), (12, 19), (13, 18)]
    assert pseudoknots == 0


def test_parse_leaves_out_pseudoknot_brackets() -> None:
    pairs, pseudoknots = parse('((..[[..))..]]')
    assert pairs == [(1, 10), (2, 9)]
    assert pseudoknots == 2


@pytest.mark.parametrize('broken', ['((..)', '(..))', '(..x.)'])
def test_parse_rejects_broken_structure(broken: str) -> None:
    with pytest.raises(ValueError):
        parse(broken)


def test_dot_bracket_round_trip() -> None:
    db = '((((...))))...((...))'
    pairs, _ = parse(db)
    assert to_dot_bracket(pairs, len(db)) == db


def test_repair_removes_only_pairs_lilp_cannot_have() -> None:
    # (1,9) G-C ok, (2,8) G-C ok, (4,6) A-A non-canonical, (11,15) G-U closing 3 ok,
    # (12,14) G-U closing only 1 nucleotide
    pairs = [(1, 9), (2, 8), (4, 6), (11, 15), (12, 14)]
    seq = 'GGAAAACCCAGGAUUCCCUU'
    kept, removed = repair(pairs, seq)
    assert kept == [(1, 9), (2, 8), (11, 15)]
    assert removed == {'noncanonical': 1, 'short_loop': 1}


def test_repair_removes_pair_closing_less_than_three_nucleotides() -> None:
    seq = 'GAACAAAAAA'
    kept, removed = repair([(1, 4)], seq)          # G-C closing 2 nucleotides
    assert kept == [] and removed['short_loop'] == 1


def test_start_values_are_only_pair_variables() -> None:
    assert start_values([(1, 9), (2, 8)]) == {'P_1_9': 1, 'P_2_8': 1}


def test_base_pair_distance() -> None:
    assert base_pair_distance([(1, 9), (2, 8)], [(1, 9), (3, 7)]) == 2


def test_prepare_from_file_takes_prediction_time_from_budget(tmp_path, monkeypatch) -> None:
    tsv = tmp_path / 'pred.tsv'
    tsv.write_text('name\tstructure\tseconds\n'
                   f'seq1\t{"(((...)))..((....))."}\t3.25\n')
    monkeypatch.setenv('LILP_START_FILE', str(tsv))
    ws = run.prepare('file', SEQ, 'seq1', 'ws-mx-tnn', 900, str(tmp_path))
    assert ws['start_seconds'] == 3.25
    assert ws['time_limit'] == pytest.approx(900 - 3.25)
    assert ws['pairs_in'] == 5 and ws['start_pairs'] == 5
    assert ws['start_values'] == {f'P_{i}_{j}': 1 for i, j in [(1, 9), (2, 8), (3, 7), (12, 19), (13, 18)]}


def test_prepare_rejects_structure_of_other_length(tmp_path, monkeypatch) -> None:
    tsv = tmp_path / 'pred.tsv'
    tsv.write_text('name\tstructure\tseconds\nseq1\t(((...)))\t1\n')
    monkeypatch.setenv('LILP_START_FILE', str(tsv))
    with pytest.raises(ValueError):
        run.prepare('file', SEQ, 'seq1', 'ws-mx-tnn', 900, str(tmp_path))


def test_unknown_source_is_an_error(monkeypatch) -> None:
    monkeypatch.setenv('LILP_START', 'vienna')
    with pytest.raises(ValueError):
        run.source_name()


PROCESSED = ['Processing user MIP start: 0 nodes explored in subMIP, total elapsed time 4s',
             'User MIP start produced solution with objective -3589 (8.07s)',
             'User MIP start produced solution with objective -4587 (8.36s)',
             'Loaded user MIP start with objective -4587',
             'Processed MIP start in 8.49 seconds (3.65 work units)']


@pytest.mark.parametrize(('messages', 'first', 'expected'), [
    (PROCESSED, True, -3589.0),                    # the first completion of the partial start
    (PROCESSED, False, -4587.0),                   # the start the search begins with
    (['User MIP start produced solution with objective -5261 (5.60s)'], False, -5261.0),
    (['Loaded user MIP start with objective -120.5'], True, -120.5),
    (['User MIP start did not produce a new incumbent solution',
      'User MIP start violates constraint UNR_8 by 1.000000000'], False, math.nan),
    ([], True, math.nan),
])
def test_start_objective_from_gurobi_messages(messages, first, expected) -> None:
    value = run.start_objective(messages, first=first)
    assert (math.isnan(value) and math.isnan(expected)) or value == expected


def test_write_row_writes_header_once(tmp_path) -> None:
    path = tmp_path / 'WS.tsv'
    ws = {'name': 'seq1', 'length': 20, 'source': 'file', 'start_seconds': 3.2, 'time_limit': 896.8,
          'pairs_in': 5, 'pseudoknots': 0, 'removed_noncanonical': 0, 'removed_short_loop': 0,
          'start_pairs': 5, 'structure_in': '.', 'structure_start': '.'}
    info = {'start_vars_set': 5, 'start_messages': ['User MIP start produced solution with objective -10'],
            'n_sols': 2, 'status': 9, 'runtime': 896.8, 'final_obj': -12.0}
    run.write_row(str(path), ws, info)
    run.write_row(str(path), ws, info)
    lines = path.read_text().splitlines()
    assert len(lines) == 3 and lines[0].split('\t') == run.COLUMNS
    row = dict(zip(run.COLUMNS, lines[1].split('\t')))
    assert row['start_obj'] == '-10.0' and row['start_obj_first'] == '-10.0'
    assert row['n_sols'] == '2' and row['source'] == 'file'
