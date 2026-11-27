from pathlib import Path
import pandas as pd
from LILP.andronescu_parameters.parser import parse_vienna_params, flatten

# Converts the parsed ViennaRNA file into tables in the same format as
# tnn_parameters/ (same names, index and column labels), so loop classes
# can use them unchanged. Orientation of every table was checked against
# ViennaRNA evaluation with the same parameter file (tests/test_andronescu_tables.py).

BLSTAR_FILE = Path(__file__).parent / 'data' / 'parameters_BLstar_Vienna.txt'

# pair order of the rows/columns/blocks in the ViennaRNA file ('@' = non-standard pair)
VIENNA_PAIRS = ['CG', 'GC', 'GU', 'UG', 'AU', 'UA', '@']
VIENNA_NUCLEOTIDES = 'NACGU'
NUCLEOTIDES = ['A', 'C', 'G', 'U']
NUCLEOTIDE_PAIRS = [x + y for x in NUCLEOTIDES for y in NUCLEOTIDES]

# pair orders used by the tnn_parameters/ modules
STEM_PAIRS = ['AU', 'UA', 'CG', 'GC', 'GU', 'UG']
MISMATCH_PAIRS = ['AU', 'CG', 'GC', 'GU', 'UA', 'UG']
INTERNAL_PAIRS = ['AU', 'CG', 'GC', 'UA', 'GU', 'UG']
LOOPS = ['internal', 'bulge', 'hairpin']
MAX_LOOP_LENGTH = 30


def _nt(x: str) -> int:
    return VIENNA_NUCLEOTIDES.index(x)


def _block_row(pair: str, x: str) -> int:
    # sections with 5 rows per pair type (rows: N A C G U)
    return VIENNA_PAIRS.index(pair) * 5 + _nt(x)


def stack_table(sections: dict) -> pd.DataFrame:
    # LILP:   wcf_df.loc[first pair (i,j), second pair (i+1,j-1)], both read 5'->3'
    # Vienna: stack[type(i,j)][type(j-1,i+1)], the second pair is read backwards
    stack = sections['stack_energies']
    data = [[stack[VIENNA_PAIRS.index(p1)][VIENNA_PAIRS.index(p2[::-1])] for p2 in STEM_PAIRS]
            for p1 in STEM_PAIRS]
    return pd.DataFrame(data, index=STEM_PAIRS, columns=STEM_PAIRS)


def initiation_table(sections: dict) -> pd.DataFrame:
    # Vienna lists start at loop size 0, LILP tables at size 1
    columns = {loop: flatten(sections[name])[1:MAX_LOOP_LENGTH + 1]
               for loop, name in zip(LOOPS, ['internal_loop', 'bulge', 'hairpin'])}
    return pd.DataFrame(columns, index=range(1, MAX_LOOP_LENGTH + 1))[LOOPS]


def hairpin_mismatch_table(sections: dict) -> pd.DataFrame:
    # mismatch_df.loc[pair (i,j), nt i+1][nt j-1]
    rows = sections['mismatch_hairpin']
    index = pd.MultiIndex.from_product([MISMATCH_PAIRS, NUCLEOTIDES])
    data = [[rows[_block_row(pair, x)][_nt(y)] for y in NUCLEOTIDES] for pair, x in index]
    return pd.DataFrame(data, index=index, columns=NUCLEOTIDES)


def int11_table(sections: dict) -> pd.DataFrame:
    # int11_df.loc[(pair (i,j), nt i+1)][(pair (k,l), nt j-1)]
    blocks = sections['int11_energies']
    index = pd.MultiIndex.from_product([INTERNAL_PAIRS, NUCLEOTIDES])
    data = [[blocks[f'{p1}..{p2[::-1]}'][_nt(x)][_nt(y)] for p2, y in index] for p1, x in index]
    return pd.DataFrame(data, index=index, columns=index)


def int12_table(sections: dict) -> pd.DataFrame:
    # 1x2 loop: int12_df.loc[(pair (i,j), nt i+1)][(nt l+1, pair (k,l), nt j-1)]
    blocks = sections['int21_energies']
    index = pd.MultiIndex.from_product([INTERNAL_PAIRS, NUCLEOTIDES])
    columns = pd.MultiIndex.from_product([NUCLEOTIDES, INTERNAL_PAIRS, NUCLEOTIDES])
    data = [[blocks[f'{p1}.{x}..{p2[::-1]}'][_nt(n)][_nt(y)] for n, p2, y in columns] for p1, x in index]
    return pd.DataFrame(data, index=index, columns=columns)


def int22_table(sections: dict) -> pd.DataFrame:
    # 2x2 loop: int22_df.loc[(pair (i,j), nt i+1 + nt j-1)][(pair (k,l), nt k-1 + nt l+1)]
    blocks = sections['int22_energies']
    index = pd.MultiIndex.from_product([MISMATCH_PAIRS, NUCLEOTIDE_PAIRS])
    acgu = ''.join(NUCLEOTIDES)
    data = [[blocks[f'{p1}.{x[0]}{y[0]}..{p2[::-1]}'][acgu.index(y[1])][acgu.index(x[1])]
             for p2, y in index] for p1, x in index]
    return pd.DataFrame(data, index=index, columns=index)


def closure_penalty(sections: dict, pair: str) -> int:
    # In this parameter set the AU/GU closure penalty of internal loops is included in
    # mismatch_interior; the CC mismatch has no bonus, so it holds the bare penalty.
    return sections['mismatch_interior'][_block_row(pair, 'C')][_nt('C')]


def internal_mismatch_bonus_table(sections: dict) -> pd.DataFrame:
    # intnn_df.loc[pair][mismatch]: mismatch bonus only, the closure penalty is added
    # separately by InternalLoop (AU_end_penalty), so it is subtracted here
    rows = sections['mismatch_interior']
    data = [[rows[_block_row(pair, xy[0])][_nt(xy[1])] - closure_penalty(sections, pair)
             for xy in NUCLEOTIDE_PAIRS] for pair in MISMATCH_PAIRS]
    return pd.DataFrame(data, index=MISMATCH_PAIRS, columns=NUCLEOTIDE_PAIRS)


def terminal_au_penalty(sections: dict) -> int:
    # last value of '# ML_params': cu, cc, ci, TerminalAU
    return sections['ML_params'][0][3]


def asymmetry_penalty(sections: dict) -> int:
    # '# NINIO': penalty per nucleotide of asymmetry, maximum
    return sections['NINIO'][0][0]


blstar = parse_vienna_params(BLSTAR_FILE)
