import re
from LILP.lilp_config import M

# Parser for RNAfold 1.8 parameter files (e.g. parameters_BLstar_Vienna.txt
# from Andronescu et al. 2010). Values in the file are in 0.01 kcal/mol,
# the same scale as LILP energies, so they are NOT multiplied by 100.
#
# Result: dict  section name -> data, where data is
#   - list of rows (list of values)              e.g. 'stack_energies', 'hairpin'
#   - dict  block label -> list of rows          e.g. 'int11_energies' -> 'CG..GC'
#   - dict  loop sequence -> value               'Tetraloops', 'Triloops'

SECTION_RE = re.compile(r'^#\s*(\S+)')
# block labels: 'CG..GC' (int11), 'CG.A..CG' (int21), 'CG.AA..CG' (int22),
# '@' marks a non-standard pair and is padded with a space, e.g. '/*  @..AU */'
BLOCK_LABEL_RE = re.compile(r'^/\*\s*([ACGU@ .]+\.\.[ACGU@ ]+?)\s*\*/$')
COMMENT_RE = re.compile(r'/\*.*?\*/')
LOOP_SEQUENCE_SECTIONS = ('Tetraloops', 'Triloops')


def parse_value(token: str) -> int:
    # 'INF' (element not allowed) and 'NST' (pair '@' that cannot stack) both
    # become M, the value LILP already uses for forbidden elements
    if token in ('INF', 'NST'):
        return M
    return int(token)


def parse_vienna_params(filepath: str) -> dict:
    sections = {}
    name = None
    block = None

    with open(filepath) as fp:
        for line in fp:
            line = line.strip()
            if not line or line.startswith('##'):
                continue

            if line.startswith('#'):
                name = SECTION_RE.match(line).group(1)
                if name == 'END':
                    break
                sections[name] = {} if name in LOOP_SEQUENCE_SECTIONS else []
                block = None
                continue

            if name is None:
                continue

            label = BLOCK_LABEL_RE.match(line)
            if label:
                if not isinstance(sections[name], dict):
                    sections[name] = {}
                block = label.group(1).replace(' ', '')
                sections[name][block] = []
                continue

            line = COMMENT_RE.sub('', line).strip()
            if not line:
                continue

            tokens = line.split()
            if name in LOOP_SEQUENCE_SECTIONS:
                sections[name][tokens[0]] = parse_value(tokens[1])
            elif block is not None:
                sections[name][block].append([parse_value(t) for t in tokens])
            else:
                sections[name].append([parse_value(t) for t in tokens])

    return sections


def flatten(rows: list) -> list:
    return [value for row in rows for value in row]
