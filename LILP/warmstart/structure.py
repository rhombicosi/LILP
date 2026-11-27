"""Secondary structure <-> base pairs, and repair of a structure before it becomes a MIP start."""
from lilp_config import MIN_D, VALID_PAIRS

OPEN, CLOSE = '(', ')'
PSEUDOKNOT_CLOSE = ']}>'        # other bracket types = pseudoknot pairs, not representable in LILP
UNPAIRED = '.,:-_'


def parse(dot_bracket):
    """Dot-bracket -> (nested pairs (i, j) 1-based with i < j, number of pseudoknot pairs left out)."""
    stack, pairs, pseudoknots = [], [], 0
    for pos, c in enumerate(dot_bracket, 1):
        if c == OPEN:
            stack.append(pos)
        elif c == CLOSE:
            if not stack:
                raise ValueError(f"unbalanced ')' at position {pos}")
            pairs.append((stack.pop(), pos))
        elif c in PSEUDOKNOT_CLOSE:
            pseudoknots += 1
        elif c not in UNPAIRED and c not in '[{<':
            raise ValueError(f'unknown character {c!r} at position {pos}')
    if stack:
        raise ValueError(f"unbalanced '(' at position {stack[-1]}")
    return sorted(pairs), pseudoknots


def to_dot_bracket(pairs, length):
    fold = ['.'] * length
    for i, j in pairs:
        fold[i - 1], fold[j - 1] = OPEN, CLOSE
    return ''.join(fold)


def repair(pairs, seq):
    """Remove the pairs for which the LILP model has no variable at all.

    A pair exists in LILP only if it is canonical or G-U (VALID_PAIRS) and closes at
    least MIN_D nucleotides (BasePair.is_valid). Loops longer than LILP allows are NOT
    repaired here: removing a pair only makes a loop longer; Gurobi may add pairs when
    it completes the partial start.
    Returns (kept pairs, {reason: number of removed pairs}).
    """
    kept, removed = [], {'noncanonical': 0, 'short_loop': 0}
    for i, j in pairs:
        if seq[i - 1] + seq[j - 1] not in VALID_PAIRS:
            removed['noncanonical'] += 1
        elif j - i - 1 < MIN_D:
            removed['short_loop'] += 1
        else:
            kept.append((i, j))
    return kept, removed


def start_values(pairs):
    """Partial MIP start: P_i_j = 1 for every pair of the structure, nothing else."""
    return {f'P_{i}_{j}': 1 for i, j in pairs}


def base_pair_distance(pairs_a, pairs_b):
    """Number of pairs in exactly one of the two structures."""
    return len(set(pairs_a) ^ set(pairs_b))
