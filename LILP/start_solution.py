from utils.constants_paths import *
from utils.prepro_run import *
from utils.prepro_utils import *
from utils.sol_converter import *
import forgi.graph.bulge_graph as cgb

def generate_start_sol(dotbracket_file, sol_start):

    with open(dotbracket_file, "r") as f:
        gen_brackets = f.read().strip()

    bg = cgb.BulgeGraph.from_dotbracket(gen_brackets)
    print("Dot-Bracket Notation:")
    print(bg.to_dotbracket_string())
    print("Bulge String:")
    print(bg.defines)

    loop_vars = []

    for loo, val in bg.defines.items():
        # print(loo)
        if loo.startswith('s'):
            for i, j in zip(range(val[0], val[1]), range(val[3], val[2], -1)):
                    loop_vars.append(f'STEM_{i}_{j}_{i + 1}_{j - 1}')
            for i, j in zip(range(val[0], val[1] + 1), range(val[3], val[2] - 1, -1)):
                    loop_vars.append(f'P_{i}_{j}')
        if loo.startswith('h'):
            loop_vars.append(f'HAIRPIN_{val[0] - 1}_{val[1] + 1}')
            loop_vars.append(f'P_{val[0] - 1}_{val[1] + 1}')
            for i in range(val[0], val[1] + 1):
                loop_vars.append(f'X_{i}')
        if loo.startswith('i'):
            left = bg.get_flanking_handles(loo, side=0)
            right = bg.get_flanking_handles(loo, side=1)
            # print(left)
            # print(right)
            if left[0] == left[1] - 1 or right[0] == right[1] - 1:
                loop_vars.append(f'BULGE_{left[0]}_{right[1]}_{left[1]}_{right[0]}')
                for i in range(left[0] + 1, left[1]):
                    loop_vars.append(f'X_{i}')
                for i in range(right[0] + 1, right[1]):
                    loop_vars.append(f'X_{i}')
                loop_vars.append(f'P_{left[0]}_{right[1]}')
                loop_vars.append(f'P_{left[1]}_{right[0]}')
            else:
                loop_vars.append(f'INTERNAL_{left[0]}_{right[1]}_{left[1]}_{right[0]}')
                for i in range(left[0] + 1, left[1]):
                    loop_vars.append(f'X_{i}')
                for i in range(right[0] + 1, right[1]):
                    loop_vars.append(f'X_{i}')
                loop_vars.append(f'P_{left[0]}_{right[1]}')
                loop_vars.append(f'P_{left[1]}_{right[0]}')

    juncs = bg.junctions[:-1]
    # print(juncs)

    for j in juncs:
        pairs = []
        for loo in j:
            # print(loo)
            handle = bg.get_flanking_handles(loo, side=0)
            # print(handle)
            pairs.append((handle[0], handle[1]))
            for i in range(handle[0] + 1, handle[1]):
                loop_vars.append(f'X_{i}')

        multi_indices = []
        multi_indices.append((pairs[0][0], pairs[-1][1]))
        loop_vars.append(f'P_{pairs[0][0]}_{pairs[-1][1]}')

        for ( _, end_curr), (start_next, _ ) in zip(pairs[:-1], pairs[1:]):
            multi_indices.append((end_curr, start_next))
            loop_vars.append(f'P_{end_curr}_{start_next}')

        # CBRANCH and BRANCH
        cb_start, cb_end = multi_indices[0], multi_indices[-1]
        i1, j1 = cb_start
        i2, j2 = cb_end

        cb_name = "CBRANCH_{}_{}_{}_{}".format(*multi_indices[0], *multi_indices[-1])
        loop_vars.append(cb_name)

        for p1, p2 in zip(multi_indices[1:-1], multi_indices[2:]):
            b_name = "BRANCH_{}_{}_{}_{}".format(*p1, *p2)
            loop_vars.append(b_name)

        for p in multi_indices[1:]:
            print(p)
            loop_vars.append(f'BP_{p[0]}_{p[1]}')

    # result = "_" + "_".join(str(n) for p in multi_indices for n in p)
    # loop_vars.append(f'MULTI{result}')

    # print(set(loop_vars))

    vars_to_one = set(loop_vars)
    # sol_path = "data\solution.sol"

    with open(sol_start, "w") as f:
        f.write("# Optimal\n")  # standard start for a CPLEX/MIP .sol file
        for var in sorted(vars_to_one):
            f.write(f"{var} 1\n")
        f.write("End\n")
