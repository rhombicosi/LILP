from gurobipy import GRB
from collections import defaultdict

def make_callback():

    best_obj = float("inf")
    best_obj_time = None
    max_intinf = -1

    # store variables once, not repeatedly calling model.getVars()
    vars_container = {}

    def callback(model, where):

        nonlocal best_obj, best_obj_time, max_intinf

        # Store variables on first callback call
        if "vars" not in vars_container:
            vars_container["vars"] = model.getVars()

        vars_list = vars_container["vars"]

        # Track feasible solutions
        if where == GRB.Callback.MIPSOL:

            obj = model.cbGet(GRB.Callback.MIPSOL_OBJ)

            if obj < best_obj:
                best_obj = obj
                best_obj_time = model.cbGet(GRB.Callback.RUNTIME)


        # Inspect LP relaxation at nodes
        elif where == GRB.Callback.MIPNODE:

            status = model.cbGet(GRB.Callback.MIPNODE_STATUS)

            if status != GRB.OPTIMAL:
                return

            vals = model.cbGetNodeRel(vars_list)

            intinf = 0
            family_count = defaultdict(int)

            for v, val in zip(vars_list, vals):

                if v.VType not in (GRB.BINARY, GRB.INTEGER):
                    continue

                if abs(val - round(val)) > 1e-5:

                    intinf += 1

                    family = v.VarName.split("_")[0]
                    family_count[family] += 1


            if intinf <= max_intinf:
                return

            max_intinf = intinf

            node = int(model.cbGet(GRB.Callback.MIPNODE_NODCNT))
            runtime = model.cbGet(GRB.Callback.RUNTIME)

            print("\n" + "="*80)
            print("NEW MAXIMUM IntInf")
            print(f"Time   : {runtime:.2f}s")
            print(f"Node   : {node}")
            print(f"IntInf : {intinf}")

            print("\nFractional variables by family:")

            for fam, cnt in sorted(
                family_count.items(),
                key=lambda x: x[1],
                reverse=True
            ):
                print(f"{fam:20s} {cnt}")

            print("="*80)

    return callback, lambda: (best_obj, best_obj_time)