"""Warm start of the LILP model.

Any source of a secondary structure (the simplified LILP model, a file written by any
other program, e.g. MXfold2) gives a structure in dot-bracket notation; the structure is
repaired (pairs that LILP cannot have at all are removed) and handed to Gurobi as a
partial MIP start: only the base-pair variables P_i_j = 1, Gurobi completes the rest
(loop variables, and pairs it has to add, e.g. to split a loop longer than LILP allows).

Common time budget: the cost of producing the start (prediction, or building and
solving the simplified model) is taken from the time limit of the full model,
    TimeLimit(full) = LILP_TIME_LIMIT - cost of the start.
Building the full model is not counted, the same as without a start.

Environment (set by run_bench_ws.sh):
    LILP_START       none | lilp_init | file
    LILP_START_FILE  TSV with columns name, structure, seconds (source 'file')
    LILP_T_INIT      time limit of the simplified model in seconds (source 'lilp_init')
"""
