from LILP.andronescu_parameters.tables import blstar, terminal_au_penalty

# The BL* parameter file has no separate special-hairpin terms (GG first mismatch,
# special GU closure, poly-C loops); they are set to 0 so HairpinLoop adds nothing.
hp_mismatch = {'GA': 0, 'GG': 0, 'UU': 0}
spec_GU_clos = 0
c3_loop = 0
all_C_loop_A = 0
all_C_loop_B = 0
AU_end_penalty = terminal_au_penalty(blstar)
GU_end_penalty = terminal_au_penalty(blstar)
