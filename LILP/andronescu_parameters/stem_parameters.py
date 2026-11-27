import numpy as np
import pandas as pd
from LILP.andronescu_parameters.tables import blstar, stack_table, terminal_au_penalty

# BL* stacking energies (Andronescu et al. 2010), in 0.01 kcal/mol, same layout as tnn_parameters
wcf_df = stack_table(blstar)

# terminal AU/GU penalty of the BL* set ('TerminalAU' in ML_params)
wcf_AU_end_penalty = terminal_au_penalty(blstar)
