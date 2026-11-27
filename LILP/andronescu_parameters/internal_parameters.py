import numpy as np
import pandas as pd
from LILP.andronescu_parameters.tables import (blstar, int11_table, int12_table, int22_table,
                                               internal_mismatch_bonus_table, closure_penalty,
                                               asymmetry_penalty)

# BL* internal loop parameters, same layout as tnn_parameters/internal_parameters.py
asymmetry = asymmetry_penalty(blstar)

# closure penalty for AU/GU pairs; in the parameter file it is part of mismatch_interior
AU_end_penalty = closure_penalty(blstar, 'AU')
GU_end_penalty = closure_penalty(blstar, 'GU')

int11_df = int11_table(blstar)
int12_df = int12_table(blstar)
int22_df = int22_table(blstar)

# mismatch bonuses (AG, GA, UU) for generic and 2x3 loops; the parameter file uses one
# table (mismatch_interior) for both
intnn_df = internal_mismatch_bonus_table(blstar)
int23_df = intnn_df
