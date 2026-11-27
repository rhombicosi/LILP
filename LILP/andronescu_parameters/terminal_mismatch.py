import numpy as np
import pandas as pd
from LILP.andronescu_parameters.tables import blstar, hairpin_mismatch_table

# BL* terminal mismatches in hairpin loops: mismatch_df.loc[pair, nt i+1][nt j-1]
mismatch_df = hairpin_mismatch_table(blstar)
