import numpy as np
import pandas as pd
from LILP.andronescu_parameters.tables import blstar, initiation_table

# BL* loop initiation for sizes 1..30; sizes above 9 (hairpin) / 6 (bulge, internal)
# are extrapolated in the parameter file; forbidden sizes are M
initiation_df = initiation_table(blstar)
