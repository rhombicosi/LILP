import pandas as pd
import numpy as np

scale = 100

# intermolecular initiation parameter
wcf_initiation = 4.09 * scale

# penalty per each AU(5'3')/UA(3'5') end
wcf_AU_end_penalty = 0.45 * scale

# Although the free energy change parameter 
# for GU followed by GU fit to +0.47 kcal/mol,
# there is a large error and this parameter 
# is set to –0.5 kcal/mol for secondary structure prediction to optimize accuracy.

# GU followed by UG is generally unfavorable (+0.47 kcal/mol),
# but is favorable in the one context shown (GC preceding pair and CG following pair).
# For the favorable case, the single reported parameter (–4.12 kcal/mol)
# is used for the total of three stacks.

# 5' GGUC 3'
# 3' CUGG 5'

# 5' GGUC 3'
# 3' CUGG 5'

wcf_pairs = ['AU', 'UA', 'CG', 'GC', 'GU', 'UG']


# NOTE : Although the free energy change parameter for GU followed by GU fit to +0.47 kcal/mol,
# there is a large error and this parameter is set to –0.5 kcal/mol for 
# secondary structure prediction to optimize accuracy.
energies = np.array([[-0.93, -1.10, -2.24, -2.08, -0.55, -1.36],
                     [-1.33, -0.93, -2.35, -2.11, -1.00, -1.27],
                     [-2.11, -2.08, -3.26, -2.36, -1.41, -2.11],
                     [-2.35, -2.24, -3.42, -3.26, -1.53, -2.51],
                     [-1.27, -1.36, -2.51, -2.11, -0.50, 1.29],
                     [-1.00, -0.55, -1.53, -1.41, 0.30, -0.50]])

# TODO: NOTE: GU followed by UG is generally unfavorable (+0.47 kcal/mol), 
# but is favorable in the one context shown (GC preceding pair and CG following pair).
# For the favorable case, the single reported parameter (–4.12 kcal/mol) is used 
# for the total of three stacks.

data = energies * scale
rows = wcf_pairs
cols = wcf_pairs

wcf_df = pd.DataFrame(
    data=data,
    index=rows,
    columns=cols
)

# print('\nDataFrame:', wcf_df)

# print(wcf_df.loc['AU','AU'])