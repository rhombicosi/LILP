import pandas as pd
import numpy as np

c = 100

# intermolecular initiation parameter
wcf_initiation = 4.09 * c

# penalty per each AU(5'3')/UA(3'5') end
wcf_AU_end_penalty = 0.45 * c

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

wcf_pairs = ['AU','UA','CG','GC', 'GU', 'UG']


# NOTE : Although the free energy change parameter for GU followed by GU fit to +0.47 kcal/mol,
# there is a large error and this parameter is set to –0.5 kcal/mol for 
# secondary structure prediction to optimize accuracy.
energies = np.array([[-0.69, -0.84, -1.29, -1.39, -0.14, -0.81],
                     [-0.68, -0.69, -1.23, -1.32, -0.02, -0.58],
                     [-1.32, -1.39, -2.07, -1.33, -0.37, -1.46],
                     [-1.23, -1.29, -2.05, -2.07, -0.99, -1.50],
                     [-0.58, -0.81, -1.50, -1.46, -0.67, -0.22],
                     [-0.02, -0.14, -0.91, -0.37, -0.38, -0.67]])

# TODO: NOTE: GU followed by UG is generally unfavorable (+0.47 kcal/mol), 
# but is favorable in the one context shown (GC preceding pair and CG following pair).
# For the favorable case, the single reported parameter (–4.12 kcal/mol) is used 
# for the total of three stacks.

data = energies * c
rows = wcf_pairs
cols = wcf_pairs

wcf_df = pd.DataFrame(
    data=data,
    index=rows,
    columns=cols
)

print('\nDataFrame:', wcf_df)

# print(wcf_df.loc['AU','AU'])