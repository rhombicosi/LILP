import pandas as pd
import numpy as np

# Terminal mismatches are non-canonical pairs adjacent to helix ends.
c = 100

bps = np.array(['AU','CG','GC','GU','UA','UG'])

nucleotides = np.array(['A', 'C', 'G', 'U'])

mismatch = np.array([[[0.41, 0.77, 0.65, 1.17], [-0.02, 0.08, 0.44, 0.43], [-0.31, 0.23, 0.15, 0.60], [0.29, -0.06, 0.32, -0.10]],
                     [[-0.35, -0.13, 0.03, 0.13], [-0.12, -0.17, -0.98, -0.22], [-0.51, -0.24, -0.01, 0.03], [-0.16, -0.57, -0.97, -0.38]],
                     [[0.01, -0.22, -0.65, 1.10], [0.11, -0.27, -0.17, 0.00], [-0.67, 0.01, -0.45, 0.36], [0.08, -0.39, -0.55, -0.92]],
                     [[0.60, 0.75, 0.43, 1.80], [0.50, 0.44, 0.70, -0.01], [-0.92, 0.15, 0.39, 1.14], [0.25, 0.07, 0.67, 0.36]],
                     [[0.17, 0.52, 0.76, 0.68], [0.53, 0.36, 0.37, 0.11], [-0.39, 0.46, 0.16, 0.68], [0.50, -0.16, 0.46, 0.29]],
                     [[0.43, 0.52, 0.87, 1.09], [0.52, 0.35, -0.10, 0.40], [-0.26, -0.28, 0.05, 0.92], [0.42, -0.38, -0.14, -0.30]]])


mismatch_data = mismatch.reshape(bps.size*nucleotides.size, nucleotides.size) * c

# print(nucleotides)
# print(bps)
# print(mismatch_data)

midx = pd.MultiIndex.from_product([bps, nucleotides])

mismatch_df = pd.DataFrame(mismatch_data, index = midx, columns=nucleotides) 

# print(mismatch_df)

# print(mismatch_df.loc['CG', 'G']['A'])
