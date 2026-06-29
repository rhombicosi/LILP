from enum import Enum

VALID_PAIRS = ['AU','UA','CG','GC','GU','UG']
MIN_D = 3
MULTI_SIZE = 15
MULTI_MIN_D = 7
MULTI_BOUND = 5
MFE = -1000
SCALE = 100
M = 10000

BRANCH_START = 6 #FOR 80-90nts
# BRANCH_D = 6
# CBRANCH_D = 12

# R = 1.9872036 × 10-3	kcal.K-1.mol-1 is the gas constant and T is the absolute temperature, 310.15 K
RT = 0.616
Cbulge = -0.9

# multiloop params
# Mathews
# A = 8.4 # intitiation
# B = -0.8 # branches
# C = 0.0 # unpaired nucleotides
A = 2.0#9.3#3.4#10.1 # intitiation
B = 0.4#0.6#-0.6#0.4#-0.3 # branches
C = -0.3 # unpaired nucleotides
# Logarithmic
# A = 10.1 # intitiation
# B = -0.3 # branches
# C = -0.3 # unpaired nucleotides
# A = 3.4 # penalty for closing the multiloop
# B = 0.9 # penalty per branch (base pair)
# C = 0.4 # penalty per unpaired nucleotide
# A = 3.4 # penalty for closing the multiloop
# B = 0.91 # average symmetry (number of unpaired nucleotides)
# C = 0.63 # number of branching helices

class LoopType(Enum):
    HAIRPIN = "hairpin"
    INTERNAL = "internal"    
    STEM = "stem"
    BULGE = "bulge"
    MULTI = "multi"
    BRANCH = "branch"
    BRANCHPAIR = "branchpair"
    CBRANCH = "cbranch"

class InternalType(Enum):
    INT11 = "int11"
    INT12 = "int12"
    INT21 = "int21"
    INT1N = "int1n"
    INT22 = "int22"
    INT23 = "int23"
    INT32 = "int32"
    INTGEN = "general"

MAX_LOOP_SIZES = {
        LoopType.HAIRPIN: 10,
        LoopType.INTERNAL: 12,
        LoopType.BULGE: 5,
        LoopType.MULTI: 8,
        LoopType.BRANCH: 14,
        LoopType.BRANCHPAIR: 7,
        LoopType.CBRANCH: 5
    }

MIN_LOOP_SIZES = {
        LoopType.HAIRPIN: 3,
        LoopType.INTERNAL: 2,
        LoopType.BULGE: 1,
        LoopType.MULTI: 14,
        LoopType.BRANCH: 7,
        LoopType.BRANCHPAIR: 7,
        LoopType.CBRANCH: 14
    }

MAX_NUM_OF_LOOPS = {
        LoopType.HAIRPIN: 1, 
        LoopType.INTERNAL: 2,
        LoopType.BULGE: 5,
        LoopType.MULTI: 1,
        LoopType.BRANCH: 4
    }
