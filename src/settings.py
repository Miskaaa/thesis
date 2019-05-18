"""
This file contains settings
"""

"""
Name of the reference otu -> ko profile file
"""
KO_PROFILE_FILENAME = "../data/reference_data/ko_counts/GG_ko_counts_90.tab"

"""
Name of the input biom table 

"""
INPUT_FILENAME = "../data/test/sample1.tab"

"""
Name of the similarity score matrix
"""
SCORE_FILENAME = "../data/reference_data/distance_matrices/distance_matrix2.tab"

"""
Method to be used for dealing with otus with unknown ko profiles
Available methods:
 - alignment_simple, random, tree, treshold, weighted, regression_by_sequence
"""
UNKNOWN_METHOD = "alignment_simple"

"""
How many similar otus to uknown otus we want to take into account
"""
SIMILAR = 4

"""
result filename
"""
RESULT_FILENAME = "../data/test/result.tab"

"""
File with 16S counts
"""
COUNT_16S_FILENAME = "../data/reference_data/ko_counts/GG_16S_counts.tab"

"""
treshold for the treshold method
"""
TRESHOLD = 0.97

"""
max diff between sequences to be classified as smilar for weighted method
"""
MAX_DIFF = 200

"""
max diff between sequences to be classified as smilar for weighted method
"""
RESTR = "70"

"""
max diff between sequences to be classified as smilar for weighted method
"""
TREE_FILE = "../data/reference_data/phylogenetic_trees/upmga.csv"
