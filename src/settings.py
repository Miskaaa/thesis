"""
This file contains settings
"""

"""
Name of the reference otu -> ko profile file
"""
KO_PROFILE_FILENAME = "D:\DP\data\\test\\GG_ko_counts_new.tab"

"""
Name of the input biom table 
@TODO - make a parameter from command line?
"""
INPUT_FILENAME = "D:\DP\data\\test\sample1.tab"

"""
Name of the similarity score matrix
"""
SCORE_FILENAME = "D:\DP\data\distance_matrix_known.tab"

"""
Method to be used for dealing with otus with unknown ko profiles
Available methods:
 - alignment_simple
"""
UNKNOWN_METHOD = "phylogenetic_tree"

"""
How many similar otus to uknown otus we want to take into account
"""
SIMILAR = 10

"""
result filename
"""
RESULT_FILENAME = "D:\DP\data\\test\\result_ref80_sample1_test1.tab"

"""
File with 16S counts
"""
COUNT_16S_FILENAME = "D:\DP\data\\GG_16S_counts.tab"