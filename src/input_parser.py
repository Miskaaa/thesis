"""
Input biom table parser
"""
import subprocess

def parse_input(filename):
    otu_dict = {}
    file = open(filename, "r")
    for line in file:
        if line.startswith("#"):
            continue
        splitted = line.strip().split()
        otu_dict[splitted[0]] = {"staggered": splitted[1], "even": splitted[2]}
    return otu_dict
