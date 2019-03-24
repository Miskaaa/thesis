import random
import os

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def generate_new_reference_table(restriction, original_reference, new_reference):
    # settings
    random.seed()
    # count number of lines
    number_of_lines = file_len(original_reference) - 1

    # open reference file
    reference_file = open(original_reference, "r")

    # determine which lines will stay in the resulting table
    all_lines = list(range(number_of_lines))
    chosen_lines = random.sample(all_lines, int(number_of_lines - restriction*number_of_lines/100.0))

    # generate new reference file
    index = 0
    new_reference_file = open(new_reference, "w")

    line = reference_file.readline()
    new_reference_file.write(line)
    while index < number_of_lines:
        line = reference_file.readline()
        if index in chosen_lines:
            new_reference_file.write(line)
        index += 1
