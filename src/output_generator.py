"""
Will get otu profiles, combine them ang generate output file
"""

import settings

def make_output(known_otus, unknown_otus):
    expected_result = open(settings.RESULT_FILENAME, "w")
    expected_result.write("#KO ID\tstaggered\teven\n")
    index = 0
    for ko_name in known_otus:
        to_write = ko_name.strip() + "\t"
        to_write += "{0:.1f}".format(known_otus[ko_name]["staggered"] + unknown_otus[ko_name]["staggered"]) + "\t"
        to_write += "{0:.1f}".format(known_otus[ko_name]["even"] + unknown_otus[ko_name]["even"]) + "\n"
        expected_result.write(to_write)
        index += 1
    expected_result.close()
