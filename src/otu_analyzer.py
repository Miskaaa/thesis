"""
Analyzes otus from input
"""
import common_functions
import settings


def analyze_otus(otus, counts_16s):
    ko_profile = common_functions.parse_ko_profile()
    result_dict = common_functions.get_result_dict()
    unknown_list = {}

    known_otus = list(ko_profile.keys())

    for otu in otus:
        if otu in ko_profile:
            even = float(otus[otu]["even"])
            staggered = float(otus[otu]["staggered"])
            abundance = 1
            if otu in counts_16s:
                abundance = counts_16s[otu]

            for ko_name in ko_profile[otu]:
                result_dict[ko_name]["even"] += even * ko_profile[otu][ko_name] / abundance
                result_dict[ko_name]["staggered"] += staggered * ko_profile[otu][ko_name] / abundance
        else:
            unknown_list[otu] = otus[otu]

    return unknown_list, known_otus, result_dict
