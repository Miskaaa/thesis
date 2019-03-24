"""
Deals with otus for which their ko profile is unknown
"""
import settings
import numpy as np
import common_functions
import random


def get_similar(otus, known):
    file = open(settings.SCORE_FILENAME, "r")
    line = file.readline()
    head = line.split("\t")[1:]
    matrix = {}

    while True:
        line = file.readline()
        if not line:
            break
        splitted = line.split("\t")
        otu_name = splitted[0]

        if otu_name in otus:
            matrix[otu_name] = []
            values = splitted[1:-1]
            increase_size = 1
            while len(matrix[otu_name]) < settings.SIMILAR:
                limit = settings.SIMILAR+increase_size
                most_similiar = np.argpartition(values, limit)[:limit]

                print(most_similiar)

                diff = []
                for i in most_similiar:
                    j = int(i)
                    similar = head[j]
                    if similar in known and similar not in matrix[otu_name]:
                        diff.append(values[j])
                        matrix[otu_name].append(similar)
                    if len(matrix[otu_name]) >= settings.SIMILAR:
                        break
                increase_size += 10

    return matrix


def get_similar_random(otus, known):
    file = open(settings.SCORE_FILENAME, "r")
    line = file.readline()
    head = line.split("\t")[1:]
    matrix = {}

    while True:
        line = file.readline()
        if not line:
            break
        splitted = line.split("\t")
        otu_name = splitted[0]

        if otu_name in otus:
            matrix[otu_name] = []
            values = splitted[1:-1]
            increase_size = 1
            while len(matrix[otu_name]) < settings.SIMILAR:
                most_similiar = random.sample(values, settings.SIMILAR)

                print(most_similiar)

                diff = []
                for i in most_similiar:
                    j = int(i)
                    similar = head[j]
                    if similar in known and similar not in matrix[otu_name]:
                        diff.append(values[j])
                        matrix[otu_name].append(similar)
                    if len(matrix[otu_name]) >= settings.SIMILAR:
                        break
                increase_size += 10

    return matrix


def get_similar_phylogenetic_tree(otus, known):
    tree_file = open("D:\DP\data\\tree_normal.txt","r")
    tree_lines = tree_file.readlines()
    similar = {}

    for otu in otus:
        for line in tree_lines:
            slots = line.split()
            codes = slots[3].split("-")
            if otu in codes:
                codes.remove(otu)
                similar[otu] = codes
                break

    return similar


def resolve_unknown(method, otus, known, counts_16s):
    ko_profile = common_functions.parse_ko_profile()
    result_dict = common_functions.get_result_dict()

    similarity = []
    if method == "alignment_simple":
        similarity = get_similar(list(otus.keys()), known)
    elif method == "random":
        similarity = get_similar_random(list(otus.keys()), known)
    elif method == "phylogenetic_tree":
        similarity = get_similar_phylogenetic_tree(list(otus.keys()), known)

    for unknown in similarity:
        similar = similarity[unknown]
        tmp_result = common_functions.get_result_dict()
        even = float(otus[unknown]["even"])
        staggered = float(otus[unknown]["staggered"])

        # scitam vsetky podobne ko profily
        similar_with_profile = 0
        known_otu = None
        for otu in similar:
            if otu in ko_profile:
                known_otu = otu
                similar_with_profile += 1
                for ko_name in ko_profile[otu]:
                    abundance = 1
                    if otu in counts_16s:
                        abundance = counts_16s[otu]
                    tmp_result[ko_name]["even"] += even * ko_profile[otu][ko_name] / abundance
                    tmp_result[ko_name]["staggered"] += staggered * ko_profile[otu][ko_name] / abundance

        if not known_otu:
            continue

        # vydelim ich poctom podobnych
        for ko_name in ko_profile[known_otu]:
            tmp_result[ko_name]["even"] /= similar_with_profile
            tmp_result[ko_name]["staggered"] /= similar_with_profile

        # dam do celkoveho vysledku
        for ko_name in ko_profile[known_otu]:
            result_dict[ko_name]["even"] += tmp_result[ko_name]["even"]
            result_dict[ko_name]["staggered"] +=  tmp_result[ko_name]["staggered"]

    return result_dict
