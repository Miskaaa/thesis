"""
Deals with otus for which their ko profile is unknown
"""
import settings
import numpy as np
import common_functions
import random
from sklearn.externals import joblib

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

def get_similar_weighted(otus, known):
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

                diff = []
                value = 0
                for i in most_similiar:
                    j = int(i)
                    similar = head[j]
                    value = int(values[j])
                    if value > settings.MAX_DIFF and len(matrix[otu_name]) > 0:
                        break
                    if similar in known and similar not in matrix[otu_name]:
                        diff.append(values[j])
                        matrix[otu_name].append((similar,value))
                    if len(matrix[otu_name]) >= settings.SIMILAR:
                        break

                if value > settings.MAX_DIFF and len(matrix[otu_name]) > 0:
                    break
                increase_size += 10

    return matrix


def get_similar_treshold(otus, known):
    file = open(settings.SCORE_FILENAME, "r")
    line = file.readline()
    head = line.split("\t")[1:]
    matrix = {}

    FIRST_LIMIT = 10
    SEQ_LENGTH = 7682
    MAX_SIMILAR = 30

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

            sim_ratio = 1
            while sim_ratio >= settings.TRESHOLD:
                limit = FIRST_LIMIT+increase_size
                most_similiar = np.argpartition(values, limit)[:limit]

                diff = []
                for i in most_similiar:
                    j = int(i)
                    similar = head[j]
                    value = int(values[j])

                    sim_ratio = (SEQ_LENGTH - value) / SEQ_LENGTH

                    if sim_ratio < settings.TRESHOLD:
                        break

                    if similar in known and similar not in matrix[otu_name]:
                        diff.append(values[j])
                        matrix[otu_name].append(similar)
                if len(matrix[otu_name]) > MAX_SIMILAR:
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

                #print(most_similiar)

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
    tree_file = open(settings.TREE_FILE,"r")
    tree_file.readline()

    tree = []
    #indexes
    parent = 0
    desc1 = 1
    desc2 = 2

    while True:
        line = tree_file.readline()
        if not line:
            break
        line = line.replace(" ","").strip()
        tree.append(line.split(",")[:3])

    similar = {}
    target_number = settings.SIMILAR

    for otu in otus:
        look_for = [otu]
        searched = []
        similar[otu] = []
        while len(similar[otu]) < target_number:
            found = []
            for l in look_for:
                found += search_in_tree(tree, l)
                searched.append(l)
            look_for = []
            for f in found:
                #skontrolovat, ci sme uz nenasli, pridat pre dalsie hladanie
                if f[desc1] in known:
                    similar[otu].append(f[desc1])
                elif f[desc1] not in searched and f[desc1] not in look_for:
                    look_for.append(f[desc1])

                if f[desc2] in known:
                    similar[otu].append(f[desc2])
                elif f[desc2] not in searched and f[desc2] not in look_for:
                    look_for.append(f[desc2])

                if f[parent] not in searched and f[parent] not in look_for:
                    look_for.append(f[parent])

    return similar


def search_in_tree(tree, target):
    result = []
    for t in tree:
        if target in t:
            result.append(t)
    return result


def get_similar_regression(otus, result_dict, counts_16s):
    file = open("D:\DP\data\\seq_source.txt", "r")
    ko_list = result_dict.keys()

    models = {}
    for ko in ko_list:
        model = joblib.load("D:\DP\data\\models\\" + str(settings.RESTR) + "\\seq_" + ko.strip() + ".pkl")
        models[ko.strip()] = model

    combinations = []
    nucl = ["G", "A", "C", "T"]
    for i in nucl:
        for j in nucl:
            for k in nucl:
                combinations.append(i + j + k)

    while True:
        line = file.readline()
        if not line:
            break
        splitted = line.split("\t")
        otu_name = splitted[0]
        otu_seq = splitted[1]

        if otu_name in otus:
            counts = []
            for comb in combinations:
                counts.append(float(otu_seq.count(comb)))
            counts = [np.array(counts)]

            abundance = 1
            if otu_name in counts_16s:
                abundance = counts_16s[otu_name]

            for ko in ko_list:
                value = models[ko.strip()].predict(counts).tolist()[0]
                result_dict[ko]["even"] += value/abundance
                result_dict[ko]["staggered"] += value/abundance

    return result_dict


def combine(result_dict, ko_profile, otus, similarity, counts_16s):
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
                    tmp_result[ko_name]["even"] +=  even * ko_profile[otu][ko_name] / abundance
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
            result_dict[ko_name]["staggered"] += tmp_result[ko_name]["staggered"]

    return result_dict


def combine_weighted(result_dict, ko_profile, otus, similarity, counts_16s):
    for unknown in similarity:
        similar = similarity[unknown]
        tmp_result = common_functions.get_result_dict()
        even = float(otus[unknown]["even"])
        staggered = float(otus[unknown]["staggered"])

        # scitam vsetky podobne ko profily
        similar_with_profile = 0
        known_otu = None
        weight_sum = 0
        for sim_item in similar:
            otu = sim_item[0]
            if otu in ko_profile:
                multiplier = settings.MAX_DIFF + 1 - int(sim_item[1])
                if multiplier <= 0:
                    multiplier = 1
                weight_sum += multiplier
                known_otu = otu
                similar_with_profile += 1
                for ko_name in ko_profile[otu]:
                    abundance = 1
                    if otu in counts_16s:
                        abundance = counts_16s[otu]
                    tmp_result[ko_name]["even"] += multiplier * even * ko_profile[otu][ko_name] / abundance
                    tmp_result[ko_name]["staggered"] += multiplier * staggered * ko_profile[otu][ko_name] / abundance

        if not known_otu:
            continue

        # vydelim ich poctom podobnych
        for ko_name in ko_profile[known_otu]:
            tmp_result[ko_name]["even"] /= weight_sum
            tmp_result[ko_name]["staggered"] /= weight_sum

        # dam do celkoveho vysledku
        for ko_name in ko_profile[known_otu]:
            result_dict[ko_name]["even"] += tmp_result[ko_name]["even"]
            result_dict[ko_name]["staggered"] += tmp_result[ko_name]["staggered"]

    return result_dict


def predict_picrust(otus, result_dict, counts_16s):
    ko_profile = parse_ko_profile()

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

    return result_dict

def parse_ko_profile():
    file = open(settings.PICRUST_FILENAME, "r")
    line = file.readline()
    ko_ids = line.split("\t")[1:]
    ko_dict = {}

    while True:
        line = file.readline()
        if not line:
            break
        if line.isspace():
            continue

        line_parts = line.split("\t")
        id = line_parts[0]
        profile = {}

        for i in range(len(ko_ids)):
            ko_id = ko_ids[i]
            value = float(line_parts[i+1])
            profile[ko_id] = value

        ko_dict[id] = profile

    return ko_dict


def resolve_unknown(method, otus, known, counts_16s):
    ko_profile = common_functions.parse_ko_profile()
    result_dict = common_functions.get_result_dict()

    if method == "alignment_simple":
        similarity = get_similar(list(otus.keys()), known)
        result_dict = combine(result_dict, ko_profile, otus, similarity, counts_16s)
    elif method == "random":
        similarity = get_similar_random(list(otus.keys()), known)
        result_dict = combine(result_dict, ko_profile, otus, similarity, counts_16s)
    elif method == "tree":
        similarity = get_similar_phylogenetic_tree(list(otus.keys()), known)
        result_dict = combine(result_dict, ko_profile, otus, similarity, counts_16s)
    elif method == "treshold":
        similarity = get_similar_treshold(list(otus.keys()), known)
        result_dict = combine(result_dict, ko_profile, otus, similarity, counts_16s)
    elif method == "weighted":
        similarity = get_similar_weighted(list(otus.keys()), known)
        result_dict = combine_weighted(result_dict, ko_profile, otus, similarity, counts_16s)
    elif method == "regression_by_sequence":
        result_dict = get_similar_regression(list(otus.keys()), result_dict, counts_16s)
    elif method == "picrust":
        result_dict = predict_picrust(otus, result_dict, counts_16s)

    return result_dict
