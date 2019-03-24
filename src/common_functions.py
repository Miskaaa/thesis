import settings


def parse_ko_profile():
    file = open(settings.KO_PROFILE_FILENAME, "r")
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


def get_result_dict():
    file = open(settings.KO_PROFILE_FILENAME, "r")
    line = file.readline()
    ko_ids = line.split("\t")[1:]
    result_dict = {}

    for ko_id in ko_ids:
        result_dict[ko_id] = {"even": 0.0, "staggered": 0.0}
    file.close()

    return result_dict


def get_16s_counts():
    file = open(settings.COUNT_16S_FILENAME, "r")
    line = file.readline()
    count_dict = {}

    while True:
        line = file.readline()
        if not line:
            break
        if line.isspace():
            continue

        line_parts = line.strip().split("\t")
        id = line_parts[0]
        count = int(line_parts[1])

        count_dict[id] = count

    return count_dict
