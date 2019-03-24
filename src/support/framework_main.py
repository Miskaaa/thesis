from restrict_reference import generate_new_reference_table
from count_correlation import get_correlation
import random
import os
import subprocess

def get_16s_counts():
    file = open("D:\DP\data\\GG_16S_counts.tab", "r")
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



new_reference = "D:\\DP\\data\\test\\GG_ko_counts_new.tab"

#settings
REF_FILE_NAME = "D:\\DP\\data\\GG_ko_counts.tab"
KOS_FILE_NAME = "D:\\DP\\data\\distance_matrix_known.tab"
NUMBER_OF_SAMPLES = 2  # kolko sa ma vygenerovat a testovat vzoriek
MAX_ABUNDANCE = 15  # obmedzuje maximalnu nahodne genrovanu abundanciu bakterie
REFERENCE_SIZE_TO_TEST = [90,80]  # ake percento zachovania referencnej tabulky sa ma testovat
NUMBER_OF_TESTS = 2  # kolkokrat sa ma kazde nastavenie otestovat
CORRELATION_FILE = "D:\\DP\\data\\correlations_presence_01.tab"
GENERATE = False  # ci sa maju aj generovat nove testove subory
COMPUTE = True  # ci sa maju aj pocitat nove veci (nie ked menime len vysledok korelacie)
ko_tresholds = ["100","90","80"]



print("Welcome in the Picrust testing framework! Let's begin.")
print("*******************************")
if GENERATE:
    print("Generating " + str(NUMBER_OF_SAMPLES) + " samples with expected results:")
    # generate samples
    sample_index = 0
    counts_16s = get_16s_counts()
    for sample_index in range(NUMBER_OF_SAMPLES):
        random.seed()
        sample_filename = "D:\\DP\\data\\test\\sample" + str(sample_index) + ".tab"
        sample_file = open(sample_filename, "w")
        # header
        sample_file.write("#OTU ID\tstaggered\teven\n")
        # sample
        ids_file = open(KOS_FILE_NAME, "r")
        line = ids_file.readline()
        # get kos dictionary
        ids_with_profile = line.split()

        ref_file = open(REF_FILE_NAME, "r")
        line = ref_file.readline()
        kos = line.split()[1:]
        kos_score = []
        for _ in range(len(kos)+1):
            kos_score.append((0.0, 0.0))

        while True:
            line = ref_file.readline()
            if not line:
                break
            if line.isspace():
                continue
            splitted_line = line.split()
            bacteria_code = splitted_line[0]
            if bacteria_code not in ids_with_profile:
                continue
            # generate random abundance
            staggered = float(random.randint(0, MAX_ABUNDANCE))
            even = float(random.randint(0, MAX_ABUNDANCE))
            sample_line = bacteria_code + "\t"
            sample_line += str(staggered) + "\t"
            sample_line += str(even)
            sample_file.write(sample_line + "\n")
            # count the resulting abundance
            index = 0
            abundance = 1
            if bacteria_code in counts_16s:
                abundance = counts_16s[bacteria_code]
            for ko_count in splitted_line[1:]:
                old_tuple = list(kos_score[index])
                old_tuple[0] += float(ko_count) * staggered / abundance
                old_tuple[1] += float(ko_count) * even / abundance
                new_tuple = (old_tuple[0], old_tuple[1])
                kos_score[index] = new_tuple
                index += 1
        sample_file.close()
        """
        # convert sample to biom (picrust needs biom input)
        filename = "test_data/sample" + str(sample_index)
        a = subprocess.Popen("biom convert -i "+filename+".txt -o "+filename+".biom --to-json --table-type='OTU table'",shell=True)
        a.wait()"""

        # write expected result
        expected_result = open("D:\\DP\\data\\test\\expected"+str(sample_index)+".tab", "w")
        expected_result.write("#KO ID\tstaggered\teven\n")
        index = 0
        for ko_name in kos:
            to_write = ko_name + "\t"
            to_write += "{0:.1f}".format(kos_score[index][0]) + "\t"
            to_write += "{0:.1f}".format(kos_score[index][1]) + "\n"
            expected_result.write(to_write)
            index += 1
        print("Sample " + str(sample_index+1) + " done! " + str(NUMBER_OF_SAMPLES - sample_index - 1) + " to go!")

print("*******************************")
print("Testing begins!")
# test picrustu
# pre kazdu velkost povodnej referencnej tabulky
correl_file = open(CORRELATION_FILE, "w")
correl_file.write("method\tscore_matrix\ttest_number\tsample\treference\tko_treshold\tstaggered\teven\n")
for new_size in REFERENCE_SIZE_TO_TEST:
    print("________________________________")
    print("We are testing " + str(new_size) + "% of the reference table")
    for i in range(NUMBER_OF_TESTS):
        print("_____ Test " + str(i+1) + " of " + str(NUMBER_OF_TESTS) + ":")
        print("Generating new reference table")
        if COMPUTE:
            generate_new_reference_table(100 - new_size, REF_FILE_NAME, new_reference)
        for j in range(NUMBER_OF_SAMPLES):

            print("Testing sample " + str(j))
            # test sample in
            sample_filename = "D:\DP\data\\\\test\sample"+str(j)+".tab"
            expected_filename = "D:\DP\data\\\\test\expected" + str(j) + ".tab"
            result_filename = "D:\DP\data\\\\test\\\\result_ref" + str(new_size)+"_sample"+str(j)+"_test"+str(i)+".tab"

            if COMPUTE:
                settings_file = open("../settings.py", "r")
                lines = settings_file.readlines()
                lines[13] = "INPUT_FILENAME = \"" + sample_filename + "\"\n"
                lines[35] = "RESULT_FILENAME = \"" + result_filename + "\"\n"
                settings_file.close()

                settings_file = open("../settings.py", "w")
                settings_file.write("".join(lines))
                settings_file.close()

                a = subprocess.Popen(
                    "python D:\DP\main.py",
                    shell=True)
                a.wait()
                print("Processed")

            for t in ko_tresholds:
                coef_staggered, coef_even = get_correlation(expected_filename, result_filename, t)

                correl_file.write("normal_average\talignment_difference\t"+str(i+1)
                                  +"\t"+str(j)+"\t"+str(new_size)+"\t" + t + "\t"
                                  +str(coef_staggered)+"\t"+str(coef_even) + "\n")
                print("St: " + str(coef_staggered) + " ev: " + str(coef_even))
    correl_file.close()
    correl_file = open(CORRELATION_FILE, "a")