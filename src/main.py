# import parts from the library
import settings
import input_parser
import otu_analyzer
import unknown_otu_resolver
import output_generator
import common_functions


from importlib import reload

reload(input_parser)
reload(otu_analyzer)
reload(unknown_otu_resolver)
reload(output_generator)

debug = False

# @Todo:
# .biom compatibility

# load input otu table
parsed_otus = input_parser.parse_input(settings.INPUT_FILENAME)
if debug:
    print("Parsed otus: ")
    print(parsed_otus)

counts_16s = common_functions.get_16s_counts()
if debug:
    print(counts_16s)

# find known ko profiles in reference table
unknown_otus, known_otus, result_dic = otu_analyzer.analyze_otus(parsed_otus, counts_16s)
if debug:
    print("result after known:")
    print(result_dic)
    print("unknown otus:")
    print(unknown_otus)

# get the most similiar known otus to unknown otus
unknown_result = unknown_otu_resolver.resolve_unknown(settings.UNKNOWN_METHOD, unknown_otus, known_otus, counts_16s)
if debug:
    print("result after known:")
    print(result_dic)
# generate output
output_generator.make_output(result_dic, unknown_result)