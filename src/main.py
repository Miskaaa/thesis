# import parts from the library
import settings
import input_parser
import otu_analyzer
import unknown_otu_resolver
import output_generator
import common_functions
import sys

from importlib import reload

reload(input_parser)
reload(otu_analyzer)
reload(unknown_otu_resolver)
reload(output_generator)

import optparse

parser = optparse.OptionParser()

parser.add_option('-i', '--input',
    action="store", dest="input_file",
    help="input file filename", default=settings.INPUT_FILENAME)
parser.add_option('-o', '--output',
    action="store", dest="output_file",
    help="output file filename", default=settings.RESULT_FILENAME)
parser.add_option('-m', '--method',
    action="store", dest="method",
    help="method for functional profile prediction", default=settings.UNKNOWN_METHOD)

options, args = parser.parse_args()

# load input otu table
parsed_otus = {}
parsed_otus = input_parser.parse_input(options.input_file)

counts_16s = common_functions.get_16s_counts()

# find known ko profiles in reference table
unknown_otus, known_otus, result_dic = otu_analyzer.analyze_otus(parsed_otus, counts_16s)

# get the most similiar known otus to unknown otus
unknown_result = {}
unknown_result = unknown_otu_resolver.resolve_unknown(options.method, unknown_otus, known_otus, counts_16s)
	
# generate output
output_generator.make_output(result_dic, unknown_result, options.output_file)

