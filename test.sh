#!/bin/bash

if [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
	echo "Script for running framework for functional profile prediction

	Format: ./test.sh python3-alias [-m method-name] [-i input-file] [-o output-file]

	Params:
		python3-alias = python3 alias on yout computer (usually python3 or python)
		method-name = which method tu use for functional profile prediction (available: alignment_simple, random, tree, treshold, weighted, regression_by_sequence)
		input-file = path to file containing input otu table
		output-file = desired path to the output file
		"
	exit
fi

cd src
ALIAS=$1
if [ "$1" == "" ]; then
   echo "No python3 alias specified!"
   exit
fi

shift

echo "Params: $@"
echo "Computing..."
$ALIAS main.py $@
echo "Finished!"
