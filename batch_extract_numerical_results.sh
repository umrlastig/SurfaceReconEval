#!/bin/bash

# Created on 22 March 2022
# Purpose: automation of the process of extracting files assessing surface reconstruction from real data


# store arguments
args=("$@")

extract_cmd_line_params() {
	EvalDir="${args[0]}"
}

print_evaluation_information() {
	for file in ${EvalDir}/*_NUM.txt; do
		echo -e "\n\n > EVALUATED MESH: ${file}\n"
		# cmd="cat ${file} | grep '=> Mean distance'"; eval $cmd
		cat ${file} | grep '=> Mean distance'
	done
}


# main #
extract_cmd_line_params
print_evaluation_information


exit 0

