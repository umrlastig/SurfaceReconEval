#!/bin/bash

# Created on 07 February 2022
# Purpose: automation of the process of surface reconstruction evaluation from real data


# store arguments
args=("$@")
executable="./evaluate_reconstruction_real_data"
visualEvalFlag="--visual-eval"

print_help(){
	echo -e "\n\n      ===> BATCH REAL-DATA EVALUATION <===\n\n"
	echo -e "--ground-truth-dir , -GT    |    Specify GROUND-TRUTH directory\n"
	echo -e "--recon-dir , -Rec          |    Specify RECONSTRUCTIONS directory\n"
	echo -e "--results-dir , -Res        |    Specify where to store NUMERICAL result files\n"
	echo -e "--visual-eval-dir , -Visu   |    Specify where to store VISUAL result files\n"
}

extract_cmd_line_params() {
	nbArgs=${#args[@]}

	# Default values:
	param_file="./eval.txt" # file to set evaluation parameters
	verboseFlag=""
	exampleFlag=""
	ext=".ply"

  # Browse arguments:
	for (( i=0; i<$nbArgs; i++ )); do
		if [ "${args[${i}]}" = "--ground-truth-dir" ] || [ "${args[${i}]}" = "-GT" ]; then
			((i=i+1))
			GTDir="${args[${i}]}"
		elif [ "${args[${i}]}" = "--recon-dir" ] || [ "${args[${i}]}" = "-Rec" ]; then
			((i=i+1))
			ReconDir="${args[${i}]}"			
		elif [ "${args[${i}]}" = "--results-dir" ] || [ "${args[${i}]}" = "-Res" ]; then
			((i=i+1))
			ResultsDir="${args[${i}]}"
		elif [ "${args[${i}]}" = "--visual-eval-dir" ] || [ "${args[${i}]}" = "-Visu" ]; then
			((i=i+1))
			VisualResDir="${args[${i}]}"
		elif [ "${args[${i}]}" = "--help" ] || [ "${args[${i}]}" = "-h" ]; then
			print_help
			exit 0
		else
			echo "Error: argument '${args[${i}]}' not valid"
		fi
	done
}

set_gt_files_string() {
	firstFile="1"
	for file in ${GTDir}/*; do
		if [ "${firstFile}" = "1" ]; then
			GTFiles="${file}"
			firstFile=0
		else
			GTFiles="${GTFiles}***${file}"
		fi
	done
}

process_meshes(){
	verboseFlag="--verbose"
	for mesh in ${ReconDir}/*; do
		meshBasename=$(basename -- "$mesh") # removes path
		meshBasename="${meshBasename%.*}" # removes extension
		cmd="${executable}\
		  --in-gtPcd ${GTFiles}\
		  --in-recon ${mesh}\
		  --out-pcd ${VisualResDir}/${meshBasename}_VISU.ply\
		  ${visualEvalFlag}\
		  ${verboseFlag}\
		  2>&1 | tee ${ResultsDir}/${meshBasename}_NUM.txt"
		# echo -e "---> executing: ${cmd}\n"
		eval ${cmd}
	done
}

print_general_information(){
	echo -e "\n\n      ===> BATCH REAL-DATA EVALUATION <===\n\n"
	echo -e " - Ground-truth directory: ${GTDir}\n"
	echo -e " - Reconstructions directory: ${ReconDir}\n"
	echo -e " - Numerical results directory: ${ResultsDir}\n"
	echo -e " - Visual results directory: ${VisualResDir}\n\n"
	echo -e " - Ground-truth files: ${GTFiles}\n\n"
}


# main #
extract_cmd_line_params
set_gt_files_string
print_general_information
process_meshes


exit 0

