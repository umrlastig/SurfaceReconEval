#!/bin/bash

# Created on 08 January 2021
# Purpose: automation of the process of surface reconstruction evaluation

param_file="./example_eval.txt" # file to set evaluation parameters

echo -e "\n---> STARTED EVALUATION PROTOCOL <---\n"

# Directories
declare -a DIRS
DIR_GT_LiDAR="GT_LiDAR/"; DIRS+=($DIR_GT_LiDAR)
DIR_GT_mesh="GT_mesh/"; DIRS+=($DIR_GT_mesh)
DIR_GT_mesh_alpha="GT_mesh_alpha/"; DIRS+=($DIR_GT_mesh_alpha)
DIR_GT_mesh_alpha_PDS="GT_mesh_alpha_PDS/"; DIRS+=($DIR_GT_mesh_alpha_PDS)
DIR_GT_mesh_alpha_PDS_processed="GT_mesh_alpha_PDS_processed/"; DIRS+=($DIR_GT_mesh_alpha_PDS_processed)
DIR_RECON_mesh="RECON_mesh/"; DIRS+=($DIR_RECON_mesh)
DIR_RECON_mesh_alpha="RECON_mesh_alpha/"; DIRS+=($DIR_RECON_mesh_alpha)
DIR_RECON_sampling="RECON_sampling/"; DIRS+=($DIR_RECON_sampling)
DIR_EVALUATIONS="evaluations"
##############
# PARAMETERS #
##############
# Check if parameters file exists and read from it
if [[ -f ${param_file} ]]
then
  echo "Reading parameters from: '${param_file}'"
else
  echo -e "ERROR: Evaluation file '${param_file}' not found!\nProgram ended unsuccessfully."
  exit 1
fi

while IFS= read -r line
do
  #echo "${line}"

  # read ground-truth mesh file:
  if [[ ${line} =~ "GT_FILE" ]]; then
    echo "  > ${line}"
    GT_FILE=$(echo $line| cut -d'=' -f 2)
  fi

  # read LiDAR scan file:
  if [[ ${line} =~ "GT_LiDAR" ]]; then
    echo "  > ${line}"
    GT_LiDAR=$(echo $line| cut -d'=' -f 2)
  fi

  # read alpha values:
  if [[ ${line} =~ "ALPHA" ]]; then
    echo "  > ${line}"
    ALPHAS=$(echo $line| cut -d'=' -f 2)
    IFS='|' read -r -a TAB_ALPHAS <<< "$ALPHAS"
    # print individual values of alpha:
    for index in "${!TAB_ALPHAS[@]}"
    do
      echo "    - alpha_${index} => ${TAB_ALPHAS[index]}"
    done
  fi

  # read poisson-disk sampling radius:
  if [[ ${line} =~ "PDS_RADIUS" ]]; then
    echo "  > ${line}"
    PDS_RADIUS=$(echo $line| cut -d'=' -f 2)
  fi
done < "${param_file}"
##############
##############

# Check existence of files
echo -e "\nChecking important files:"
for file in ${GT_FILE} ${GT_LiDAR}
do
  ./check_file_exists.sh "${file}";
  if [ $? -ne 0 ]; then exit 1; fi # exit on failure of child script
done

# Setup sub-directories for a new evaluation:
echo -e "\nSetting up sub-directories for a new evaluation:"
EVALUATION=$(echo $GT_FILE| cut -d'.' -f 1) # removes extension from $GT_FILE
EVALUATION=$(basename -- "$EVALUATION") # removes path
dt=$(date '+D%d_M%m_Y%Y_T%H_%M_%S'); EVALUATION+="_${dt}" # adds date and time
EVALUATION="${DIR_EVALUATIONS}/${EVALUATION}"
echo $EVALUATION
mkdir -p ${EVALUATION}; echo -e " + NEW EVALUATION directory: '${EVALUATION}/'"
for dir in "${DIRS[@]}"; do
  mkdir "${EVALUATION}/${dir}";
  echo "   + sub-directory '${dir}'"
done

# Copy important files to newly created environment:
echo -e "\nCopying input files to newly created environment:"
cp ${GT_FILE} "${EVALUATION}/${DIR_GT_mesh}"; echo " + copied ${GT_FILE}"
cp ${GT_LiDAR} "${EVALUATION}/${DIR_GT_LiDAR}"; echo " + copied ${GT_LiDAR}"

echo -e "\nCheckings done. Everything has been setup successfully :)"
echo -e "---------------------------------------------------------\n\n\n"
echo -e "IMPORTANT! Please COPY THE RESULTS of each surface reconstruction algorithm you have run on '${GT_LiDAR}' into directory:"
echo -e "\n   ===> ${EVALUATION}/${DIR_RECON_mesh} <===\n"
read -p "Done? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1

nb_of_results=$(ls ${EVALUATION}/${DIR_RECON_mesh} | wc -l)
echo "${nb_of_results} files have been found."
read -p "Continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1

echo "Started the process"




