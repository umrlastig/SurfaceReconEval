#!/bin/bash

# Created on 08 January 2021
# Purpose: automation of the process of surface reconstruction evaluation


# Files
param_file="./example_eval.txt" # file to set evaluation parameters
pds_script="PDS.mlx" # Poisson disk sampling script (meshlabserver)
## c++ executables:
mesh_alpha_exec="./../mesh_alpha"
remove_points_too_far_from_P="./../remove_points_too_far_from_P"
mean_and_max_distance_from_P_to_mesh="./../mean_and_max_distance_from_P_to_mesh"

strAlpha="_alpha_"
ext=".ply"

echo -e "\n---> STARTED EVALUATION PROTOCOL <---\n"


###############
# Directories #
###############
declare -a DIRS
DIR_LiDAR_scan="LiDAR_scan"; DIRS+=($DIR_LiDAR_scan)
DIR_mesh="mesh"; DIRS+=($DIR_mesh)
DIR_mesh_alpha="mesh_alpha"; DIRS+=($DIR_mesh_alpha)
DIR_mesh_alpha_PDS="mesh_alpha_PDS"; DIRS+=($DIR_mesh_alpha_PDS)
DIR_mesh_alpha_PDS_processed="mesh_alpha_PDS_processed"; DIRS+=($DIR_mesh_alpha_PDS_processed)
DIR_EVALUATIONS="evaluations"
###############



#############
# Functions #
#############
func_process_mesh() {
  meshFile=$1

  echo -e "\n    ----> PROCESSING: $(basename -- "$meshFile") <----"

  evalBasename=$(basename -- "$meshFile") # removes path
  evalBasename="${evalBasename%.*}" # removes extension


  for alpha in "${TAB_ALPHAS[@]}"; do
    echo -e "\nalpha : ${alpha}"

    fileName="${evalBasename}${strAlpha}${alpha}${ext}"

    # 1. Compute mesh_alpha for all alphas
    echo -e "\n--> COMPUTING MESH_ALPHA"
    outMesh="${EVALUATION}/${DIR_mesh_alpha}/${fileName}"
    cmd="${mesh_alpha_exec}\
      ${meshFile}\
      ${GT_LiDAR}\
      ${alpha}\
      ${outMesh}"
    # echo "Executing ${cmd}"
    eval ${cmd}


    # 2. Poisson disk sampling
    echo -e "\n--> COMPUTING POISSON DISK SAMPLING"
    inMesh=${outMesh}
    outPointSet="${EVALUATION}/${DIR_mesh_alpha_PDS}/${fileName}"
    cmd="LC_ALL=C meshlabserver\
      -i ${inMesh}\
      -o ${outPointSet}\
      -s ${pds_script}"
    # echo "Executing ${cmd}"
    eval ${cmd}


    # 3. Remove points from PDS that are too far from LiDAR scan
    echo -e "\n--> POST-PROCESSING: REMOVING POINTS TOO FAR FROM LiDAR SCAN\n"
    inPointSet=${outPointSet}
    outPointSet="${EVALUATION}/${DIR_mesh_alpha_PDS_processed}/${fileName}"
    cmd="${remove_points_too_far_from_P}\
      ${inPointSet}\
      ${GT_LiDAR}\
      ${alpha}\
      ${outPointSet}"
    # echo "executing: ${cmd}"
    eval ${cmd}
  done
}

func_run_distance_computation() {
  pcdFile=$1
  meshFile=$2
  outFile=$3

  # check existence of files
  ./check_file_exists.sh "${pcdFile}"; if [ $? -ne 0 ]; then exit 1; fi
  ./check_file_exists.sh "${meshFile}"; if [ $? -ne 0 ]; then exit 1; fi

  cmd="${mean_and_max_distance_from_P_to_mesh}\
    ${pcdFile}\
    ${meshFile}\
    ${outFile}"
  # echo "executing ${cmd}"
  eval ${cmd}
}

func_compute_distances() {
  file=$1 # the one from ${DIR_mesh}
  outFile="someFile.gnumeric"
  echo -e "\n\n--> ASSESSMENT OF: '$(basename -- "$file")'"


  evalBasename=$(basename -- "$file") # removes path
  evalBasename="${evalBasename%.*}" # removes extension

  for alpha in "${TAB_ALPHAS[@]}"; do
    echo -e "\nalpha : ${alpha}"
    fileName="${evalBasename}${strAlpha}${alpha}${ext}"

    # PRECISION: pcd=RECON_PDS / mesh=GT_alpha
    echo "[PRECISION]"
    pcdFile="${EVALUATION}/${DIR_mesh_alpha_PDS_processed}/${fileName}"
    meshFile="${EVALUATION}/${DIR_mesh_alpha}/${fileNameGT}${strAlpha}${alpha}${ext}"
    func_run_distance_computation ${pcdFile} ${meshFile} ${outFile}

    # RECALL: pcd=GT_PDS / mesh=RECON_alpha
    echo "[RECALL]"
    pcdFile="${EVALUATION}/${DIR_mesh_alpha_PDS_processed}/${fileNameGT}${strAlpha}${alpha}${ext}"
    meshFile="${EVALUATION}/${DIR_mesh_alpha}/${fileName}"
    func_run_distance_computation ${pcdFile} ${meshFile} ${outFile}
  done

}
#############


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
  if [[ ${line} =~ "GT_MESH" ]]; then
    echo "  > ${line}"
    GT_MESH=$(echo $line| cut -d'=' -f 2)
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


# Check existence of files
echo -e "\nChecking important files:"
for file in ${GT_MESH} ${GT_LiDAR}; do
  ./check_file_exists.sh "${file}";
  if [ $? -ne 0 ]; then exit 1; fi # exit on failure of child script
done


# Setup sub-directories for a new evaluation:
echo -e "\nSetting up sub-directories for a new evaluation:"
EVALUATION=$(basename -- "$GT_MESH") # removes path
EVALUATION="${EVALUATION%.*}" # removes extension
dt=$(date '+D%d_M%m_Y%Y_T%H_%M_%S'); EVALUATION+="_${dt}" # adds date and time
EVALUATION="${DIR_EVALUATIONS}/${EVALUATION}"
# echo $EVALUATION
mkdir -p ${EVALUATION}; echo -e " + NEW EVALUATION directory: '${EVALUATION}/'"
for dir in "${DIRS[@]}"; do
  mkdir "${EVALUATION}/${dir}";
  echo "   + sub-directory '${dir}/'"
done


# Copy important files to newly created environment:
## GT :
echo -e "\nCopying input files to newly created environment:"
extensionGT="${GT_MESH##*.}"
fileNameGT="ground_truth"
new_GT_MESH="${EVALUATION}/${DIR_mesh}/${fileNameGT}.${extensionGT}"
cp ${GT_MESH} ${new_GT_MESH}
echo -e " + copied ${GT_MESH}\n    > and renamed it to: '${fileNameGT}.${extensionGT}'"
GT_MESH=${new_GT_MESH} # update
## LiDAR :
extensionLiDAR="${GT_LiDAR##*.}"
fileNameLiDAR="lidar_scan"
new_LiDAR="${EVALUATION}/${DIR_LiDAR_scan}/${fileNameLiDAR}.${extensionLiDAR}"
cp ${GT_LiDAR} ${new_LiDAR}
echo -e " + copied ${GT_LiDAR}\n    > and renamed it to: '${fileNameLiDAR}.${extensionLiDAR}'"
GT_LiDAR=${new_LiDAR} # update



# Prompt the user for copying the meshes to be evaluated in the proper directory
echo -e "\nCheckings done. Everything has been setup successfully :)"
echo -e "---------------------------------------------------------\n\n\n"
echo -e "IMPORTANT! Please COPY THE RESULTS of each surface reconstruction algorithm\
  you have run on '${GT_LiDAR}' and would like to assess into directory:"
echo -e "\n   ===> ${EVALUATION}/${DIR_mesh}/ <===\n"
read -p "Done? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1

nb_of_results=$(ls ${EVALUATION}/${DIR_mesh} | wc -l)
echo "${nb_of_results} files have been found, including ground-truth"
read -p "Continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1

echo -e "Started the process"


#############
# EXECUTION #
#############
for file in "${EVALUATION}/${DIR_mesh}"/*; do
  func_process_mesh ${file}
done

for file in "${EVALUATION}/${DIR_mesh}"/*; do
  if [ ${file} != ${GT_MESH} ]; then # do not evaluate ground-truth
    func_compute_distances ${file}
  fi
done
#############

exit 0

