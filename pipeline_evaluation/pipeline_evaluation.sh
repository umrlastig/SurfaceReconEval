#!/bin/bash

# Created on 08 January 2021
# Purpose: automation of the process of surface reconstruction evaluation


# store arguments
args=("$@")
declare -a FILES; declare -a DIRS


extract_cmd_line_params() {
  args=$1
  nbArgs=${#args[@]}

  # Default values:
  param_file="./eval.txt" # file to set evaluation parameters
  verboseFlag=""
  exampleFlag=""
  ext=".ply"

  # Browse arguments:
  for (( i=0; i<$nbArgs; i++ )); do
    if [ "${args[${i}]}" = "--example" ] || [ "${args[${i}]}" = "-e" ]; then
      echo "Running example"
      exampleFlag="--example"
      param_file="./example_eval.txt" # file to set evaluation parameters
    elif [ "${args[${i}]}" = "--verbose" ] || [ "${args[${i}]}" = "-v" ]; then
      echo "verbose flag active"
      verboseFlag="-v"
    else
      echo "Error: argument '${args[${i}]}' not valid"
    fi
  done
}

initiate_files_and_dirs() {
  ## Poisson disk sampling script (meshlabserver)
  pds_script="PDS.mlx"
  ## c++ executables:
  mesh_alpha_exec="./../compute_mesh_alpha"
  remove_points_too_far_from_P="./../remove_points_too_far_from_P"
  mean_and_max_distance_from_P_to_mesh="./../mean_and_max_distance_from_P_to_mesh"

  strAlpha="_alpha_"

  # Directories
  DIR_LiDAR_scan="LiDAR_scan"
  DIR_mesh="mesh"
  DIR_mesh_alpha="mesh_alpha"
  DIR_mesh_alpha_PDS="mesh_alpha_PDS"
  DIR_mesh_alpha_PDS_processed="mesh_alpha_PDS_processed"
  DIR_results="results"
  DIR_logs="logs"
  DIR_EVALUATIONS="evaluations"

  FILES+=($pds_script)
  FILES+=($mesh_alpha_exec)
  FILES+=($remove_points_too_far_from_P)
  FILES+=($mean_and_max_distance_from_P_to_mesh)

  DIRS+=($DIR_LiDAR_scan)
  DIRS+=($DIR_mesh)
  DIRS+=($DIR_mesh_alpha)
  DIRS+=($DIR_mesh_alpha_PDS)
  DIRS+=($DIR_mesh_alpha_PDS_processed)
  DIRS+=($DIR_results)
  DIRS+=($DIR_logs)
}

update_PDS_radius() {
  newR=$1
  left="<Param type=\"RichAbsPerc\" value=\""
  right="\" min=\"0\" name=\"Radius\" max=\"356.693\"\/>"
  cmd="sed -i 's/${left}.*${right}/${left}${newR}${right}/g' ${pds_script}"
  eval ${cmd}
}

generate_log_res_file() {
  file=$1

  evalBasename=$(basename -- "$file") # removes path
  evalBasename="${evalBasename%.*}" # removes extension

  logName="_LOG"; extLog=".txt"
  resName="_RESULTS"; extRes=".txt"

  logFile="${EVALUATION}/${DIR_logs}/${evalBasename}${logName}${extLog}"
  resFile="${EVALUATION}/${DIR_results}/${evalBasename}${resName}${extRes}"
}

process_current_mesh() {
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
      --input-mesh ${meshFile}\
      --point-cloud ${GT_LiDAR}\
      --alpha ${alpha}\
      --output-mesh ${outMesh}\
      ${verboseFlag}"
    # echo "Executing ${cmd}"
    eval ${cmd}


    # 2. Poisson disk sampling
    echo -e "\n--> COMPUTING POISSON DISK SAMPLING\n"
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
      --input-sampling ${inPointSet}\
      --input-P ${GT_LiDAR}\
      --alpha ${alpha}\
      --output-file ${outPointSet}\
      ${verboseFlag}"
    # echo "executing: ${cmd}"
    eval ${cmd}
  done
}

run_distance_computation() {
  pcdFile=$1
  meshFile=$2

  # check existence of files
  ./check_file_exists.sh "${pcdFile}"; if [ $? -ne 0 ]; then exit 1; fi
  ./check_file_exists.sh "${meshFile}"; if [ $? -ne 0 ]; then exit 1; fi

  cmd="${mean_and_max_distance_from_P_to_mesh}\
    --point-cloud ${pcdFile}\
    --mesh ${meshFile}\
    ${verboseFlag}"
  # echo "executing ${cmd}"
  eval ${cmd}
}

compute_distances() {
  file=$1 # the one from ${DIR_mesh}

  # generate a file name for results storing:
  resName="_RESULTS"; extRes=".txt"
  extFile="${file##*.}"
  evalBasename=$(basename -- "$file") # removes path
  evalBasename="${evalBasename%.*}" # removes extension

  echo -e "\n\n--> ASSESSMENT OF: '$(basename -- "$file")'"



  for alpha in "${TAB_ALPHAS[@]}"; do
    echo -e "\nalpha : ${alpha}"
    fileName="${evalBasename}${strAlpha}${alpha}${ext}"

    # PRECISION: pcd=RECON_PDS / mesh=GT_alpha
    echo "[PRECISION]"
    pcdFile="${EVALUATION}/${DIR_mesh_alpha_PDS_processed}/${fileName}"
    meshFile="${EVALUATION}/${DIR_mesh_alpha}/${fileNameGT}${strAlpha}${alpha}${ext}"
    run_distance_computation ${pcdFile} ${meshFile}

    # RECALL: pcd=GT_PDS / mesh=RECON_alpha
    echo "[RECALL]"
    pcdFile="${EVALUATION}/${DIR_mesh_alpha_PDS_processed}/${fileNameGT}${strAlpha}${alpha}${ext}"
    meshFile="${EVALUATION}/${DIR_mesh_alpha}/${fileName}"
    run_distance_computation ${pcdFile} ${meshFile}
  done
  echo -e "\n# assessment done"
}

read_parameters() {
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
    # read ground-truth mesh file:
    if [[ ${line} =~ "GT_MESH" ]]; then
      echo "  > ${line}"
      GT_MESH=$(echo $line| cut -d'=' -f 2)
      FILES+=($GT_MESH)
    fi

    # read LiDAR scan file:
    if [[ ${line} =~ "GT_LiDAR" ]]; then
      echo "  > ${line}"
      GT_LiDAR=$(echo $line| cut -d'=' -f 2)
      FILES+=($GT_LiDAR)
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

    # read poisson-disk sampling radius and update PDS script:
    if [[ ${line} =~ "PDS_RADIUS" ]]; then
      echo "  > ${line}"
      PDS_RADIUS=$(echo $line| cut -d'=' -f 2)
    fi
    update_PDS_radius ${PDS_RADIUS}
  done < "${param_file}"
}

set_up_environment() {
  # Check existence of files
  echo -ne "\nChecking important files: "
  for file in "${FILES[@]}" ; do
    ./check_file_exists.sh "${file}";
    if [ $? -ne 0 ]; then exit 1; fi # exit on failure of child script
  done
  echo -ne "OK\n"


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
  cp ${GT_MESH} ${new_GT_MESH} && echo -e " + copied ${GT_MESH}\n    > and renamed it to: '${fileNameGT}.${extensionGT}'"
  GT_MESH=${new_GT_MESH} # update
  ## LiDAR :
  extensionLiDAR="${GT_LiDAR##*.}"
  fileNameLiDAR="lidar_scan"
  new_LiDAR="${EVALUATION}/${DIR_LiDAR_scan}/${fileNameLiDAR}.${extensionLiDAR}"
  cp ${GT_LiDAR} ${new_LiDAR} && echo -e " + copied ${GT_LiDAR}\n    > and renamed it to: '${fileNameLiDAR}.${extensionLiDAR}'"
  GT_LiDAR=${new_LiDAR} # update
}

chose_evaluated_meshes() {
  # Prompt the user for copying the meshes to be evaluated in the proper directory
  echo -e "\nCheckings done. Everything has been setup successfully :)"
  echo -e "---------------------------------------------------------\n\n\n"
  if [[ ${exampleFlag} == "--example" ]]; then
    echo "Running example"
    mesh1="example/PoissonRecon_example.ply"
    mesh2="example/SSDRecon_example.ply"
    cmd="cp ${mesh1} ${mesh2} ${EVALUATION}/${DIR_mesh}/"
    # echo "executing: ${cmd}"
    eval ${cmd}
    nb_of_results=$(ls ${EVALUATION}/${DIR_mesh} | wc -l)
    echo "$((${nb_of_results}-1)) file(s) have been found"
    read -p "Continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1
  else
    echo -e "IMPORTANT! Please COPY THE RESULTS of each surface reconstruction algorithm\
    you have run on '${GT_LiDAR}' and would like to assess into directory:"
    echo -e "\n   ===> $(pwd)/${EVALUATION}/${DIR_mesh}/ <===\n"
    read -p "Done? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1

    nb_of_results=$(ls ${EVALUATION}/${DIR_mesh} | wc -l)
    echo "$((${nb_of_results}-1)) file(s) have been found" # -1 because ground-truth was also copied automatically
    read -p "Continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1
  fi
  echo -e "Started the process"
}

mesh_processing() {
  ## Mesh processing
  for file in "${EVALUATION}/${DIR_mesh}"/*; do
    generate_log_res_file ${file}
    process_current_mesh ${file} 2>&1 | tee -a ${logFile}
  done
}

assessment() {
  ## Assessment
  for file in "${EVALUATION}/${DIR_mesh}"/*; do
    if [ ${file} != ${GT_MESH} ]; then
      generate_log_res_file ${file}
      compute_distances ${file} | tee ${resFile}
    fi
  done
}

main(){
  echo -e "\n---> STARTED EVALUATION PROTOCOL <---\n"
  extract_cmd_line_params ${args}
  initiate_files_and_dirs
  read_parameters
  set_up_environment
  chose_evaluated_meshes
  mesh_processing
  assessment
}

main

exit 0

