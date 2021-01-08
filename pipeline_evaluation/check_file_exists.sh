#!/bin/bash

if [ "$#" -eq 1 ]; then
  file=$1
  # check if file exists:
  if [[ -f ${file} ]]; then
    echo " OK: '${file}' does exist"
  else
    echo -e "ERROR: '${file}' not found!\nProgram ended unsuccessfully :("
    exit 1
  fi

else
  echo "Please specify only one file to check"
  exit 1

fi
