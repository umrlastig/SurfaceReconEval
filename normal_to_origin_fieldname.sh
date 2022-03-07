#!/bin/bash

file=$1
cmdX="sed -i 's/nx/x_origin/' ${file}"
cmdY="sed -i 's/ny/y_origin/' ${file}"
cmdZ="sed -i 's/nz/z_origin/' ${file}"

echo $cmdX; eval $cmdX
echo $cmdY; eval $cmdY
echo $cmdZ; eval $cmdZ

