#!/bin/bash

file=$1
ext=$(echo $file | rev | cut -d\/ -f1 | cut -d. -f1 | rev)
name=$(echo $file | rev | cut -d\/ -f1 | cut -d. -f2- | rev)

mkdir -p ./results/$name
output_path="./results/$name/"
output_prefix="$output_path/$name"

echo "processing $name.$ext"

./install/bin/ttkRipsMultiParameterSamplingCmd -i $file -o $output_prefix -fitType 0 -P 0.5

echo "creating paraview state file"
display_file=$output_path/display_$name.pvsm
cp ./src/display.pvsm $display_file


sed -i "s/NAME/$name/" $display_file
