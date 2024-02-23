#!/bin/bash

# Assuming your Python script is named "your_python_script.py"
python_script="/nucleus/projects/sraeisid/Bionano_vitual_kt/OMKar/report_KT.py"

# Loop from 1 to 40
for ((i=1; i<=9; i++))
do
    # Construct the input file name
    input_file="/nucleus/projects/sraeisid/Paul_Dremsek/00${i}_pipeline_results.zip-dir/${i}.txt"

    # Check if the input file exists
    if [ -e "$input_file" ]; then
        # echo "Running iteration $i with input file: $input_file"
        # Run the Python script with the current input file
        # input_file="-i /nucleus/projects/sraeisid/Paul_Dremsek/0${i}_pipeline_results.zip-dir/${i}.txt -centro /home/sraeisid/bionano_utill2/hg38_centro.txt -genes ../genes.bed -dec ../DDG2P_14_11_2023.csv"
        python "$python_script" -i "$input_file" -centro /home/sraeisid/bionano_utill2/hg38_centro.txt -genes ../genes.bed -dec ../DDG2P_14_11_2023.csv
    # else
        # echo "Input file $input_file not found."
    fi
done
