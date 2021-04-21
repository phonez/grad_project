#!/bin/bash

# Check params file if there is "N    CA   C    O "
# This is a confusing bug in molfile_to_params_polymer.py

cd $(dirname $0)
cd ..
params_path=$PWD"/output/params/"
cd ${params_path}

chi_num=1
wrong_chi='N    CA   C    O  '
has_wrong_chi=0
has_next_chi=0

for params in *.params
do
    tmp=`grep "${wrong_chi}" ${params}`
    result=$?

    if (( ${result}==${has_wrong_chi} ))
    then
        echo "Found wrong chi in ${params}!"
        
        chi_num=`echo ${tmp:4:1}`
        sed -i '/'"${wrong_chi}"'/d' ${params} # delete wrong chi

        (( chi_num++ ))
        next_chi="CHI ${chi_num}"
        tmp=`grep "${next_chi}" ${params}`
        result=$?

        while (( ${result}==${has_next_chi} ))
        do
            correct_chi_num=`expr ${chi_num} - 1`
            sed -i 's/'"${next_chi}"'/'"CHI ${correct_chi_num}"'/g' ${params} # set correct number of other chis

            (( chi_num++ ))
            next_chi="CHI ${chi_num}"
            tmp=`grep "${next_chi}" ${params}`
            result=$?
        done

        echo "Error fixed."
    else
        echo "Check complete: no wrong chi found in ${params}"
    fi
    chi_num=1
done