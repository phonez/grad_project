#!/bin/bash

# Check params file if there is "N    CA   C    O " (probably not)
# Dump pdb file from params file to check if it can be imported or if it records the correct topology info of ncaa 

cd $(dirname $0)
cd ..
params_path=$PWD"/output/params/"
cd ${params_path}

wrong_chi='N    CA   C    O  '
has_wrong_chi=0

for params in *.params
do
    tmp=`grep "${wrong_chi}" ${params}`
    result=$?

    if (( ${result}==${has_wrong_chi} ))
    then
        echo "Found wrong chi in ${params}!"
    else
        echo "Check complete: no wrong chi found in ${params}"
        python ../../scripts/step_05_dump_pdb.py ${params//.params}
    fi
done