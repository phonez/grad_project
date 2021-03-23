#!/bin/bash

# check params file if there is "CHI 1  N    CA   C    O "

params_filedir="/home/zhangf/grad_project/pipeline/params/"

cd ${params_filedir}

chi_num=1
wrong_chi='CHI 1  N    CA   C    O  '
has_wrong_chi=0
has_next_chi=0

for paramsfile in *.params
do
tmp=`grep "${wrong_chi}" ${paramsfile}`

if [ $?==${has_wrong_chi} ]
then
    sed -i '/'"${wrong_chi}"'/d' ${paramsfile} # delete wrong chi
    
    (( chi_num++ ))
    next_chi="CHI ${chi_num}"
    tmp=`grep "${next_chi}" ${paramsfile}`
    result=$?
    
    while (( ${result}==${has_next_chi} ))
    do
        correct_chi_num=`expr ${chi_num} - 1`
        sed -i 's/'"${next_chi}"'/'"CHI ${correct_chi_num}"'/g' ${paramsfile} # set correct number of other chis

        (( chi_num++ ))
        next_chi="CHI ${chi_num}"
        tmp=`grep "${next_chi}" ${paramsfile}`
        result=$?
    done
fi
chi_num=1
done