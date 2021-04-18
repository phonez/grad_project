#!/bin/bash

# Run a batch of gaussian jobs.
# This script refers to http://sobereva.com/258.

cd $(dirname $0)
cd ..
gjf_path=$PWD"/output/gjf/"
log_path=$PWD"/output/log/"
cd ${gjf_path}

cur=0
total=`ls ./*.gjf|wc -l`
for gjf_file in *.gjf
do
((cur++))
echo Running ${gjf_file} ... \($cur of $total\)
time g09 < ${gjf_file} > ${gjf_file//gjf/log}
mv ${gjf_file//gjf/log} ${log_path}
echo ${gjf_file} has finished
done

