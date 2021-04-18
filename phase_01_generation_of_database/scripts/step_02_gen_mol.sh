#!/bin/bash

# Convert gaussian output log file to mol file.

cd $(dirname $0)
cd ..
log_path=$PWD"/output/log/"
mol_path=$PWD"/output/mol/"
cd ${log_path}

for log_file in *.log
do
obabel -i g09 ${log_file} -o mol -O ${log_file//log/mol}
mv ${log_file//log/mol} ${mol_path}
done
  