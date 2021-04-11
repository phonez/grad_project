#!/bin/bash

# Convert modified mol file to params file.

cd $(dirname $0)
cd ..
modified_mol_path=$PWD"/out/modified_mol/"
params_path=$PWD"/out/params/"
cd ${modified_mol_path}

# name the ncaa with three characters previously
for modified_mol_file in *.mol
do
log_file=${params_path}${modified_mol_file//mol/log}
python2 /home/zhangf/grad_project/gen_params_protocol/modified_mol_file_to_params_polymer.py\
 --clobber --polymer --no-pdb --name ${modified_mol_file//.mol} ${modified_mol_file} > ${log_file}

mv ${modified_mol_file//mol/params} ${params_path}
rm ${log_file}
done