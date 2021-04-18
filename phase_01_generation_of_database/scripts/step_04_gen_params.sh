#!/bin/bash

# Convert modified mol file to params file.

cd $(dirname $0)
cd ..
modified_mol_path=$PWD"/output/modified_mol/"
params_path=$PWD"/output/params/"
script=$PWD"/scripts/step_04_molfile_to_params_polymer.py"
cd ${modified_mol_path}

# name the ncaa with three characters previously
for modified_mol_file in *.mol
do
log_file=${params_path}${modified_mol_file//mol/log}
python2 ${script} --clobber --polymer --no-pdb --name ${modified_mol_file//.mol} ${modified_mol_file} > ${log_file}

mv ${modified_mol_file//mol/params} ${params_path}
done