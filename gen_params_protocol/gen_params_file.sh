#!/bin/bash

# modified mol file --> params file

modified_mol_filedir="/home/zhangf/grad_project/pipeline/modified_mol/"
params_filedir="/home/zhangf/grad_project/pipeline/params/"

cd ${modified_mol_filedir}

for molfile in *.mol
do
python2 /home/zhangf/grad_project/gen_params_protocol/molfile_to_params_polymer.py\
 --clobber --polymer --no-pdb --name ${molfile//.mol} -k ${molfile//mol/kin} ${molfile} # name the ncaa with three characters previously

mv ${molfile//mol/kin} ${params_filedir}
mv ${molfile//mol/params} ${params_filedir} 
done