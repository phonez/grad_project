#!/bin/bash

# Incorparate all the scripts except gen_rotlib_in.py
# Convert mol2 file to params file.

cd $(dirname $0)

echo "Generating gjf files..."
python3 ./step_00_gen_gjf.py
if [ $? == 0 ]
then
    echo "Done."
else
    echo "Failed in some cases."
	exit 1
fi

echo "Running gaussian..."
if [ $? == 0 ]
then
    echo "Done."
else
    echo "Failed in some cases."
	exit 1
fi

echo "Generating mol files..."
bash ./step_02_gen_mol.sh
if [ $? == 0 ]
then
    echo "Done."
else
    echo "Failed in some cases."
    exit 1
fi

echo "Generating modified mol files..."
python3 ./step_03_gen_modified_mol.py
if [ $? == 0 ]
then
    echo "Done."
else
    echo "Failed in some cases."
    exit 1
fi

echo "Generating params files..."
bash ./step_04_gen_params.sh
if [ $? == 0 ]
then
    echo "Done."
else
    echo "Failed in some cases."
    exit 1
fi

echo "Checking chis in params files..."
bash ./step_05_check_chi.sh
if [ $? == 0 ]
then
    echo "Done."
else
    echo "Failed in some cases."
    exit 1
fi