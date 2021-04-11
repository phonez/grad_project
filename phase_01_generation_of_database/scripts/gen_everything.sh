#!/bin/bash

# Incorparate all the scripts except gen_rotlib_in.py
# Convert mol2 file to params file.

cd $(dirname $0)

echo "Generating gjf files..."
python3 ./gen_gjf.py
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
bash ./gen_mol.sh
if [ $? == 0 ]
then
    echo "Done."
else
    echo "Failed in some cases."
    exit 1
fi

echo "Generating modified mol files..."
python3 ./gen_modified_mol.py
if [ $? == 0 ]
then
    echo "Done."
else
    echo "Failed in some cases."
    exit 1
fi

echo "Generating params files..."
bash ./gen_params.sh
if [ $? == 0 ]
then
    echo "Done."
else
    echo "Failed in some cases."
    exit 1
fi

echo "Checking chis in params files..."
bash ./check_chi.sh
if [ $? == 0 ]
then
    echo "Done."
else
    echo "Failed in some cases."
    exit 1
fi