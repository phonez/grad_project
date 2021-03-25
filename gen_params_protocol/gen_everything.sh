#!/bin/bash

:<<BLOCK
echo "Generating gjf files..."
python3 /home/zhangf/grad_project/gen_params_protocol/gen_gjf.py
if [ $? == 0 ]
then
    echo "Done."
else
    echo "Failed in some cases."
fi

echo "Running gaussian..."
# bash gen_log.sh
if [ $? == 0 ]
then
    echo "Done."
else
    echo "Failed in some cases."
fi
BLOCK

echo "Generating mol files..."
bash /home/zhangf/grad_project/gen_params_protocol/gen_mol.sh
if [ $? == 0 ]
then
    echo "Done."
else
    echo "Failed in some cases."
    exit 1
fi

echo "Generating modified mol files..."
python3 /home/zhangf/grad_project/gen_params_protocol/gen_modified_mol.py
if [ $? == 0 ]
then
    echo "Done."
else
    echo "Failed in some cases."
    exit 1
fi

echo "Generating params files..."
bash /home/zhangf/grad_project/gen_params_protocol/gen_params.sh
if [ $? == 0 ]
then
    echo "Done."
else
    echo "Failed in some cases."
    exit 1
fi

echo "Checking params files..."
bash /home/zhangf/grad_project/gen_params_protocol/check_chi.sh
if [ $? == 0 ]
then
    echo "Done."
else
    echo "Failed in some cases."
    exit 1
fi