#!/bin/bash

# gaussian output log file --> mol file

log_filedir="/home/zhangf/grad_project/pipeline/log/"
mol_filedir="/home/zhangf/grad_project/pipeline/mol/"

cd ${log_filedir}

for logfile in *.log
do
obabel -i g09 ${logfile} -o mol -O ${logfile//log/mol}
mv ${logfile//log/mol} ${mol_filedir}
done
  