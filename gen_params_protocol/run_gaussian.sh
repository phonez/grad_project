#!/bin/bash
cur=0
total=`ls ./*.gjf|wc -l`
for inf in *.gjf
do
((cur++))
echo Running ${inf} ... \($cur of $total\)
g09 < ${inf} > ${inf//gjf/log}
echo ${inf} has finished
done

