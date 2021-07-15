#!/bin/bash
#Author: Charles Henry

for FILE in trees/*
do

echo $FILE

python3 real_test.py $FILE event150_event70.newick

done
