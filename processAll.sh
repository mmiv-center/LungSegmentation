#!/bin/bash

input="/Users/hauke/Documents/WORK/LIDC-IDRI/LIDC-IDRI"
output="/tmp"
dd=`pwd`
cd "$input"
dirs=$(find . -mindepth 1 -maxdepth 1 -type d -print)
echo "${#dirs[@]} files"
cd ${dd}

# remove the newlines (spaces are not allowed in the path)
dirs=`echo $dirs | tr '\n' ' '`
# run 6 segmentations in parallel on the current hardware
parallel -j 4 ./LungSegmentation ${input}/{} ${output}/{} ::: $dirs
