#!/bin/bash

dir1="/data/xwang/Load/FASTQ1"
dir2="/home/xwang/Dropbox/GitHub/Load/shell"

files=`find $dir1 -name '*X_R1.fastq'`

for name1 in $files; do
  name2=`basename $name1`
  name3=${name2/_R1.fastq/}
  echo $name3
  qsub -v arg=$name3 $dir2/rsem.sh
done
  
# files=`find $dir2 -name '*.sh.e*'`
# 
# for name1 in $files; do
#   name2=`basename $name1`
#   echo $name2 >> 1
#   grep "at least" $name2 >> 1
# done 

