#!/bin/sh

#cd /data/xwang/AD/howell_2013
#files=`find ./ -name '*.fasta'`
#
#for name1 in $files; do
#  name2=`basename $name1`
#  name3=${name2/.fasta/}
#  echo $name3
#  mv $name2 ${name2/.fasta/.fastq}
#done

cd /data/xwang/AD/howell_new_trim
files=`find ./ -not -name '*unpaired*'`

for name1 in $files; do
  name2=`basename $name1`
  name3=${name2/_GES*/}
  echo $name3
done

