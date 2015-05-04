#!/bin/bash

dir=${HOME}/Dropbox/GitHub/Load/shell/log

trims=`find $dir -name 'trim.sh.e*'`
rsems=`find $dir -name 'rsem.sh.e*'`

for file in $trims; do
  name=`basename $file`
# echo $name >> 1
  grep "Both Surviving" $name >> trim.log
done 

for file in $rsems; do
  name=`basename $file`
# echo $name >> 1
  grep "at least" $name >> rsem.log
done 

