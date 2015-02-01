#!/bin/sh

file=$1

qsub -l nodes=1:ppn=20,walltime=59:59:59 $file

