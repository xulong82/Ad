#!/bin/sh
#PBS -l nodes=1:ppn=20,walltime=59:59:59

par=${index}
/data/xwang/AD/SH/myrsem.sh $par

