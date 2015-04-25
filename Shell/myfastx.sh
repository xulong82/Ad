#!/bin/sh
#
# RNA-seq data analysis pipeline 
# Author: XuLong Wang (xulong.wang@jax.org)

# Quality control with FASTX toolkit
# Compute quality statistics

echo $0
begin=`date +%h`

module load fastx/0.13

dir_sou="/hpcdata/xwang/AD/howell_2014/retina"
dir_des="/hpcdata/xwang/AD/FASTX/retina"
#
files=`find $dir_sou -name *.fastq`

for filename1 in $files; do
  filename2=$(basename $filename1)
  echo $filename2
  fastx_quality_stats -Q33 -i "$dir_sou"/$filename2 -o "$dir_des"/${filename2/.fastq/.txt}
  # Quality boxplot for each cycle
  fastq_quality_boxplot_graph.sh -i "$dir_des"/${filename2/.fastq/.txt} -o "$dir_des"/${filename2/.fastq/_Qbox.png}
  # Nucleotide distribution for each cycle
  fastx_nucleotide_distribution_graph.sh -i "$dir_des"/${filename2/.fastq/.txt} -o "$dir_des"/${filename2/.fastq/_Ndist.png}
done

end=`date +%h`
echo $((end-begin))

