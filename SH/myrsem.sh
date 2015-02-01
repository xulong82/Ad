#!/bin/sh
#
# RNA-seq data analysis pipeline
# Author: XuLong Wang (xulong.wang@jax.org)

# RSEM: alignment and summation

echo $0
begin=`date +%h`
 
module load bowtie/1.0.0

dir_sou="/hpcdata/xwang/AD/howell_2013_trim"
dir_des="/hpcdata/xwang/AD/RSEM_2013"
myref="/hpcdata/xwang/Rsem/GRCm38"

#index=$1
#files=`find $dir_sou -name *_L00"$index"_R1.fastq`

files=`find $dir_sou -name '*[0-9]_R1.fastq'`

for filename1 in $files; do
  filename2=`basename $filename1`
  filename3=${filename2/_R1.fastq/}
  echo $filename3
  rsem-calculate-expression -p 20 \
  			    --phred33-quals \
  			    --forward-prob 0.5 \
  			    --paired-end \
                            "$dir_sou"/"$filename3"_R1.fastq \
                            "$dir_sou"/"$filename3"_R2.fastq \
  			    "$myref" \
  			    "$dir_des"/"$filename3"
done

end=`date +%h`
echo $((end-begin))

# Options
# Make the genome bam (Caution: This take 10 hours per sample)
#        		    --output-genome-bam 
# Generating a Wiggle file
#rsem-bam2wig sorted_bam_input wig_output.wig wiggle_name

# Prepare reference sequences with RSEM
#rsem-prepare-reference --gtf ./GTF/mm10knownGene.gtf \
#		       --transcript-to-gene-map ./GTF/mm10knownIsoforms.txt \
#                       ./Genome/mm10.fa \
#                       ./RSEM/mm10rsem
#
#rsem-prepare-reference --gtf ./GTF/Mus_musculus.GRCm38.73.gtf \
#                       ./Genome/GRCm38.73 \
#                       ./RSEM/NCBIM37.59
#
#rsem-prepare-reference --gtf ./GTF/Mus_musculus.NCBIM37.59.gtf \
#                       ./Genome/mm9/regular \
#                       ./RSEM/NCBIM37
#
