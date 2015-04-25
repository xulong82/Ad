#!/bin/sh
#
# Author: XuLong Wang (xulong.wang@jax.org)

module load bowtie
module load samtools

echo $0
begin=`date +%h`

index="/data/xwang/AD/B6xC3H/B6xC3H.transcriptome"
dir1="/data/xwang/AD/howell_2014/brain_trim"
dir2="/data/xwang/AD/emase"

files=`find $dir1 -name '*LaneALL_R1.fastq'`
for filename1 in $files; do
  filename2=`basename $filename1`
  filename3=${filename2/_R1.fastq/}
  echo "$filename3"
  echo "running bowtie"
  bowtie -p 20 \
         -q \
         --phred33-quals \
         --chunkmbs 512 \
         --sam \
         -a --best --strata \
         "$index" \
         -1 "$dir1"/"$filename3"_R1.fastq \
         -2 "$dir1"/"$filename3"_R2.fastq \
         "$dir2"/"$filename3".sam
  echo "running samtools"
  samtools view -bSF 4 \
                -o "$dir2"/"$filename3".bam \
                "$dir2"/"$filename3".sam
  echo "delete the sam file"
  rm "$dir2"/"$filename3".sam
done

end=`date +%h`
echo $((end-begin))

