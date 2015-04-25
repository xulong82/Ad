#!/bin/sh
# Author: XuLong Wang (xulong.wang@jax.org)

# input: bam files from mybowtie.sh
# output: HDF5 file and expression estimation by EMASE

module load bowtie
module load python 
module load samtools
module load rsem/1.2.15

echo $0
begin=`date +%h`

REF="/data/xwang/AD/B6xC3H/B6xC3H.transcriptome"
DIR="/data/xwang/AD/emase"
SRCDIR="/data/xwang/AD/SH/empack"
IDFILE="/data/xwang/AD/B6xC3H/ENSMUST.ids"
LENDIR="/data/xwang/AD/B6xC3H"
GRPFILE="/data/xwang/AD/B6xC3H/ENSMUSG_TO_ENSMUST.tsv"
INDEX_DTYPE="uint32" 
DATA_DTYPE="uint8" 
MINUNIQ=0
PSEUDOCOUNT=0.0
MAXITERS=500
TOLERANCE=0.01

files=`find $DIR -name '*LaneALL.bam'`
for filename1 in $files; do
  filename2=`basename $filename1`
  filename3=${filename2/.bam/}
  echo "$filename3"
  echo "Creating HDF5 file"
  ${SRCDIR}/AlignmentProcessingFactories.py -a "$DIR"/"$filename3".bam \
                                    -i ${IDFILE} \
                                    -s B6,C3H \
                                    -o "$DIR"/"$filename3".h5 \
                                    -t ${INDEX_DTYPE} \
                                    -T ${DATA_DTYPE}
  echo "Expression estimation by EMASE"
   ${SRCDIR}/emase.py -i "$DIR"/"$filename3".h5 \
              -g ${GRPFILE} \
              -d ${LENDIR} \
              -b ${MINUNIQ} \
              -p ${PSEUDOCOUNT} \
              -m ${MAXITERS} \
              -t ${TOLERANCE} \
              -o "$DIR"/${filename3}.emase
# echo "Expression estimation by RSEM"
# rsem-calculate-expression -p 20 \
#  			    --phred33-quals \
#  			    --forward-prob 0.5 \
#  			    --paired-end \
#                            "$DIR1"/"$filename"_R1.fastq \
#                            "$DIR1"/"$filename"_R2.fastq \
#  			    "$REF" \
#  			    "$DIR2"/"$filename"
done

end=`date +%h`
echo $((end-begin))

