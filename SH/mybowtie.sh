#!/bin/sh
#
# ChIP-seq data analysis pipeline
# Author: XuLong Wang (xulong.wang@jax.org)

# Usage: sh mybowtie.sh > mylog

echo $0
echo "begin: `date`"

myindex="/data/xwang/Genome/GRCm38.73/GRCm38.73"
mydir="/data/xwang/Prdm9/Encode"

for filename in H3k04me1MAdult8wksC57bl6StdRawDataRep1 \
                H3k04me1MAdult8wksC57bl6StdRawDataRep2 \
		InputMAdult8wksC57bl6StdRawDataRep1 \
		InputMAdult8wksC57bl6StdRawDataRep2; do
  echo "Processing $filename..."
  bowtie -p 32 \
         -q \
	 -I 36 \
         --phred64-quals \
	 --best \
	 --chunkmbs 512 \
	 --sam \
	 "$myindex" \
	 "$mydir"/wgEncodeLicrHistoneTestis"$filename".fastq \
	 "$mydir"/wgEncodeLicrHistoneTestis"$filename".sam
done

echo "end: `date`"

# Build index with bowtie-build
# cd /hpcdata/xwang/Genome/GRCm38.73
# bowtie-build Mus_musculus.GRCm38.73.dna.chromosome.1.fa,Mus_musculus.GRCm38.73.dna.chromosome.2.fa,Mus_musculus.GRCm38.73.dna.chromosome.3.fa,Mus_musculus.GRCm38.73.dna.chromosome.4.fa,Mus_musculus.GRCm38.73.dna.chromosome.5.fa,Mus_musculus.GRCm38.73.dna.chromosome.6.fa,Mus_musculus.GRCm38.73.dna.chromosome.7.fa,Mus_musculus.GRCm38.73.dna.chromosome.8.fa,Mus_musculus.GRCm38.73.dna.chromosome.9.fa,Mus_musculus.GRCm38.73.dna.chromosome.10.fa,Mus_musculus.GRCm38.73.dna.chromosome.11.fa,Mus_musculus.GRCm38.73.dna.chromosome.12.fa,Mus_musculus.GRCm38.73.dna.chromosome.13.fa,Mus_musculus.GRCm38.73.dna.chromosome.14.fa,Mus_musculus.GRCm38.73.dna.chromosome.15.fa,Mus_musculus.GRCm38.73.dna.chromosome.16.fa,Mus_musculus.GRCm38.73.dna.chromosome.17.fa,Mus_musculus.GRCm38.73.dna.chromosome.18.fa,Mus_musculus.GRCm38.73.dna.chromosome.19.fa,Mus_musculus.GRCm38.73.dna.chromosome.X.fa,Mus_musculus.GRCm38.73.dna.chromosome.Y.fa,Mus_musculus.GRCm38.73.dna.chromosome.MT.fa GRCm38.73

# echo "." | mail -s "Task in Rockhopper have completed!" emailofx@gmail.com 

