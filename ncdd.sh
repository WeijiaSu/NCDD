#!/bin/bash

read=$1
readName=$(basename $read);
NanoStat --fastq $read -n $readName".stats"
porechop -i $read --extra_end_trim 0 --discard_middle -o $readName".pre.fastq"
ref=$2
refName=$(basename $ref);
minimap2 -ax map-ont $ref $readName".pre.fastq" -Y -t 16 | samtools view -bS | samtools sort > $readName"-"$refName".bam";
samtools index $readName"-"$refName".bam";
samtools stats $readName"-"$refName".bam" | head -n 30 > $readName"-"$refName".bam.stats";

python3 /data/zhanglab/Weijia_Su/Git/NCDD/NCDD.py -bam $readName"-"$refName".bam" -ref $ref -Prefix $readName"_"$refName;
python3 /data/zhanglab/Weijia_Su/Git/NCDD/ncdd_res.py $readName"_"$refName".candi.tsv_circleAnalyze.txt" $readName"_"$refName".Type.txt";

bamCoverage -b $readName"-"$refName".bam" -o $readName"-"$refName".bam.bw";
