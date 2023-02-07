#!/bin/bash

## Modify this job script accordingly and submit with the command:
##    sbatch HPC.sbatch
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # 16 processor core(s) per node
#SBATCH --job-name='NCDD'
#SBATCH --partition="interactive"
#SBATCH --mem=100000
#SBATCH --mail-user=weijia.su@duke.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output="NCDD-%j.out"
#SBATCH --error="NCDD-%j.err"
## LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

dir="/data/zhanglab/Weijia_Su/Nanopore_Raw_Data/20220926-human_ecDNA/"
Hg38="/data/zhanglab/Weijia_Su/Genomes/Human/hg38.fa"

#for i in 07 08 09 10 11 12;
#for i in unclassified;
#do
#reads=$dir"barcode"$i".fastq";
#reads=$dir$i".fastq";
#for ref in $Hg38;
#do
#/data/zhanglab/Weijia_Su/anaconda3/bin/time -v bash ncdd.sh $reads $ref;
#done
#done

minimap2 -ax map-ont /data/zhanglab/Weijia_Su/Genomes/Human/hg38.fa BoudaryReads_barcode09.fa -Y | samtools view -bS -F 4| samtools sort > BoudaryReads_barcode09_hg39.bam
