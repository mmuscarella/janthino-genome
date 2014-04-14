#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=50gb,walltime=10:00:00
#PBS -M mmuscare@indiana.edu
#PBS -m abe
#PBS -j oe
source ~/khmerEnv/bin/activate
module load java
module load fastqc
module load bioperl
module load velvet
module load VelvetOptimiser
PATH=$PATH:/N/soft/mason/galaxy-apps/fastx_toolkit_0.0.13
PATH=$PATH:/usr/local/share/khmer/scripts/
cd /N/dc2/projects/Lennon_Sequences/Janthino_Test
interleave-reads.py ./711_ATTCCT_L007_R1_001.fastq ./711_ATTCCT_L007_R2_001.fastq -o ./janthino.interleaved.fastq
fastq_quality_filter -Q33 -q 30 -p 50 -i ./janthino.interleaved.fastq > ./janthino.interleaved-trim.fastq
extract-paired-reads.py ./janthino.interleaved-trim.fastq
normalize-by-median.py -k 20 -C 20 -N 4 -x 2e9 -p --savehash normC20k20.kh ./janthino.interleaved-trim.fastq.pe
extract-paired-reads.py ./janthino.interleaved-trim.fastq.pe.keep 
VelvetOptimiser.pl -s 31 -e 61 -f '-shortPaired -fastq ./janthino.interleaved-trim.fastq.pe.keep.pe' -t 4 -k 'n50'
deactivate

