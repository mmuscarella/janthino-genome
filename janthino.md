# Lennon Lab Genome Analysis: *Janthinobacter* #
---

## Project Goal: *De Novo* Genome Assembly and Annotation ##

IU has a computer cluster specifically designed for genomic analysis: [Mason](https://kb.iu.edu/data/bbhh.html)

Log-on to Mason from unix terminal

    > ssh [username]@mason.indiana.edu

Log-on to Mason from Putty (using cmd or PowerShell)

    > putty -ssh [username]@mason.indiana.edu
    OR (if you've saved settings for Mason)
    > putty -load mason

The original sequence files (fastq.gz) can be found at:

    /N/dc2/projects/Lennon_Sequences/Janthino_Genome

### Sequence Quality Checks ###

Before you assemble and annotate the *Janthino* genome, you first need to assess the quality of the raw data

**A. Copy files into a working directory** 

    > cd /N/dc2/projects/Lennon_Sequences/Janthino_Genome  
    > mkdir ../janthino_test
    > cp ./*.fastq.gz ../janthino_test
    
    Run as a background process (recommended)
    > nohup cp ./*.fastq.gz ../janthino_test &

    Did nohub finish?
    > ps -U[username]

**B. Unzip compressed files** (for more information about gunzip see: [Linux / Unix Command: gzip](http://linux.about.com/od/commands/l/blcmdl1_gzip.htm))

    > cd ../janthino_test
    > gunzip -fv ./*.gz

**C. Check raw sequence quality with *FastQC*** (for more information see: [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

    > module load java
    > module load fastqc
    > for file in ./*fastq
    > do
    > fastqc $file >> results.out
    > done

### Sequence Pre-Processing ###

** D. Remove Adapters with **

    Add the FastX-toolkit to path
    > PATH=$PATH:/N/soft/mason/galaxy-apps/fastx_toolkit_0.0.13
    Remove R1 Adapters
    > fastx_clipper -v -a GCTCTTCCGATCT -i [infile] -o ./trimmed/[outfile]
    Remove R2 Adapters
    > fastx_clipper -v -a AGATCGGAAGAGC -i [infile] -o ./trimmed/[outfile]



./cutadapt -f fastq -O $stringency -q 20 -a AGATCGGAAGAGC input_file.fastq
        

## Stopping Point 3/20/2014 ##


E. Interleave Paired End Reads


R1 and R2 (paired end sequencing)
        >module load velvet
        > sufflesSequences_fastq.pl ./trimmed/"sequence_R1" ./trimmed"sequence_R2" ./interleaved/"output_filename"
        
        
        I'm not sure how to make this more automated yet, so probably just sequence by sequence
        
        There is also an interleave command in Khmer
        
/usr/local/share/khmer/scripts/interleave-reads.py s?_pe > combined.fq
Remove Low Quality Reads
> 
fastq_quality_filter -Q33 -q 30 -p 50 -i combined.fq > combined-trim.fq
Remove Orphanned Reads
        
        
    7. Re-Check Sequence Quality with fastqc (Mason Cluster)
    
    8. Here is the point where we need to do some quality filtering. 
        
    7. Use Velvet to Assemble sequences
    Mason automatically kills memory intensive jobs that aren't submitted via qsub
    However; you can start an interactive session using the following
    > qsub -I -q shared -l nodes=1:ppn=4,vmem=10gb,walltime=4:00:00
    Otherwise, use a qsub script
    > velveth auto 31,45,2 -fastq -shortPaired1 "interleaved_input_file"
    Velveth reads in these sequence files and simply produces a hashtable  and two output files (Roadmaps and Sequences) which are necessary for  the subsequent program, velvetg. 
    > velvetg auto_33 -exp_cov auto
    
Digital Normalizatoin
source khmerEnv/bin/activate
qsub -I -q shared -l nodes=1:ppn=4,vmem=10gb,walltime=4:00:00
normalize-by-median.py -x 2e9 -k 21 711_interleaved.fastq
strip-and-split-for-assembly.py ecoli_ref.fq.gz.keep.abundfilt.keep
% for i in {19..51..2}; do
    velveth ecoli.kak.$i $i -fasta -short ecoli*.se -shortPaired ecoli*.pe;
    velvetg ecoli.kak.$i -exp_cov auto -cov_cutoff auto -scaffolding no;
  done
PATH=$PATH:/N/soft/mason/galaxy-apps/fastx_toolkit_0.0.13