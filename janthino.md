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

    We also need a few folders for out outputs
    > mkdir ./trimmed
    > mkdir ./interleaved
    > mkdir ./quality
    > mkdir ./processing

    I will explain what each of these are for later.

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

**D. Remove Adapters with with *FastX-toolkit***

    Add the FastX-toolkit to path
    > PATH=$PATH:/N/soft/mason/galaxy-apps/fastx_toolkit_0.0.13
    Remove R1 Adapters
    > fastx_clipper -v -a GCTCTTCCGATCT -i ./janthino.R1.fastq -o ./trimmed/janthino.trim.R1.fastq
    Remove R2 Adapters
    > fastx_clipper -v -a AGATCGGAAGAGC -i ./janthino.R2.fastq -o ./trimmed/janthino.trim.R2.fastq

**OR Remove Adapters with *cutadapt***

    > module load cutadapt
    > cutadapt -a GCTCTTCCGATCT ./janthino.trim.R1.fastq -o ./trimmed/janthiono.trim.R1.fastq
    > cutadapt -a AGATCGGAAGAGC ./janthino.trim.R2.fastq -o ./trimmed/janthino.trim.R2.fastq

Method | Pros | Cons 
:--------: |:-----: |:------: 
*FastX-toolkit*|                     |                     
*cutadapt* |                           |                  

**E. Interleave Paired End Reads**  
Paired end sequencing (from HiSeq or MiSeq) yield two files per sample: R1 and R2.
To assemble the raw reads into larger contigs, aligning software needs paired reads in the same file and in the correct order. 
The process used to do this is called *interleaving*. There are a few tools out there for interleaving paired end sequences.
The software package *Velvet* includes a perl script: *sufflesSequences_fastq.pl*.
The software package *Khmer* includes a python script: *interleave-reads.py*.
Both should work, but I haven't actually tested this.

***Velvet Method***

    >module load velvet
    > sufflesSequences_fastq.pl ./trimmed/"sequence_R1" ./trimmed/"sequence_R2" ./interleaved/"output_filename"

***Khmer Methods***

    Add the Khmer to path
    > PATH=$PATH:/usr/local/share/khmer/scripts/
    > interleave-reads.py ./trimmed/"sequence_R1" ./trimmed/"sequence_R2" ./interleaved/"output_filename"
    From Khmer-protocols:
    > interleave-reads.py s?_pe > combined.fq

**F. Remove Any Low Quality Reads**

If you look at the FastQC output, you will see that there is quite a bit of variation in the quality (Phred quality scores) for each file. Just as a refresher, let's run *FastQC* again:

    > fastqc ./interleaved/*.fastq >> results.out

***Remove Low Quality Reads***

    > fastq_quality_filter -Q33 -q 30 -p 50 -i ./interleaved/janthino.combined.fastq > ./quality/janthino.combined-trim.fastq

***Check Quality Again

    > fastqc ./quality/janthino.combined-trim.fastq

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