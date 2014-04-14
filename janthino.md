# Lennon Lab Genome Analysis: *Janthinobacter* #
---

## Project Goal: *De Novo* Genome Assembly and Annotation ##

IU has a computer cluster specifically designed for genomic analysis: [Mason](https://kb.iu.edu/data/bbhh.html)

Log-on to Mason from unix terminal

    > ssh [username]@mason.indiana.edu

Log-on to Mason from Putty (using cmd or PowerShell)

    > putty -ssh [username]@mason.indiana.edu

*OR (if you've saved settings for Mason)*

    > putty -load mason

Some of the following may take >20 minutes. 
Mason automatically kills memory/time intensive jobs.
However, you can start an interactive session using the following:

    > qsub -I -q shared -l nodes=1:ppn=4,vmem=10gb,walltime=4:00:00

Or just write the entire work flow as a script and submit with torque

***Time to get the data and start some fun***  
The original sequence files (fastq.gz) can be found at:

    /N/dc2/projects/Lennon_Sequences/Janthino_Genome

## Assembly ##
### Sequence Quality Checks ###

Before you assemble and annotate the *Janthino* genome, you first need to assess the quality of the raw data

**A. Copy files into a working directory** 

    > cd /N/dc2/projects/Lennon_Sequences/Janthino_Genome  
    > mkdir ../Janthino_Test
    > cp ./*.fastq.gz ../Janthino_Test
    > cd ../Janthino_Test
    
Run as a background process (recommended)

    > nohup cp ./*.fastq.gz ../Janthino_Test &

Did a nohub process finish?

    > ps -U[username]
    OR
    > top -U[username]

**B. Unzip compressed files** (for more information about gunzip see: [Linux / Unix Command: gzip](http://linux.about.com/od/commands/l/blcmdl1_gzip.htm))

    > cd ../Janthino_Test
    > gunzip -fv ./*.gz

**C. Check raw sequence quality with *FastQC*** (for more information see: [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

    > module load java
    > module load fastqc
    > for file in ./*fastq
    > do
    > fastqc $file >> results.out
    > done

Because the output is in HTML, you will need to move this file to your local machine before you can open it. Moving files works best in Quarry so open another ssh session

    > cd /N/dc2/projects/Lennon_Sequences/Janthino_Test/
    > cp ./*.zip /afs/iu.edu/home/[Path to RFS]

Or you can just look at the summary

    > less ./711_ATTCCT_L007_R*_001_fastqc/summary.txt
    You can switch between the two files using 
    > :n
    OR
    > :p
    Exit with
    > q

What you will notice is that there are both warnings and failures for both files. Hopefully, we can deal with these by some pre-processing.

### Sequence Pre-Processing ###

**D. Interleave Paired End Reads**  
Paired end sequencing (from HiSeq or MiSeq) yield two files per sample: R1 and R2.
To assemble the raw reads into larger contigs, aligning software needs paired reads in the same file and in the correct order. 
The process used to do this is called *interleaving*. 
There are a few tools out there for interleaving paired end sequences.
The software package *Velvet* includes a perl script: *sufflesSequences_fastq.pl*.
The software package *Khmer* includes a python script: *interleave-reads.py*.
Both should work, but I haven't actually tested this.

***Velvet Method***

    > module load velvet
    > shuffleSequences_fastq.pl  ./711_ATTCCT_L007_R1_001.fastq ./711_ATTCCT_L007_R2_001.fastq ./janthino.interleaved.fastq

The Khmer tools do not recognize these as paired reads due to changes to files names (still exploring). Use the Khmer script for interleaving.

***Khmer Methods***

Start the Khmer env (must have this downloaded prior. See: [Directions](https://khmer.readthedocs.org/en/latest/install.html)) 

    > source ~/khmerEnv/bin/activate
    > interleave-reads.py ./711_ATTCCT_L007_R1_001.fastq ./711_ATTCCT_L007_R2_001.fastq -o ./janthino.interleaved.fastq

*The Khmer method takes a long time - go get a beer!!!!*

**E. Remove Any Low Quality Reads**

If you look at the FastQC output, you will see that there is quite a bit of variation in the quality (Phred quality scores) for each file. Just as a refresher, let's run *FastQC* again:

    > fastqc ./janthino.interleaved.fastq >> results.out
    > less ./janthino.interleaved_fastqc/summary.txt

***Remove Low Quality Reads***

    > PATH=$PATH:/N/soft/mason/galaxy-apps/fastx_toolkit_0.0.13
    > fastq_quality_filter -Q33 -q 30 -p 50 -i ./janthino.interleaved.fastq > ./janthino.interleaved-trim.fastq

***Check Quality Again***

    > fastqc ./janthino.interleaved-trim.fastq >> results.out
    > less ./janthino.interleaved-trim_fastqc/summary.txt

**F. Remove Orphaned Reads**

Orphaned reads are sequences in a pair end project that have for what ever reason, lost the corresponding pair. This actually causes issues when assembling the sequences.

    > PATH=$PATH:/usr/local/share/khmer/scripts/
    > extract-paired-reads.py ./janthino.interleaved-trim.fastq
        
***Re-Check Sequence Quality with *FastQC****  
Wondering what the data look like now?

    > fastqc ./janthino.interleaved-trim.fastq.?e

**G. Digital Normalization**  
The largest issue in genome sequencing is coverage.
We want to make sure that we have good coverage across the entire genome.
However, data shows that we get unequal coverage and though the median coverage may be 50X some areas will have up to 10 times that.
To deal with this tools have been developed to normalize the data.
This benefits the assembly in multiple ways: it equalizes the coverage across the genome, it lowers the error rate, and it greatly reduces the file sizes.
All of these end up benefiting the assembly process.

**Here, we will use a digital normalization process from the *Khmer* package.**
**This method was developed by Titus Brown and Colleagues.**
**They recommend using a three step normalization method, but I've found that the best results are achieved by using only the first step.**

***Normalization***  

The following code will normalize everything to a coverage of 20; keep pairs using ‘-p’:

    > normalize-by-median.py -k 20 -C 20 -N 4 -x 2e9 -p --savehash normC20k20.kh ./janthino.interleaved-trim.fastq.pe

This produces a set of ‘.keep’ files, as well as a normC20k20.kh database file. 

**The following are steps that the Khmer protocol recommends but I find unhelpful (produce worse assemblies):**  

Use ‘filter-abund’ to trim off any k-mers that are abundance-1 in high-coverage reads. The -V option is used to make this work better for variable coverage data sets:

    > filter-abund.py -V normC20k20.kh ./*.keep

Finally, do the second round of normalization to C=5

    > normalize-by-median.py -k 20 -C 5 -N 4 -x 2e9 -p  ./*.keep.abundfilt

**If you find those steps useful, feel free to use.**
**However, you will need to adjust the following code a bit.**
**Otherwise, continue as directed.**

***Check for orphaned reads***  

The process of error trimming could have orphaned reads, so split the PE file into still-interleaved and non-interleaved reads:

    > extract-paired-reads.py ./janthino.interleaved-trim.fastq.pe.keep

**H. Genome Assembly: Using Velvet Optimization**

    > module load bioperl
    > module load velvet
    > module load VelvetOptimiser

    > VelvetOptimiser.pl -s 31 -e 61 -f '-shortPaired -fastq ./janthino.interleaved-trim.fastq.pe.keep.pe' -t 4 -k 'n50'

This process may take a while. 
It will go through all hash values from 31 to 61. 
On my trials, it will choose 57 for the final assembly. 
This might take a while, but at the end you should get less than 100 contigs. 

Results that I've gotten:
Velvet hash value: 57
Roadmap file size: 113711413
Total number of contigs: 51
n50: 766333
length of longest contig: 1316573
Total bases in contigs: 6088622
Number of contigs > 1k: 18
Total bases in contigs > 1k: 6081198
Paired Library insert stats:
Paired-end library 1 has length: 247, sample standard deviation: 85
Paired-end library 1 has length: 247, sample standard deviation: 86


## Assembly work flow as a torque script: ##

    > #!/bin/bash
    > #PBS -k o
    > #PBS -l nodes=1:ppn=8,vmem=50gb,walltime=10:00:00
    > #PBS -M mmuscare@indiana.edu
    > #PBS -m abe
    > #PBS -j oe
    > source ~/khmerEnv/bin/activate
    > module load java
    > module load fastqc
    > module load bioperl
    > module load velvet
    > module load VelvetOptimiser
    > PATH=$PATH:/N/soft/mason/galaxy-apps/fastx_toolkit_0.0.13
    > PATH=$PATH:/usr/local/share/khmer/scripts/
    > cd /N/dc2/projects/Lennon_Sequences/Janthino_Test
    > interleave-reads.py ./711_ATTCCT_L007_R1_001.fastq ./711_ATTCCT_L007_R2_001.fastq -o ./janthino.interleaved.fastq
    > fastq_quality_filter -Q 33 -q 30 -p 50 -i ./janthino.interleaved.fastq > ./janthino.interleaved-trim.fastq
    > extract-paired-reads.py ./janthino.interleaved-trim.fastq
    > normalize-by-median.py -k 20 -C 20 -N 4 -x 2e9 -p --savehash normC20k20.kh ./janthino.interleaved-trim.fastq.pe
    > extract-paired-reads.py ./janthino.interleaved-trim.fastq.pe.keep
    > VelvetOptimiser.pl -s 31 -e 61 -f '-shortPaired -fastq ./janthino.interleaved-trim.fastq.pe.keep.pe' -t 4 -k 'n50'
    > deactivate







## Gene Predication ##

###A: Using prodigal ###

###B: Using RAST ###




---

## Other code that may be of interest: ##

**Remove Adapters with with *FastX-toolkit* (This is my preferred method)**

Add the FastX-toolkit to path

    > PATH=$PATH:/N/soft/mason/galaxy-apps/fastx_toolkit_0.0.13

Remove R1 Adapters

    > fastx_clipper -v -Q 33 -a GCTCTTCCGATCT -i ./711_ATTCCT_L007_R1_001.fastq -o ./trimmed/janthino.trim.R1.fastq

Remove R2 Adapters

    > fastx_clipper -v -Q 33 -a AGATCGGAAGAGC -i ./711_ATTCCT_L007_R2_001.fastq -o ./janthino.trim.R2.fastq

**Remove Adapters with *cutadapt***

    > module load cutadapt
    > cutadapt -a GCTCTTCCGATCT ./711_ATTCCT_L007_R1_001.fastq -o ./trimmed/janthiono.trim2.R1.fastq
    > cutadapt -a AGATCGGAAGAGC ./711_ATTCCT_L007_R2_001.fastq -o ./trimmed/janthino.trim2.R2.fastq

Method | Pros | Cons 
:--------: |:-----: |:------: 
*FastX-toolkit*|Removes adapters & low quality/short seqs|Takes a while
*cutadapt* |Removes adapters and adapters with minor errors|Takes a longer time

