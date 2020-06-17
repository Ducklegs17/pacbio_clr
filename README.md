# PacBio de novo CLR assembly and QA workflow

A reproducible Snakemake pipeline for the de-novo assembly and quality assessment of PacBio continuous long reads.

## Getting Set Up

This workflow is set up for execution on a Linux HPC with SLURM scheduling. It also requires an installation of Snakemake, Conda and Singularity. 
To run, follow these steps:

1. Clone the repo
2. Install tools as detailed below
3. Create a folder called `reference` in the base directory and save the reference genome of interest in there under the name `genome.fasta`. Repeat for mitochondria and chloroplast references, naming them `mitochondria.fasta` and `chloroplast.fasta`. 
4. Create a folder called `0_raw` that contains the raw CLR reads in .subreads.bam format. 
5. Enter the prefix of this file into the `PREFIXES` variable at the very top of the Snakefile.
6. Activate Snakemake
7. Enter the desired file path (found in the `rule all:` section) with wildcards substituted for actual values 
8. Multiple values can be specified for each wildcard

eg. To create a 15x and 20x Canu assembly using bbduk parameters of 17 for kmer length and 0.7 for kmer read coverage, the following command would suffice
`3_genome_assembly/canu/m64015_90510_20042/random_17_0.7_{15,20}/assembly.fasta`

The workflow is set up to assemble nuclear genomes, chloroplast and mitochondria so `genome` must be specified as the {ASS_TYPE} wildcard

## Tools requiring installation

### Wtdbg2
```
git clone https://github.com/ruanjue/wtdbg2
cd wtdbg2 && make 
```
### Mummer v4.0
```
wget https://github.com/mummer4/mummer/releases/download/v4.0.0beta2/mummer-4.0.0beta2.tar.gz
tar -xf mummer-4.0.0beta2.tar.gz
rm mummer-4.0.0beta2.tar.gz
cd mummer-4.0.0beta2/
./configure
make
```
### LTR_FINDER_Parallel (requires downloading only)
```
git clone https://github.com/oushujun/LTR_FINDER_parallel.git
```

### Filtlong (only requires installation if you want to use the additional scripts)
```
git clone https://github.com/rrwick/Filtlong.git
cd Filtlong
make -j
bin/filtlong -h
```

