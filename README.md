# pacbio_clr

A reproducible Snakemake pipeline for the de-novo assembly of PacBio continuous long reads.

Steps to running on ronin
Change profile/slurm/config.yaml "phoenix" to ronin (references the cluster-configs/ronin.yaml

### Wtdbg2
git clone https://github.com/ruanjue/wtdbg2
cd wtdbg2 && make 

### Mummer v4.0
wget https://github.com/mummer4/mummer/releases/download/v4.0.0beta2/mummer-4.0.0beta2.tar.gz
tar -xf mummer-4.0.0beta2.tar.gz
rm mummer-4.0.0beta2.tar.gz
cd mummer-4.0.0beta2/
./configure
make

### LTR_FINDER_Parallel (requires downloading only)
git clone https://github.com/oushujun/LTR_FINDER_parallel.git
