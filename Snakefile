PREFIXES = ["m64015_90510_20042"]
FRACS = ["0.01"]
MAX_THREADS = 32
LENGTH_CUTOFF = 5000
LENGTH_CUTOFF_PR = 12000

ECOLI_NUMS = ["1", "2", "3"]

#Falcon is run locally as it handles its own job submission
localrules:
	all

rule all:
	input:
		"SequelToolsResults/summaryTable.txt",
		expand("0_raw/{PREFIX}_{FRAC}.subreads.bam", PREFIX = PREFIXES, FRAC = FRACS),
		expand("1_fasta/{PREFIX}_{FRAC}.fasta", PREFIX = PREFIXES, FRAC = FRACS),
		"2_assembly/assembly.fasta",

#Check read quality and length stats
rule sequelTools:
	input:
		"0_raw/bamList.txt",
	output:
		"SequelToolsResults/summaryTable.txt",
	log:
		"logs/sequelTools/sequelTools.log"
	benchmark:
		"benchmarks/sequelTools.tsv"
	threads:
		1
	conda:
		"envs/sequeltools.yaml",
	shell:
		"""
		(bash SequelTools.sh -t Q -v -u {input}) 2> {log} 
		"""

#Randmly subset the bam file
rule sambamba_random:
	input:
		"0_raw/{prefix}.subreads.bam",
	output:
		"0_raw/{prefix}_{frac}.subreads.bam",
	log:
		"logs/sambamba_random/{prefix}_{frac}.log",
	benchmark:
		"benchmarks/sambamba_random/{prefix}_{frac}.tsv",
	threads:
		MAX_THREADS
	conda:
		"envs/sambamba.yaml",
	shell:
		"""
		(sambamba view -h -t {threads} -f bam --subsampling-seed=123 -s {wildcards.frac} {input} -o {output}) 2> {log}
		"""

#Convert from bam to fasta
rule samtools_fasta:
	input:
		"0_raw/{prefix}_{frac}.subreads.bam",
	output:
		"1_fasta/{prefix}_{frac}.fasta",
	log:
		"logs/samtools_fasta/{prefix}_{frac}.log",
	benchmark:
		"benchmarks/samtools_fasta/{prefix}_{frac}.tsv",
	threads:
		1
	conda:
		"envs/samtools.yaml",
	shell:
		"""
		samtools fasta -0 {output} {input} -@ {threads}
		"""

#Generate Flye assembly
rule flye:
	input:
		file = expand("1_fasta/ecoli.{num}.fasta", num = ECOLI_NUMS),
#		file = "1_fasta/ecoli.1.fasta"
	output:
#		dir = directory("2_assembly"),
		ass = "2_assembly/assembly.fasta",
	log:
		"logs/flye/assembly.log",
	benchmark:
		"benchmarks/flye/assembly.tsv",
	threads:
		5
	params:
		out = "2_assembly"
	conda:
		"envs/flye.yaml",
	shell:
		"""
		(flye --pacbio-raw {input.file} --genome-size 5.5m --out-dir {params.out} --threads {threads}) 2> {log}
		"""


#Create the file of file names .txt for falcon
#rule make_fofn:
#	input:
#		expand("1_fasta/ecoli.{num}.fasta", num = ECOLI_NUMS),
#	output:
#		"fofn.txt",
#	log:
#		"logs/make_fofn/make_fofn.log",
#	benchmark:
#		"benchmarks/make_fofn/make_fofn.tsv",
#	threads:
#		1
#	shell:
#		"""
#		ls ./1_fasta/ecoli* > {output}
#		"""


#Generate falcon assembly
#rule pb_falcon:
#	input:
#		file = expand("1_fasta/ecoli.{num}.fasta", num = ECOLI_NUMS),
#		fofn = "fofn.txt",
#		config = "fc_run.cfg",
#	output:	
#		primary_contigs = "2_LC_{length_cut}_LCPR_{length_cut_pr}/2-asm-falcon/contig/p_ctg.fa",
#		associated_contigs = "2_LC_{length_cut}_LCPR_{length_cut_pr}/2-asm-falcon/a_ctg.fa",
#	log:
#		"logs/pb_falcon/{length_cut}_{length_cut_pr}.log",
#	benchmark:       
#		"../benchmarks/pb_falcon/{length_cut}_{length_cut_pr}.tsv",
#	threads:
#		1     
#	conda:
#		"envs/falcon.yaml",
#	shell:  
#		"""
#		mkdir -p 2_LC_{wildcards.length_cut}_LCPR_{wildcards.length_cut_pr}
#		cp -u {input.config} 2_LC_{wildcards.length_cut}_LCPR_{wildcards.length_cut_pr}/	
#		cp -u {input.fofn} 2_LC_{wildcards.length_cut}_LCPR_{wildcards.length_cut_pr}/
#		cd 2_LC_{wildcards.length_cut}_LCPR_{wildcards.length_cut_pr}		
#		(fc_run {input.config}) 2> {log} 
		"""

# Quast checks quality statistics of an assembly
#rule quast:
#	input:
#
#	output:
#
#	log:
#
#	benchmark:
#	
#	threads:
#
#	conda:
#
#	shell:
#		"""
#
#		"""

# BUSCO checks for the presence of single copy orthologs
#rule busco:
#	input:
#
#	output:
#
#	log:
#
#	benchmark:
#
#	threads:
#
#	conda:
#
#	shell:
#		"""
#
#		"""

