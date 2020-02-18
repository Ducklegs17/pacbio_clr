PREFIXES = ["m64015_90510_20042"]
FRACS = ["0.01"]
MAX_THREADS = 32
LENGTH_CUTOFF = 5000
LENGTH_CUTOFF_PR = 12000
#BUSCO_TAX = [""]
BBDUK_KMER_LENGTH = ["17"]
BBDUK_MIN_COVERAGE = ["0.7"]
ECOLI_NUMS = ["1", "2", "3"]
LENGTH_CHLOROPLAST = ["134502"]
LENGTH_GENOME = ["5.5m"]

localrules: 
	all

rule all:
	input:
		"SequelToolsResults/summaryTable.txt",
		expand("1_subset/{PREFIX}_{FRAC}.subreads.bam", PREFIX = PREFIXES, FRAC = FRACS),
		expand("2_fastq/frac_{PREFIX}_{FRAC}.fastq", PREFIX = PREFIXES, FRAC = FRACS),
		"3_assembly/assembly.fasta",
		"4_quast/quast_test/quast_done",
		expand("0_chloroplast/{PREFIX}_{KMER}_{COV}.subreads.fasta", PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE),
		expand("2_fastq/all_{PREFIX}.fastq", PREFIX = PREFIXES),
		expand("3_flye/chloroplast/{PREFIX}/kmer_{KMER}_cov_{COV}/assembly.fasta", PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE),
		expand("3_raven/chloroplast/{PREFIX}/kmer_{KMER}_cov_{COV}/assembly.fasta", PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE),
		expand("3_raven/chloroplast/{PREFIX}/kmer_{KMER}_cov_{COV}/assembly.gfa", PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE),

#Check read quality and length stats
rule sequelTools:
	input:
		"0_raw/bamList.txt",
	output:
		"SequelToolsResults/summaryTable.txt",
	log:
		"logs/sequelTools/sequelTools.log",
	benchmark:
		"benchmarks/sequelTools.tsv",
	threads:
		1
	conda:
		"envs/sequeltools.yaml",
	shell:
		"""
		(bash SequelTools.sh -t Q -v -u {input}) 2> {log} 
		"""

#Randomly subset the bam file
rule sambamba_random:
	input:
		"0_raw/{prefix}.subreads.bam",
	output:
		"1_subset/{prefix}_{frac}.subreads.bam",
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

#Convert from bam to fastq
rule samtools_fastq:
	input:
		"1_subset/{prefix}_{frac}.subreads.bam",
	output:
		"2_fastq/frac_{prefix}_{frac}.fastq",
	log:
		"logs/samtools_fastq/{prefix}_{frac}.log",
	benchmark:
		"benchmarks/samtools_fastq/{prefix}_{frac}.tsv",
	threads:
		1
	conda:
		"envs/samtools.yaml",
	shell:
		"""
		(samtools fastq -0 {output} {input} -@ {threads}) 2> {log}
		"""

#Generate Flye assembly
rule flye:
	input:
		file = expand("2_fastq/ecoli.{num}.fastq", num = ECOLI_NUMS),
	output:
		ass = "3_assembly/assembly.fasta",
	log:
		"logs/flye/assembly.log",
	benchmark:
		"benchmarks/flye/assembly.tsv",
	threads:
		5
	params:
		out = "3_assembly", size = LENGTH_GENOME,
	conda:
		"envs/flye.yaml",
	shell:
		"""
		(flye --pacbio-raw {input.file} --genome-size {params.size} --out-dir {params.out} --threads {threads}) 2> {log}
		"""

rule quast:
	input:
		fasta = "3_assembly/assembly.fasta",
	output:
		"4_quast/quast_test/quast_done",
	log:
		"logs/quast/assembly.log",
	benchmark:
		"benchmarks/quast/assembly.tsv",
	threads:
		1
	params:
		out = "4_quast/quast_test"
	conda:
		"envs/quast.yaml",
	shell:
		"""
		(quast {input.fasta} --threads {threads} -o {params.out}) 2> {log}
		touch {output} 
		"""


#BUSCO checks for the presence of single copy orthologs
#rule busco:
#	input:
#		fasta = "2_assembly/assembly.fasta",
#	output:
#		"5_busco/missing_busco_list.tsv",
#	log:
#		"logs/busco/ass.log",
#	benchmark:
#		"benchmarks/busco/ass.tsv",
#	threads:
#		1
#	params:
#		out = "4_busco",
#		lineage = BUSCO_LINEAGE,
#	conda:
#		"envs/busco.yaml",
#	shell:
#		"""
#		busco -m genome -i {input} -o {params.out} -l {params.lineage} -c {threads} 
#		"""

#=====================================================
# Chloroplast assembly
#=====================================================

#convert raw bam files to fastq
rule samtools_raw_fastq:
	input:
		"0_raw/{prefix}.subreads.bam",
	output:
		"2_fastq/all_{prefix}.fastq",
	log:
		"logs/samtools_raw_fastq/{prefix}_all.log",
	benchmark:
		"benchmarks/samtools_raw_fastq/{prefix}_all.tsv",
	threads:
		1
	conda:
		"envs/samtools.yaml",
	shell:
		"""
		(samtools fastq -0 {output} {input} -@ {threads}) 2> {log}
		"""


#bbduk to extract chloroplast reads
rule bbduk:
	input:
		fa = "2_fastq/all_{prefix}.fastq",
		ref = "references/chloroplast.fasta",
	output:
		match = "0_chloroplast/{prefix}_{kmer}_{cov}.subreads.fasta",
	log:
		"logs/bbduk/{prefix}_{kmer}_{cov}.log",
	benchmark:
		"benchmarks/bbduk/{prefix}_{kmer}_{cov}.tsv",
	threads:
		10
	conda:
		"envs/bbmap.yaml",
	shell:
		"""
		(bbduk.sh in={input.fa} outm={output.match} ref={input.ref} threads={threads} k={wildcards.kmer} -Xmx1g mincovfraction={wildcards.cov}) 2> {log}
		"""

#Assemble chloroplast using flye.
rule flye_chloroplast:
	input:
		"0_chloroplast/{prefix}_{kmer}_{cov}.subreads.fasta",
	output:
		"3_flye/chloroplast/{prefix}/kmer_{kmer}_cov_{cov}/assembly.fasta",
	log:
		"logs/flye_chloroplast/{prefix}_{kmer}_{cov}.log",
	benchmark:
		"benchmarks/flye_chloroplast/{prefix}_{kmer}_{cov}.tsv",
	threads:
		5
	params:
		out = "3_flye/chloroplast/{prefix}/kmer_{kmer}_cov_{cov}",
		size = LENGTH_CHLOROPLAST,  
	conda:
		"envs/flye.yaml",
	shell:
		"""
		flye --pacbio-raw {input} --genome-size {params.size} --out-dir {params.out} --threads {threads}
		"""	

#Assemble chloroplast using Raven-assembler
#only takes a single input file for long reads
rule raven_chloroplast:
	input:
		"0_chloroplast/{prefix}_{kmer}_{cov}.subreads.fasta",
	output:
		fasta = "3_raven/chloroplast/{prefix}/kmer_{kmer}_cov_{cov}/assembly.fasta",
		gfa = "3_raven/chloroplast/{prefix}/kmer_{kmer}_cov_{cov}/assembly.gfa",
	log:
		"logs/raven_chloroplast/{prefix}_{kmer}_{cov}.log",	
	benchmark:
		"benchmarks/raven_chloroplast/{prefix}_{kmer}_{cov}.tsv",
	threads:
		MAX_THREADS
	conda:
		"envs/raven.yaml",
	shell:
		"""
		raven -t {threads} --polishing 1 --graphical-fragment-assembly {output.gfa} {input} > {output.fasta}
		"""

#Assemble chloroplast using canu? or miniasm?

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
#		"""
