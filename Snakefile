PREFIXES = ["m64015_90510_20042"]
FRACS = ["0.01"]
MAX_THREADS = 32
ECOLI_NUMS = ["1", "2", "3"]

localrules:
	all

rule all:
	input:
		"SequelToolsResults/summaryTable.txt",
		expand("0_raw/{PREFIX}_{FRAC}.subreads.bam", PREFIX = PREFIXES, FRAC = FRACS),
		expand("1_fasta/{PREFIX}_{FRAC}.fasta", PREFIX = PREFIXES, FRAC = FRACS),
		"fofn.txt",

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

rule make_fofn:
	input:
		expand("1_fasta/ecoli.{num}.fasta", num = ECOLI_NUMS),
	output:
		"fofn.txt",
	log:
		"logs/make_fofn/make_fofn.log",
	benchmark:
		"benchmarks/make_fofn/make_fofn.tsv",
	threads:
		1
	shell:
		"""
		ls ./1_fasta/ecoli* > {output}
		"""

rule pb_falcon:
	input:
		file = expand("1_fasta/ecoli.{num}.fasta", num = ECOLI_NUMS),
		fofn = "fofn.txt",
	output:
		primary_contigs = "2_aligned/p_ctg.fa",
		associated_contigs = "2_aligned/a_ctg.fa",		
	log:
		"logs/pb_falcon/pb_falcon.log",
	benchmark:       
		"benchmarks/pb_falcon/pb_falcon.tsv",
	threads:
		1     
	conda:
		"envs/pacbio.yaml",
	shell:  
		"""    
		fc_run fc_run.cfg 
		"""


