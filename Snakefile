PREFIXES = ["m64015_90510_20042"]
MITOCONDRIA_REF = "reference/mitochondria.fasta",
CHLOROPLASTS_REF = "reference/chloroplast.fasta",
#FRACS = ["0.01"]
READ_COVERAGE = ["100"]
MAX_THREADS = 32
#LENGTH_CUTOFF = 5000
#LENGTH_CUTOFF_PR = 12000
BBDUK_KMER_LENGTH = ["17"]
BBDUK_MIN_COVERAGE = ["0.7"]
ECOLI_NUMS = ["1", "2", "3"]
LENGTH_CHLOROPLAST = ["134502"]
LENGTH_MITOCHONDRIA = ["415805"]
LENGTH_GENOME = ["5500000"]
ASSEMBLY_TOOLS = ["flye", "raven"]
ASSEMBLY_TYPE = ["chloroplast", "mitochondria"]

localrules: 
	all,

rule all:
	input:
		"SequelToolsResults/summaryTable.txt",
		expand("1_{ASS_TYPE}_reads/{PREFIX}_{KMER}_{COV}.fasta", ASS_TYPE = ["mitochondria","chloroplast"], PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE),
		expand("2_{ASS_TYPE}_subset/{PREFIX}_{KMER}_{COV}_{DEPTH}.fasta", ASS_TYPE = ASSEMBLY_TYPE, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE, DEPTH = READ_COVERAGE),
		expand("3_{ASS_TYPE}_assembly/{TOOL}/{PREFIX}/{KMER}_{COV}_{DEPTH}/assembly.fasta", ASS_TYPE = ASSEMBLY_TYPE, TOOL = ASSEMBLY_TOOLS, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE, DEPTH = READ_COVERAGE),

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
#rule sambamba_random:
#	input:
#		"0_raw/{prefix}.subreads.bam",
#	output:
#		"1_subset/{prefix}_{frac}.subreads.bam",
#	log:
#		"logs/sambamba_random/{prefix}_{frac}.log",
#	benchmark:
#		"benchmarks/sambamba_random/{prefix}_{frac}.tsv",
#	threads:
#		MAX_THREADS
#	conda:
#		"envs/sambamba.yaml",
#	shell:
#		"""
#		(sambamba view -h -t {threads} -f bam --subsampling-seed=123 -s {wildcards.frac} {input} -o {output}) 2> {log}
#		"""

##Convert from bam to fastq
#rule samtools_fastq:
#	input:
#		"0_raw/{prefix}.subreads.bam",
#	output:
#		"1_fastq/{prefix}.fastq",
#	log:
#		"logs/samtools_fastq/{prefix}.log",
#	benchmark:
#		"benchmarks/samtools_fastq/{prefix}.tsv",
#	threads:
#		MAX_THREADS
#	conda:
#		"envs/samtools.yaml",
#	shell:
#		"""
#		(samtools fastq -0 {output} {input} -@ {threads}) 2> {log}
#		"""

#subset bam into mitochondria, chloroplast and genome.
rule bbduk:
	input:
		seq = "0_raw/{prefix}.subreads.bam",
		ass = "reference/{ass_type}.fasta",
	output:
		out = "1_{ass_type}_reads/{prefix}_{kmer}_{cov}.fasta",
	log:
		"logs/bbduk/{ass_type}/{prefix}_{kmer}_{cov}.log",
	benchmark:
		"benchmarks/bbduk/{ass_type}/{prefix}_{kmer}_{cov}.tsv",
	wildcard_constraints:
		ass_type = "(mitochondria|chloroplast)",
	threads:
		20
	conda:
		"envs/bbmap.yaml",
	shell:
		"""
		echo "Extracting {wildcards.ass_type} reads with bbduk"
		bbduk.sh in={input.seq} outm={output.out} ref={input.ass} threads={threads} k={wildcards.kmer} -Xmx1g mincovfraction={wildcards.cov}
		"""

#rule bbduk_genome:
#	input:
#		seq = "0_raw/{prefix}.subreads.bam",
#		refs = ["reference/chloroplast.fasta", "reference/mitochondria.fasta"],
#	output:
#		"1_genome_reads/{prefix}_{kmer}_{cov}.fastq",
#	log:
#		"logs/bbduk_genome/{prefix}_{kmer}_{cov}.log",
#	benchmark:
#		"benchmarks/bbduk_genome/{prefix}_{kmer}_{cov}.tsv",
#	threads:
#		32
#	conda:
#		"envs/bbmap.yaml",
#	shell:
#		"""
#		bbduk.sh in={input.seq} out={output} ref={input.refs} 
#		"""
	
	
#subset each to specified depth for each genome
rule seqtk_mitochondria:
	input:
		"1_mitochondria_reads/{prefix}_{kmer}_{cov}.fasta",
	output:
		"2_mitochondria_subset/{prefix}_{kmer}_{cov}_{depth}.fasta",
	log:
		"logs/seqtk/mitochondria/{prefix}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/seqtk/mitochondria/{prefix}_{kmer}_{cov}_{depth}.tsv",
	threads:
		1
	params:
		length = LENGTH_MITOCHONDRIA,
	conda:
		"envs/seqtk.yaml",
	shell:
		"""
		BP_READS="$(grep -v "^>" {input} | wc | awk "{{print \$3-\$1}}")"
#		sed -n '2~4p' ../1_mitochondria_reads/m64015_90510_20042_17_0.7.fastq | awk '{{tot+=length($1)}}END{{print tot}}'
		NO_READS="$(grep '>' {input} | wc -l)"
		echo "Total number of reads: ${{NO_READS}}"
		AVG="$(( ${{BP_READS}} / ${{NO_READS}} ))"
		echo "AVG length of reads: ${{AVG}}"
		NUM="$(( {params.length} * {wildcards.depth} / ${{AVG}} ))"
		echo "Number of reads required to achieve depth of {wildcards.depth}: ${{NUM}}"
		MAX_DEPTH="$(( ${{BP_READS}} / {params.length} ))"
		if [ ${{NUM}} -gt ${{NO_READS}} ]; then
			NUM=${{NO_READS}}
			echo "Number of reads required for requested coverage is greater than the number of reads available." 
			echo "All reads will be used, giving a coverage depth of ${{MAX_DEPTH}}x"
		fi
		echo "Subsetting ..."
		seqtk sample -s100 {input} ${{NUM}} > {output}
		"""

rule seqtk_chloroplast:
	input:
		"1_chloroplast_reads/{prefix}_{kmer}_{cov}.fasta",
	output:
		"2_chloroplast_subset/{prefix}_{kmer}_{cov}_{depth}.fasta",
	log:
		"logs/seqtk/chloroplast/{prefix}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/seqtk/chloroplast/{prefix}_{kmer}_{cov}_{depth}.tsv",
	threads:
		1
	params:
		length = LENGTH_CHLOROPLAST,
	conda:
		"envs/seqtk.yaml",
	shell:
		"""
		BP_READS="$(grep -v "^>" {input} | wc | awk "{{print \$3-\$1}}")"
#		sed -n '2~4p' ../1_mitochondria_reads/m64015_90510_20042_17_0.7.fastq | awk '{{tot+=length($1)}}END{{print tot}}'
		NO_READS="$(grep '>' {input} | wc -l)"
		echo "Total number of reads: ${{NO_READS}}"
		AVG="$(( ${{BP_READS}} / ${{NO_READS}} ))"
		echo "AVG length of reads: ${{AVG}}"
		NUM="$(( {params.length} * {wildcards.depth} / ${{AVG}} ))"
		echo "Number of reads required to achieve depth of {wildcards.depth}: ${{NUM}}"
		MAX_DEPTH="$(( ${{BP_READS}} / {params.length} ))"
		if [ ${{NUM}} -gt ${{NO_READS}} ]; then
			NUM=${{NO_READS}}
			echo "Number of reads required for requested coverage is greater than the number of reads available." 
			echo "All reads will be used, giving a coverage depth of ${{MAX_DEPTH}}x"
		fi
		echo "Subsetting ..."
		seqtk sample -s100 {input} ${{NUM}} > {output}
		"""

rule seqtk_genome:
	input:
		"1_genome_reads/{prefix}_{kmer}_{cov}.fastq",
	output:
		"2_genome_subset/{prefix}_{kmer}_{cov}_{depth}.fastq",
	log:
		"logs/seqtk/genome/{prefix}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/seqtk/genome/{prefix}_{kmer}_{cov}_{depth}.tsv",
	threads:
		1
	params:
		length = LENGTH_GENOME,
	conda:
		"envs/seqtk.yaml",
	shell:
		"""
		BP_READS=$(grep -v '>' {input} | wc | awk '{{print $3-$1}}')
		NO_READS=$(grep '>' {input} | wc -l)
		AVG=$((BP_READS / NO_READS))
		NUM=$(({params.length}*100/AVG))
		if [ ${{NUM}} > ${{NO_READS}} ]; then
			NUM=${{NO_READS}}
		fi
		seqtk sample -s100 {input} ${{NUM}} > {output}
		"""


#Generate Flye assembly

rule flye_genome:
	input:
		"2_genome_subset/{prefix}_{kmer}_{cov}_{depth}.fastq",
	output:
		"3_genome_assembly/flye/{prefix}/{kmer}_{cov}_{depth}/assembly.fasta"
	log:
		"logs/flye/genome/{prefix}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/flye/genome/{prefix}_{kmer}_{cov}_{depth}.tsv",
	threads:
		MAX_THREADS
	params:
		out = "3_genome_assembly/flye/{prefix}/{kmer}_{cov}_{depth}/assembly.fasta",
		size = LENGTH_GENOME,
	conda:
		"envs/flye.yaml",
	shell:
		"""
		(flye --pacbio-raw {input} --genome-size {params.size} --out-dir {params.out} --threads {threads}) 2> {log}
		"""

rule flye:
	input:
		"2_{ass_type}_subset/{prefix}_{kmer}_{cov}_{depth}.fasta",
	output:
		"3_{ass_type}_assembly/flye/{prefix}/{kmer}_{cov}_{depth}/assembly.fasta"
	log:
		"logs/flye/{ass_type}/{prefix}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/flye/{ass_type}/{prefix}_{kmer}_{cov}_{depth}.tsv",
	threads:
		10
	params:
		out = "3_{ass_type}_assembly/flye/{prefix}/{kmer}_{cov}_{depth}",
		chlor = LENGTH_CHLOROPLAST,
		mito = LENGTH_MITOCHONDRIA,
	conda:
		"envs/flye.yaml",
	shell:
		"""
		if [ {wildcards.ass_type} == 'chloroplast' ]; then
			SIZE={params.chlor}
			echo "Starting chloroplast assembly ..."
		fi

		if [ {wildcards.ass_type} == 'mitochondria' ]; then
			SIZE={params.mito}
			echo "Starting mitochondria assembly ..."
		fi

		(flye --pacbio-raw {input} --genome-size ${{SIZE}} --out-dir {params.out} --threads {threads}) 2> {log}
		"""


##Assemble chloroplast using Raven-assembler
##only takes a single input file for long reads
rule raven:
	input:
		"2_{ass_type}_subset/{prefix}_{kmer}_{cov}_{depth}.fasta",
	output:
		fasta = "3_{ass_type}_assembly/raven/{prefix}/{kmer}_{cov}_{depth}/assembly.fasta",
		gfa = "3_{ass_type}_assembly/raven/{prefix}/{kmer}_{cov}_{depth}/assembly.gfa",
	log:
		"logs/raven/{ass_type}/{prefix}_{kmer}_{cov}_{depth}.log",	
	benchmark:
		"benchmarks/raven/{ass_type}/{prefix}_{kmer}_{cov}_{depth}.tsv",
	threads:
		10
	conda:
		"envs/raven.yaml",
	shell:
		"""
		raven -t {threads} --polishing 1 --graphical-fragment-assembly {output.gfa} {input} > {output.fasta}
		"""

#Assemble chloroplast using smartDenovo - no polishing step included
rule smartdenovo:
	input:
		"2_{ass_type}_subset/{prefix}_{kmer}_{cov}_{depth}.fasta",
	output:
		fasta = "3_{ass_type}_assembly/smartdenovo/{prefix}/{kmer}_{cov}_{depth}/assembly.fasta",
	log:
		"logs/smartdenovo/{ass_type}/{prefix}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/smartdenovo/{ass_type}/{prefix}_{kmer}_{cov}_{depth}.tsv",
	threads:
		10
	params:
		prefix = "assembly",
		location = "3_{ass_type}_assembly/smartdenovo/{prefix}/{kmer}_{cov}_{depth}",
	conda:
		"envs/smartdenovo.yaml",
	shell:	
		"""
		smartdenovo.pl -p {params.prefix} -t {threads} -c 1 {input} > {params.prefix}.mak
		make -f {params.prefix}.mak
		mv {params.prefix}* {params.location}
		"""

rule quast:
	input:
		fasta = "3_{ass_type}_assembly/{tool}/{prefix}/{kmer}_{cov}_{depth}/assembly.fasta",
	output:
		"4_quast/{ass_type}/{tool}/{sample}/{kmer}_{cov}_{depth}/report.tsv",
	log:
		"logs/quast/{ass_type}/{tool}_{sample}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/quast/{ass_type}/{tool}_{sample}_{kmer}_{cov}_{depth}.tsv",
	threads:
		1
	params:
		out = "4_quast/{tool}/{ass_type}{sample}_{kmer}_{cov}_{depth}",
	conda:
		"envs/quast.yaml",
	shell:
		"""
		(quast {input.fasta} --threads {threads} -o {params.out}) 2> {log}
		"""

#rule flye_mitochondria:
#	input:
#		"2_mitochondria_subset/{prefix}_{kmer}_{cov}_{depth}.fastq",
#	output:
#		"3_mitochondria_assembly/flye/{prefix}/{kmer}_{cov}_{depth}/assembly.fasta"
#	log:
#		"logs/flye/mitochondria/{prefix}_{kmer}_{cov}_{depth}.log",
#	benchmark:
#		"benchmarks/flye/mitochondria/{prefix}_{kmer}_{cov}_{depth}.tsv",
#	threads:
#		5
#	params:
#		out = "3_mitochondria_assembly/flye/{prefix}/{kmer}_{cov}_{depth}/assembly.fasta",
#		size = LENGTH_MITOCHONDRIA,
#	conda:
#		"envs/flye.yaml",
#	shell:
#		"""
#		(flye --pacbio-raw {input} --genome-size {params.size} --out-dir {params.out} --threads {threads}) 2> {log}
#		"""

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
#rule samtools_raw_fastq:
#	input:
#		"0_raw/{prefix}.subreads.bam",
#	output:
#		"2_fastq/all_{prefix}.fastq",
#	log:
#		"logs/samtools_raw_fastq/{prefix}_all.log",
#	benchmark:
#		"benchmarks/samtools_raw_fastq/{prefix}_all.tsv",
#	threads:
#		1
#	conda:
#		"envs/samtools.yaml",
#	shell:
#		"""
#		(samtools fastq -0 {output} {input} -@ {threads}) 2> {log}
#		"""
#
#
##bbduk to extract chloroplast reads
#rule bbduk:
#	input:
#		fa = "2_fastq/all_{prefix}.fastq",
#		ref = "references/chloroplast.fasta",
#	output:
#		match = "0_chloroplast/{prefix}_{kmer}_{cov}.subreads.fasta",
#	log:
#		"logs/bbduk/{prefix}_{kmer}_{cov}.log",
#	benchmark:
#		"benchmarks/bbduk/{prefix}_{kmer}_{cov}.tsv",
#	threads:
#		10
#	conda:
#		"envs/bbmap.yaml",
#	shell:
#		"""
#		(bbduk.sh in={input.fa} outm={output.match} ref={input.ref} threads={threads} k={wildcards.kmer} -Xmx1g mincovfraction={wildcards.cov}) 2> {log}
#		"""
#
##Assemble chloroplast using flye.
#rule flye_chloroplast:
#	input:
#		"0_chloroplast/{prefix}_{kmer}_{cov}.subreads.fasta",
#	output:
#		"3_flye/chloroplast/{prefix}/kmer_{kmer}_cov_{cov}/assembly.fasta",
#	log:
#		"logs/flye_chloroplast/{prefix}_{kmer}_{cov}.log",
#	benchmark:
#		"benchmarks/flye_chloroplast/{prefix}_{kmer}_{cov}.tsv",
#	threads:
#		5
#	params:
#		out = "3_flye/chloroplast/{prefix}/kmer_{kmer}_cov_{cov}",
#		size = LENGTH_CHLOROPLAST,  
#	conda:
#		"envs/flye.yaml",
#	shell:
#		"""
#		flye --pacbio-raw {input} --genome-size {params.size} --out-dir {params.out} --threads {threads}
#		"""	
#

#
##Compare chloroplast assemblies
#rule mummer:
#	input:
#		ref = "3_raven/chloroplast/{prefix}/kmer_{kmer}_cov_{cov}/assembly.fasta",
#		query = "3_flye/chloroplast/{prefix}/kmer_{kmer}_cov_{cov}/assembly.fasta"
#	output:
#		mums = "mummer/{prefix}_{kmer}_{cov}_raven-flye.mums",
#		gp = "mummer/{prefix}_{kmer}_{cov}_raven-flye.gp",
#	log:
#		"logs/mummer/{prefix}_{kmer}_{cov}.log",
#	benchmark:
#		"benchmarks/mummer/{prefix}_{kmer}_{cov}.log",
#	threads:
#		1
#	params:
#		pref = "{prefix}_{kmer}_{cov}_raven-flye",
#		
#	conda:
#		"envs/mummer.yaml",
#	shell:
#		"""
#		mummer -mum -b -c {input.ref} {input.query} > {output.mums}
#		mummerplot --postscript --prefix={params.pref} {output.mums}
#
#		mv {params.pref}.* mummer/
#		"""
#


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
