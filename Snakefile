PREFIXES = ["m64015_90510_20042"]
MITOCONDRIA_REF = "reference/mitochondria.fasta",
CHLOROPLASTS_REF = "reference/chloroplast.fasta",
READ_COVERAGE = ["100","70","50","30","20","10"]
MAX_THREADS = 32
BBDUK_KMER_LENGTH = ["17"]
BBDUK_MIN_COVERAGE = ["0.7"]
LENGTH_CHLOROPLAST = ["134502"]
LENGTH_MITOCHONDRIA = ["415805"]
LENGTH_GENOME = ["5500000"]
ASSEMBLY_TOOLS = ["flye","raven","smartdenovo","canu"]
ASSEMBLY_TYPE = ["chloroplast"]
READ_SELECTION_METHOD = ["random","longest"]

localrules: 
	all,

rule all:
	input:
		"SequelToolsResults/summaryTable.txt",
		expand("1_{ASS_TYPE}_reads/{PREFIX}_{KMER}_{COV}.fasta", ASS_TYPE = ASSEMBLY_TYPE, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE),
		expand("1_{ASS_TYPE}_reads/sorted/{PREFIX}_{KMER}_{COV}.fasta", ASS_TYPE = ASSEMBLY_TYPE, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE),
		expand("2_{ASS_TYPE}_subset/{READ_SELECTION}/{PREFIX}_{KMER}_{COV}_{DEPTH}.fasta", ASS_TYPE = ASSEMBLY_TYPE, READ_SELECTION = READ_SELECTION_METHOD, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE, DEPTH = READ_COVERAGE),
		expand("3_{ASS_TYPE}_assembly/{TOOL}/{PREFIX}/{READ_SELECTION}_{KMER}_{COV}_{DEPTH}/assembly.fasta", ASS_TYPE = ASSEMBLY_TYPE, TOOL = ASSEMBLY_TOOLS, READ_SELECTION = READ_SELECTION_METHOD, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE, DEPTH = READ_COVERAGE),
		expand("mummer/{ASS_TYPE}_ref/{PREFIX}_{READ_SELECTION}_query_{KMERQ}_{COVQ}_{DEPTHQ}.mums", PREFIX = PREFIXES, ASS_TYPE = ASSEMBLY_TYPE, READ_SELECTION = READ_SELECTION_METHOD, KMERQ = BBDUK_KMER_LENGTH, COVQ = BBDUK_MIN_COVERAGE, DEPTHQ = READ_COVERAGE), 
		expand("quast/{ASS_TYPE}/{TOOL}/{PREFIX}/{READ_SELECTION}/{KMER}/{COV}/{DEPTH}/report.tsv", ASS_TYPE = ASSEMBLY_TYPE, TOOL = ASSEMBLY_TOOLS, PREFIX = PREFIXES, READ_SELECTION = READ_SELECTION_METHOD, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE, DEPTH = READ_COVERAGE),

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
	

rule bbmap_sort:
	input:
		"1_{ass_type}_reads/{prefix}_{kmer}_{cov}.fasta",
	output:
		"1_{ass_type}_reads/sorted/{prefix}_{kmer}_{cov}.fasta",
	log:
		"logs/bbmap_sort/{ass_type}/{prefix}_{kmer}_{cov}.log",
	benchmark:
		"benchmarks/bbmap_sort/{ass_type}/{prefix}_{kmer}_{cov}.log",
	threads:
		1
	conda:
		"envs/bbmap.yaml",
	shell:
		"""
		(sortbyname.sh in={input} out={output} name=f length=t ascending=f -Xmx1g) 2> {log}
		"""
	
#subset each to specified depth for each genome
rule select_longest:
	input:
		"1_{ass_type}_reads/sorted/{prefix}_{kmer}_{cov}.fasta",
	output:
		"2_{ass_type}_subset/longest/{prefix}_{kmer}_{cov}_{depth}.fasta",
	log:
		"logs/select_longest/{ass_type}/{prefix}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/select_longest/{ass_type}/{prefix}_{kmer}_{cov}_{depth}.tsv",
	threads:
		1
	params:
		length_mito = LENGTH_MITOCHONDRIA,
		length_chloro = LENGTH_MITOCHONDRIA,
		length_genome = LENGTH_GENOME,
	shell:
		"""
		if [ {wildcards.ass_type} == "chloroplast" ]; then
			LEN={params.length_chloro}
		elif [ {wildcards.ass_type} == "mitochondria" ]; then
			LEN={params.length_mito}
		elif [ {wildcards.ass_type} == "genome" ]; then
			LEN={params.length_genome}
		fi	
		
		BP_READS="$(grep -v "^>" {input} | wc | awk "{{print \$3-\$1}}")"
		NO_READS="$(grep '>' {input} | wc -l)"
		echo "Total number of reads: ${{NO_READS}}"
		AVG="$(( ${{BP_READS}} / ${{NO_READS}} ))"
		echo "AVG length of reads: ${{AVG}}"
		NUM="$(( ${{LEN}} * {wildcards.depth} / ${{AVG}} ))"
		echo "Number of reads required to achieve depth of {wildcards.depth}: ${{NUM}}"
		MAX_DEPTH="$(( ${{BP_READS}} / ${{LEN}} ))"
		if [ ${{NUM}} -gt ${{NO_READS}} ]; then
			NUM=${{NO_READS}}
			echo "Number of reads required for requested coverage is greater than the number of reads available." 
			echo "All reads will be used, giving a coverage depth of ${{MAX_DEPTH}}x"
		fi
		echo "Subsetting ..."
		awk "/^>/ {{n++}} n>${{NUM}} {{exit}} {{print}}" {input} > {output}
		"""

rule seqtk_mitochondria:
	input:
		"1_mitochondria_reads/{prefix}_{kmer}_{cov}.fasta",
	output:
		"2_mitochondria_subset/random/{prefix}_{kmer}_{cov}_{depth}.fasta",
	log:
		"logs/seqtk/mitochondria/{prefix}_random_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/seqtk/mitochondria/{prefix}_random_{kmer}_{cov}_{depth}.tsv",
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
		"2_chloroplast_subset/random/{prefix}_{kmer}_{cov}_{depth}.fasta",
	log:
		"logs/seqtk/chloroplast/{prefix}_random_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/seqtk/chloroplast/{prefix}_random_{kmer}_{cov}_{depth}.tsv",
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
		"2_{ass_type}_subset/{read_select}/{prefix}_{kmer}_{cov}_{depth}.fasta",
	output:
		"3_{ass_type}_assembly/flye/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta"
	log:
		"logs/flye/{ass_type}/{read_select}_{prefix}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/flye/{ass_type}/{read_select}_{prefix}_{kmer}_{cov}_{depth}.tsv",
	threads:
		10
	params:
		out = "3_{ass_type}_assembly/flye/{prefix}/{read_select}_{kmer}_{cov}_{depth}",
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

		FIRST="${{SIZE:0:1}}"
		SECOND="${{SIZE:1:1}}"

		if [ ${{#SIZE}} == "6" ]; then
		        NUSIZE="0.${{FIRST}}m"             
		elif [ ${{#SIZE}} == "7" ]; then
		        NUSIZE="${{FIRST}}.${{SECOND}}m"
		elif [ ${{#SIZE}} == "8" ]; then
		        NUSIZE="${{FIRST}}${{SECOND}}m"
		elif [ ${{#SIZE}} == "9" ]; then
		        NUSIZE="${{FIRST}}.${{SECOND}}g"
		elif ${{#SIZE}} == "10" ]; then
			NUSIZE="${{FIRST}}${{SECOND}}g"
		fi

		(flye --pacbio-raw {input} --genome-size ${{NUSIZE}} --out-dir {params.out} --threads {threads}) 2> {log}
		"""


##Assemble chloroplast using Raven-assembler
##only takes a single input file for long reads
rule raven:
	input:
		"2_{ass_type}_subset/{read_select}/{prefix}_{kmer}_{cov}_{depth}.fasta",
	output:
		fasta = "3_{ass_type}_assembly/raven/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
		gfa = "3_{ass_type}_assembly/raven/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.gfa",
	log:
		"logs/raven/{ass_type}/{read_select}_{prefix}_{kmer}_{cov}_{depth}.log",	
	benchmark:
		"benchmarks/raven/{ass_type}/{read_select}_{prefix}_{kmer}_{cov}_{depth}.tsv",
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
		"2_{ass_type}_subset/{read_select}/{prefix}_{kmer}_{cov}_{depth}.fasta",
	output:
		"3_{ass_type}_assembly/smartdenovo/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
	log:
		"logs/smartdenovo/{ass_type}/{read_select}_{prefix}_{kmer}_{cov}_{depth}.log",
	shadow:
		"shallow"
	benchmark:
		"benchmarks/smartdenovo/{ass_type}/{read_select}_{prefix}_{kmer}_{cov}_{depth}.tsv",
	threads:
		10
	params:
		prefix = "assembly",
		location = "3_{ass_type}_assembly/smartdenovo/{prefix}/{read_select}_{kmer}_{cov}_{depth}",
	conda:
		"envs/smartdenovo.yaml",
	shell:	
		"""
		smartdenovo.pl -p {params.prefix} -t {threads} -c 1 {input} > {params.prefix}.mak
		make -f {params.prefix}.mak
		mv {params.prefix}* {params.location}
		gunzip {params.location}/assembly.fa.gz
		mv {params.location}/assembly.fa {params.location}/assembly.fasta
		"""

rule canu:
	input:
		"2_{ass_type}_subset/{read_select}/{prefix}_{kmer}_{cov}_{depth}.fasta",
	output:
		"3_{ass_type}_assembly/canu/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
	log:
		"logs/canu/{ass_type}/{read_select}_{prefix}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/canu/{ass_type}_{read_select}_{kmer}_{cov}_{depth}_{prefix}.tsv"
	threads:
		1
	params:
		prefix = "assembly",
		dir = "3_{ass_type}_assembly/canu/{prefix}/{read_select}_{kmer}_{cov}_{depth}",
		chlor = LENGTH_CHLOROPLAST,
		mito = LENGTH_MITOCHONDRIA,
#	conda:
#		"envs/canu.yaml",
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
		
		canu -p {params.prefix} -d {params.dir} genomeSize=${{SIZE}} -pacbio-raw {input} useGrid=True 

		while :
		do 
			if [ ? == ? ]; then
				echo "Canu is finished!"
				break
			fi
			sleep 5m 
		done
		mv assembly.conti
		"""

rule quast:
	input:
		"3_{ass_type}_assembly/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta"
	output:
		"quast/{ass_type}/{tool}/{prefix}/{read_select}/{kmer}/{cov}/{depth}/report.tsv",
	log:
		"logs/quast/{ass_type}/{tool}/{prefix}_{read_select}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/quast/{ass_type}/{tool}/{prefix}_{read_select}_{kmer}_{cov}_{depth}.tsv",
	threads:
		1
	params:
		out = "quast/{ass_type}/{tool}/{prefix}/{read_select}/{kmer}/{cov}/{depth}",
	conda:
		"envs/quast.yaml",
	shell:
		"""
		(quast {input} --threads {threads} -o {params.out}) 2> {log}
		"""

rule mummer:
	input:
#		ref = "3_chloroplast_assembly/flye/{prefix}/{kmerr}_{covr}_{depthr}/assembly.fasta",
		chlor = "reference/chloroplast.fasta",
		mito = "reference/mitochondria.fasta",
		query = "3_{ass_type}_assembly/flye/{prefix}/{read_select}_{kmerq}_{covq}_{depthq}/assembly.fasta",
	output:
		mums = "mummer/{ass_type}_ref/{prefix}_{read_select}_query_{kmerq}_{covq}_{depthq}.mums",
		gp = "mummer/{ass_type}_ref/{prefix}_{read_select}_query_{kmerq}_{covq}_{depthq}.gp",
	log:
		"logs/mummer/{prefix}_{ass_type}/{read_select}_ref_query_{kmerq}_{covq}_{depthq}.log",
	benchmark:
		"benchmarks/mummer/{prefix}_{ass_type}/{read_select}_ref_query_{kmerq}_{covq}_{depthq}.tsv",
	threads:
		1
	shadow:
		"shallow"
	params:
		pref = "{prefix}_{read_select}_query_{kmerq}_{covq}_{depthq}",
	conda:
		"envs/mummer.yaml",
	shell:
		"""
		if [ {wildcards.ass_type} == "chloroplast" ]; then
		
			mummer -mum -b -c {input.chlor} {input.query} > {params.pref}.mums
			mummerplot --postscript --prefix={params.pref} {params.pref}.mums
			mv {params.pref}.* mummer/{wildcards.ass_type}_ref
		fi

		if [ {wildcards.ass_type} == "mitochondria" ]; then
			mummer -mum -b -c {input.mito} {input.query} > {params.pref}.mums
			mummerplot --postscript --prefix={params.pref} {params.pref}.mums
			mv {params.pref}.* mummer/{wildcards.ass_type}_ref
		fi
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
