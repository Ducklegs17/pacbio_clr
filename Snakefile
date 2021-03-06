#!/bin/bash
#Instructions:

PREFIXES = ["m64015_90510_20042"]
MITOCONDRIA_REF = "reference/mitochondria.fasta"
CHLOROPLASTS_REF = "reference/chloroplast.fasta",
READ_COVERAGE = ["100"]
MAX_THREADS = 48
BBDUK_KMER_LENGTH = ["17"]
BBDUK_MIN_COVERAGE = ["0.7"]
LENGTH_CHLOROPLAST = ["134502"]
LENGTH_MITOCHONDRIA = ["415805"]
LENGTH_GENOME = ["387500000"]
ASSEMBLY_TOOLS = ["raven","wtdbg2","canu","flye"]
ASSEMBLY_TYPE = ["genome"]
READ_SELECTION_METHOD = ["longest","random"]
POLISH_DEPTH = ["10","25","50","75","100","125","150"]
POLISH_READ_SELECT = ["10","25","50","75","100","125","150"]
CANU_BENCHMARK_FILES = ["job_finish_state","submit_time","start_time","end_time","walltime_reserved","walltime_elapsed","max_memory","max_disk_write","max_disk_read","num_cores"]
POLISHROUND = ["1","2"]

WTDBG2_PATH1 = "/home/a1761942/fast_dir/tools/wtdbg2/wtdbg2"
WTDBG2_PATH2 = "/home/a1761942/fast_dir/tools/wtdbg2/wtpoa-cns"
MUMMER_PATH = "/home/a1761942/fast_dir/tools/mummer-4.0.0beta2/"
LTR_FINDER_PATH = "/home/a1761942/fast_dir/tools/LTR_FINDER_parallel-master/LTR_FINDER_parallel"

localrules: 
	all,
	canu,
	generate_coverage_list,
	benchcanu,
	select_longest,
	filtlong,
	contig_length,

rule all:
	input:
		"SequelToolsResults/summaryTable.txt",
		expand("0_raw/fastq/{PREFIX}.fastq", PREFIX = PREFIXES),
		expand("1_{ASS_TYPE}_reads/{PREFIX}_{KMER}_{COV}.fasta", ASS_TYPE = ASSEMBLY_TYPE, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE),
		expand("1_{ASS_TYPE}_reads/sorted/{PREFIX}_{KMER}_{COV}.fasta", ASS_TYPE = ASSEMBLY_TYPE, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE),
		expand("1_{ASS_TYPE}_reads/sorted/{PREFIX}_{KMER}_{COV}_coveragetable.txt", ASS_TYPE = ASSEMBLY_TYPE, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE),
		expand("2_{ASS_TYPE}_subset/{READ_SELECTION}/{PREFIX}_{KMER}_{COV}_{DEPTH}.fasta", ASS_TYPE = ASSEMBLY_TYPE, READ_SELECTION = READ_SELECTION_METHOD, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE, DEPTH = READ_COVERAGE),
		expand("3_{ASS_TYPE}_assembly/{TOOL}/{PREFIX}/{READ_SELECTION}_{KMER}_{COV}_{DEPTH}/assembly.fasta", ASS_TYPE = ASSEMBLY_TYPE, TOOL = ASSEMBLY_TOOLS, READ_SELECTION = READ_SELECTION_METHOD, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE, DEPTH = READ_COVERAGE),
		expand("mummer/{ASS_TYPE}/prefix_{PREFIX}_assemblytool_{TOOL}_readselect_{READ_SELECTION}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}.delta", PREFIX = PREFIXES, ASS_TYPE = ASSEMBLY_TYPE, TOOL = ASSEMBLY_TOOLS, READ_SELECTION = READ_SELECTION_METHOD, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE, DEPTH = READ_COVERAGE), 
		expand("quast/assemblytype_{ASS_TYPE}_assemblytool_{TOOL}_prefix_{PREFIX}_readselect_{READ_SELECTION}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}/report.tsv", ASS_TYPE = ASSEMBLY_TYPE, TOOL = ASSEMBLY_TOOLS, PREFIX = PREFIXES, READ_SELECTION = READ_SELECTION_METHOD, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE, DEPTH = READ_COVERAGE),
		expand("quastref/assemblytype_{ASS_TYPE}_assemblytool_{TOOL}_prefix_{PREFIX}_readselect_{READ_SELECTION}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}/report.tsv", ASS_TYPE = ASSEMBLY_TYPE, TOOL = ASSEMBLY_TOOLS, PREFIX = PREFIXES, READ_SELECTION = READ_SELECTION_METHOD, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE, DEPTH = READ_COVERAGE),
		expand("ltr/harvest/assemblytype_genome_assemblytool_{TOOL}_readselect_{READ_SELECT}_prefix_{PREFIX}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}/assembly.fa.harvest.scn",TOOL = ASSEMBLY_TOOLS, READ_SELECT = READ_SELECTION_METHOD, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE, DEPTH = READ_COVERAGE),
		expand("ltr/finder/assemblytype_genome_assemblytool_{TOOL}_readselect_{READ_SELECT}_prefix_{PREFIX}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}/assembly.fasta.finder.combine.scn", TOOL = ASSEMBLY_TOOLS, READ_SELECT = READ_SELECTION_METHOD, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE, DEPTH = READ_COVERAGE),
		expand("ltr/retriever/assemblytype_genome_assemblytool_{TOOL}_readselect_{READ_SELECT}_prefix_{PREFIX}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}/LAI_score.txt", TOOL = ASSEMBLY_TOOLS, READ_SELECT = READ_SELECTION_METHOD, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE, DEPTH = READ_COVERAGE),	
		expand("benchcanu/assemblytype_{ASS_TYPE}_assemblytool_{TOOL}_prefix_{PREFIX}_readselect_{READ_SELECT}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}/{CANU_BENCHMARKS}.txt", ASS_TYPE = ASSEMBLY_TYPE, TOOL = ASSEMBLY_TOOLS, READ_SELECT = READ_SELECTION_METHOD, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE, DEPTH = READ_COVERAGE, CANU_BENCHMARKS = CANU_BENCHMARK_FILES),
		expand("busco/assemblytype_genome_assemblytool_{TOOL}_readselect_{READ_SELECT}_prefix_{PREFIX}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}/results/run_poales_odb10/full_table.tsv", TOOL = ASSEMBLY_TOOLS, READ_SELECT = READ_SELECTION_METHOD, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE, DEPTH = READ_COVERAGE),
		expand("4_{ASS_TYPE}_polished/{TOOL}/{PREFIX}/{READ_SELECTION}_{KMER}_{COV}_{DEPTH}/{POLISHDEPTH}_{POLISH_SELECT}/{ROUND}_assembly.fasta", ASS_TYPE = ASSEMBLY_TYPE, POLISH_SELECT = POLISH_READ_SELECT,TOOL = ASSEMBLY_TOOLS, READ_SELECTION = READ_SELECTION_METHOD, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE, DEPTH = READ_COVERAGE,POLISHDEPTH = POLISH_DEPTH, ROUND = POLISHROUND),
		expand("length/{ASS_TYPE}_{TOOL}_{PREFIX}_{READ_SELECT}_{KMER}_{COV}_{DEPTH}.txt", ASS_TYPE = ASSEMBLY_TYPE, PREFIX = PREFIXES, TOOL = ASSEMBLY_TOOLS, READ_SELECT = READ_SELECTION_METHOD, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE, DEPTH = READ_COVERAGE),
		expand("quastref/assemblytype_{ASS_TYPE}_assemblytool_{TOOL}_prefix_{PREFIX}_readselect_{READ_SELECTION}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}_polishselect_{POLISH_SELECT}_polishdepth_{POLISHDEPTH}_polishround_{ROUND}/report.tsv",ASS_TYPE = ASSEMBLY_TYPE, POLISH_SELECT = POLISH_READ_SELECT,TOOL = ASSEMBLY_TOOLS, READ_SELECTION = READ_SELECTION_METHOD, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE, DEPTH = READ_COVERAGE, POLISHDEPTH = POLISH_DEPTH, ROUND = POLISHROUND),
		expand("mummer/{ASS_TYPE}/prefix_{PREFIX}_assemblytool_{TOOL}_readselect_{READ_SELECTION}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}_polishselect_{POLISH_SELECT}_polishdepth_{POLISHDEPTH}_polishround_{ROUND}.delta", ASS_TYPE = ASSEMBLY_TYPE, POLISH_SELECT = POLISH_READ_SELECT,TOOL = ASSEMBLY_TOOLS, READ_SELECTION = READ_SELECTION_METHOD, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE, DEPTH = READ_COVERAGE, POLISHDEPTH = POLISH_DEPTH, ROUND = POLISHROUND),
		expand("ltr/retriever/assemblytype_genome_assemblytool_{TOOL}_readselect_{READ_SELECT}_prefix_{PREFIX}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}_polishselect_{POLISH_SELECT}_polishdepth_{POLISHDEPTH}_polishround_{ROUND}/LAI_score.txt", ASS_TYPE = ASSEMBLY_TYPE, POLISH_SELECT = POLISH_READ_SELECT,TOOL = ASSEMBLY_TOOLS, READ_SELECT = READ_SELECTION_METHOD, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE, DEPTH = READ_COVERAGE, POLISHDEPTH = POLISH_DEPTH, ROUND = POLISHROUND),
		expand("busco/assemblytype_genome_assemblytool_{TOOL}_readselect_{READ_SELECT}_prefix_{PREFIX}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}_polishselect_{POLISH_SELECT}_polishdepth_{POLISHDEPTH}_polishround_{ROUND}/results/run_poales_odb10/full_table.tsv", ASS_TYPE = ASSEMBLY_TYPE, POLISH_SELECT = POLISH_READ_SELECT,TOOL = ASSEMBLY_TOOLS, READ_SELECT = READ_SELECTION_METHOD, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE, DEPTH = READ_COVERAGE, POLISHDEPTH = POLISH_DEPTH, ROUND = POLISHROUND),

rule get_raw_length_distribution:
	input:
		"0_raw/{prefix}.subreads.bam",
	output:
		"0_raw/{prefix}.length.txt",
	threads:
		1
	conda:
		"envs/bbmap.yaml",
	shell:
		"""
		reformat.sh in={input} out=temp_{wildcards.prefix}.fa
		cat temp_{wildcards.prefix}.fa | awk '$0 ~ ">" {{if (NR > 1) {{print c;}} c=0;printf substr($0,2,100) "\t"; }} $0 !~ ">" {{c+=length($0);}} END {{ print c; }}' > {output}
		"""

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

#subset bam into mitochondria and chloroplast
rule bbduk:
	input:
		seq = "0_raw/{prefix}.subreads.bam",
		ass = "reference/{ass_type}.fasta",
	output:
		out = "1_{ass_type}_reads/{prefix}_{kmer}_{cov}.fasta",
	log:
		"logs/bbduk/{ass_type}/{prefix}_{kmer}_{cov}.log",
	benchmark:
		"benchmarks/bbduk/assemblytype_{ass_type}_prefix_{prefix}_kmer_{kmer}_cov_{cov}.tsv",
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


rule bam2fastq:
	input:
		 seq = "0_raw/{prefix}.subreads.bam",
	output:
		"0_raw/fastq/{prefix}.fastq",
	log:
		"logs/bbduk2fastq/{prefix}.log",
	benchmark:
		"benchmarks/bbduk2fastq/{prefix}.tsv",
	threads:
		MAX_THREADS
	conda:
		"envs/bbmap.yaml",
	shell:
		"""
		reformat.sh in={input.seq} out={output} aqhist=test_meanQual_histogram.txt
		"""

rule filtlong:
	input:
		"0_raw/fastq/{prefix}.fastq",
	output:
		"2_genome_subset/filtlong/{prefix}_{kmer}_{cov}_{depth}.fasta",
	log:
		"logs/filtlong/genome_{prefix}_filtlong_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/filtlong/assemblytype_genome_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		MAX_THREADS
	params:
		length = LENGTH_GENOME,
	conda:
		"envs/filtlong.yaml",
	shell:
		"""
		DEPTH="$(( {wildcards.depth} * {params.length} ))"
		filtlong --min_length {wildcards.kmer} --min_mean_q {wildcards.cov} --target_bases ${{DEPTH}} {input} | sed -n '1~4s/^@/>/p;2~4p' > {output}
#		../tools/Filtlong/scripts/read_info_histograms.sh {input} > {output}
		"""

rule filtlong_genome:
	input:
		"0_raw/fastq/{prefix}.fastq",
	output:
		"2_genome_subset/filtlong/{prefix}_{kmer}_{cov}_{depth}.fasta",
	log:
		"logs/filtlong/genome_{prefix}_filtlong_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/filtlong/assemblytype_genome_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		MAX_THREADS
	params:
		length = LENGTH_GENOME,
	conda:
		"envs/filtlong.yaml",
	shell:
		"""
#		DEPTH="$(( {wildcards.depth} * {params.length} ))"
#		filtlong --min_length {wildcards.kmer} --min_mean_q {wildcards.cov} --target_bases ${{DEPTH}} {input} | sed -n '1~4s/^@/>/p;2~4p' > {output}
		../tools/Filtlong/scripts/read_info_histograms.sh {input} > {output}
		"""

rule bbduk_genome:
	input:
		seq = "0_raw/{prefix}.subreads.bam",
		ref1 = "reference/chloroplast.fasta", 
		ref2 = "reference/mitochondria.fasta",
	output:
		"1_genome_reads/{prefix}_{kmer}_{cov}.fasta",
	log:
		"logs/bbdukgenome/assembly_type_genome_prefix_{prefix}_kmer_{kmer}_cov_{cov}.log",
	benchmark:
		"benchmarks/bbdukgenome/assemblytype_genome_prefix_{prefix}_kmer_{kmer}_cov_{cov}.tsv",
	threads:
		MAX_THREADS
	conda:
		"envs/bbmap.yaml",
	shell:
		"""
		bbduk.sh in={input.seq} out={output} ref={input.ref1},{input.ref2} threads={threads} k={wildcards.kmer} -Xmx2g mincovfraction={wildcards.cov}
		"""

rule bbmap_sort:
	input:
		"1_{ass_type}_reads/{prefix}_{kmer}_{cov}.fasta",
	output:
		fasta = "1_{ass_type}_reads/sorted/{prefix}_{kmer}_{cov}.fasta",
	log:
		"logs/bbmap_sort/{ass_type}_{prefix}_{kmer}_{cov}.log",
	benchmark:
		"benchmarks/bbmapsort/assemblytype_{ass_type}_prefix_{prefix}_kmer_{kmer}_cov_{cov}.tsv",
	threads:
		1
	wildcard_constraints:
		ass_type = "(mitochondria|chloroplast|genome)",
#	shadow:
#		"shallow"
	conda:
		"envs/bbmap.yaml",
	shell:
		"""
#increased Xmx option from 1g to 10g for sorting genome. also increased time from 1 min to 1 hr.
		(sortbyname.sh in={input} out={output.fasta} name=f length=t ascending=f -Xmx60g) 2> {log}
		"""

rule generate_coverage_list:
	input:
		"1_{ass_type}_reads/sorted/{prefix}_{kmer}_{cov}.fasta",	
	output:
		"1_{ass_type}_reads/sorted/{prefix}_{kmer}_{cov}_coveragetable.txt",
	log:
		"logs/generate_coverage_list/assemblytype_{ass_type}_prefix_{prefix}_kmer_{kmer}_cov_{cov}.log",	
	benchmark:
		"benchmarks/generatecoveragelist/assemblytype_{ass_type}_prefix_{prefix}_kmer_{kmer}_cov_{cov}.tsv",
	threads:
		1
	params:
		length_mito = LENGTH_MITOCHONDRIA,
		length_chloro = LENGTH_CHLOROPLAST,
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
		touch length.tmp && rm length.tmp
		x=0
		awk '/^>/{{if (l!="") print l; l=0; next}}{{l+=length($0)}}' {input} >> length.tmp
		
		while IFS= read -r line
		do
			x=$(( ${{x}} + ${{line}} ))
			pr=$(( ${{x}} / ${{LEN}} ))
			echo "${{pr}}" >> {output}
		done < "length.tmp"

		rm length.tmp
		"""
	
#subset each to specified depth for each genome
rule select_longest:
	input:
		fa = "1_{ass_type}_reads/sorted/{prefix}_{kmer}_{cov}.fasta",
		list = "1_{ass_type}_reads/sorted/{prefix}_{kmer}_{cov}_coveragetable.txt",
	output:
		"2_{ass_type}_subset/longest/{prefix}_{kmer}_{cov}_{depth}.fasta",
	log:
		"logs/select_longest/{ass_type}/{prefix}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/selectlongest/assemblytype_{ass_type}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		1
	params:
		length_mito = LENGTH_MITOCHONDRIA,
		length_chloro = LENGTH_CHLOROPLAST,
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
		
		NO_READS=$(grep '^>' {input.fa} | wc -l)
		NUM=$(awk '$1<{wildcards.depth}{{c++}} END{{print c+0}}' {input.list})
		NUM=$(( ${{NUM}} + 1 ))
		
		if [ ${{NUM}} -eq ${{NO_READS}} ]; then
			DEPTH=$( tail -n 1 {input.list})
			echo "Number of reads required for requested coverage is greater than the number of reads available." 
			echo "All reads will be used, giving a coverage depth of ${{DEPTH}}x"
		fi
		
		awk "/^>/ {{n++}} n>${{NUM}} {{exit}} {{print}}" {input} > {output}
		
		"""

rule seqtk:
	input:
		"1_{ass_type}_reads/{prefix}_{kmer}_{cov}.fasta",
	output:
		"2_{ass_type}_subset/random/{prefix}_{kmer}_{cov}_{depth}.fasta",
	log:
		"logs/seqtk/{ass_type}_{prefix}_random_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/seqtk/assemblytype_{ass_type}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		1
	params:
		chlor = LENGTH_CHLOROPLAST,
		mito = LENGTH_MITOCHONDRIA,
	wildcard_constraints:
		ass_type = "(mitochondria|chloroplast)",
	conda:
		"envs/seqtk.yaml",
	shell:
		"""
		if [ {wildcards.ass_type} == "chloroplast" ]; then
			LENGTH={params.chlor}
		fi

		if [ {wildcards.ass_type} == "mitochondria" ]; then
			LENGTH={params.mito}
		fi


		BP_READS="$(grep -v "^>" {input} | wc | awk "{{print \$3-\$1}}")"
#		sed -n '2~4p' ../1_mitochondria_reads/m64015_90510_20042_17_0.7.fastq | awk '{{tot+=length($1)}}END{{print tot}}'
		NO_READS="$(grep '>' {input} | wc -l)"
		echo "Total number of reads: ${{NO_READS}}"
		AVG="$(( ${{BP_READS}} / ${{NO_READS}} ))"
		echo "AVG length of reads: ${{AVG}}"
		NUM="$(( ${{LENGTH}} * {wildcards.depth} / ${{AVG}} ))"
		echo "Number of reads required to achieve depth of {wildcards.depth}: ${{NUM}}"
		MAX_DEPTH="$(( ${{BP_READS}} / ${{LENGTH}} ))"
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
		"1_genome_reads/{prefix}_{kmer}_{cov}.fasta",
	output:
		"2_genome_subset/random/{prefix}_{kmer}_{cov}_{depth}.fasta",
	log:
		"logs/seqtk/genome_{prefix}_random_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/seqtk/assemblytype_genome_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		1
	params:
		length = LENGTH_GENOME,
	conda:
		"envs/seqtk.yaml",
	shell:
		"""
		BP_READS="$(grep -v "^>" {input} | wc | awk "{{print \$3-\$1}}")"
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

#Generate Flye assembly

rule flye_genome:
	input:
		"2_{ass_type}_subset/{read_select}/{prefix}_{kmer}_{cov}_{depth}.fasta",
	output:
		"3_{ass_type}_assembly/flye/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta"
	log:
		"logs/flye/{ass_type}/{read_select}_{prefix}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/flyegenome/assemblytype_{ass_type}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		MAX_THREADS
	resources:
		time = lambda wildcards, input: (2000 if wildcards.ass_type == "genome" else 10),
		mem_mb = lambda wildcards, input: (370000 if wildcards.ass_type == "genome" else 5000),
		cpu = lambda wildcards, input: (MAX_THREADS if wildcards.ass_type == "genome" else 5),
	params:
		out = "3_{ass_type}_assembly/flye/{prefix}/{read_select}_{kmer}_{cov}_{depth}",
		size = LENGTH_GENOME,
	conda:
		"envs/flye.yaml",
	wildcard_constraints:
		ass_type = "genome",
	shell:
		"""
		SIZE={params.size}
		echo "Starting genome assembly ..."
		FIRST="${{SIZE:0:1}}"
		SECOND="${{SIZE:1:1}}"
#		ASM="--asm-coverage 30"
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
		echo "about to run flye"
		(flye --pacbio-raw {input} --genome-size ${{NUSIZE}} --out-dir {params.out} --threads {threads} --asm-coverage 40) 2> {log}
		"""

rule flye:
	input:
		"2_{ass_type}_subset/{read_select}/{prefix}_{kmer}_{cov}_{depth}.fasta",
	output:
		"3_{ass_type}_assembly/flye/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta"
	log:
		"logs/flye/{ass_type}/{read_select}_{prefix}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/flye/assemblytype_{ass_type}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		20
	params:
		out = "3_{ass_type}_assembly/flye/{prefix}/{read_select}_{kmer}_{cov}_{depth}",
		chlor = LENGTH_CHLOROPLAST,
		mito = LENGTH_MITOCHONDRIA,
	wildcard_constraints:
		ass_type = "(mitochondria|chloroplast)",
	resources:
		time = lambda wildcards, input: (60 if wildcards.ass_type == 'mitochondria' else 10),
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
		"benchmarks/raven/assemblytype_{ass_type}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		MAX_THREADS
	resources: 
		time = lambda wildcards, input: (2000 if wildcards.ass_type == "genome" else 10),
		mem_mb = lambda wildcards, input: (350000 if wildcards.ass_type == "genome" else 5000),
		cpu = lambda wildcards, input: (MAX_THREADS if wildcards.ass_type == "genome" else 5),
	conda:
		"envs/raven.yaml",
	shell:
		"""
		raven -t {threads} --polishing 1 --graphical-fragment-assembly {output.gfa} {input} > {output.fasta}
		"""

rule wtdbg2:
	input:
		"2_{ass_type}_subset/{read_select}/{prefix}_{kmer}_{cov}_{depth}.fasta",
	output:
		"3_{ass_type}_assembly/wtdbg2/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
	log:
		"logs/wtdbg2/{ass_type}/{read_select}_{prefix}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/wtdbg2/assemblytype_{ass_type}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		MAX_THREADS
	resources:
		time = lambda wildcards, input: (2000 if wildcards.ass_type == "genome" else 10),
		mem_mb = lambda wildcards, input: (350000 if wildcards.ass_type == "genome" else 5000),
		cpu = lambda wildcards, input: (MAX_THREADS if wildcards.ass_type == "genome" else 5),	
	shadow:	
		"shallow"
	params:
		wtdbg2_1 = WTDBG2_PATH1,
		wtdbg2_2 = WTDBG2_PATH2,
		mito = LENGTH_MITOCHONDRIA,
		chlor = LENGTH_CHLOROPLAST,
		genome = LENGTH_GENOME,
		prefix = "assembly",
		location = "3_{ass_type}_assembly/wtdbg2/{prefix}/{read_select}_{kmer}_{cov}_{depth}/",
	shell:
		"""

		if [ {wildcards.ass_type} == 'chloroplast' ]; then
			SIZE={params.chlor}
		fi		

		if [ {wildcards.ass_type} == 'mitochondria' ]; then
			SIZE={params.mito}
		fi		

		if [ {wildcards.ass_type} == 'genome' ]; then
			SIZE={params.genome}
		fi	
	
		{params.wtdbg2_1} -t {threads} -x sq -g ${{SIZE}} -L 5000 -i {input} -fo {params.prefix}
		{params.wtdbg2_2} -t {threads} -i {params.prefix}.ctg.lay.gz -fo {params.prefix}.ctg.fa
		mv {params.prefix}* {params.location}
		mv {params.location}{params.prefix}.ctg.fa {params.location}{params.prefix}.fasta
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
		"benchmarks/smartdenovo/assemblytype_{ass_type}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		MAX_THREADS
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

#Not sure if job names are necessary to run many jobs in parallel. As it is below, can't run for more than 1 different prefix at a time if job names are necessary. 
rule canu:
	input:
		"2_{ass_type}_subset/{read_select}/{prefix}_{kmer}_{cov}_{depth}.fasta",
	output:
		"3_{ass_type}_assembly/canu/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
	log:
		"logs/canu/{ass_type}/{read_select}_{prefix}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/canu/assemblytype_{ass_type}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}_prefix_{prefix}.tsv"
	threads:
		1
	params:
		prefix = "assembly",
		dir = "3_{ass_type}_assembly/canu/{prefix}/{read_select}_{kmer}_{cov}_{depth}",
		chlor = LENGTH_CHLOROPLAST,
		mito = LENGTH_MITOCHONDRIA,
		genome = LENGTH_GENOME,
		jobName = "clr{depth}{read_select}",
	conda:
		"envs/canu.yaml",
	shell:
		"""
		rm -f {params.dir}/failure
		cp canuFailure.sh {params.dir}
		CORMHAP=""
		REPEAT=""
		ETIME="--time=00:10:00"
		NUMREADCOR=""
		if [ {wildcards.ass_type} == 'chloroplast' ]; then
			SIZE={params.chlor}
			JTIME="--time=00:40:00"
		fi

		if [ {wildcards.ass_type} == 'mitochondria' ]; then
			SIZE={params.mito}
			JTIME="--time=00:40:00"
		fi
		
		if [ {wildcards.ass_type} == 'genome' ]; then
			SIZE={params.genome}
			JTIME="--time=10:00:00"
			ETIME="--time=10:00:00"
			echo "Starting job {params.jobName}, Canu genome assembly."
			(canu -p {params.prefix} -d {params.dir} genomeSize=${{SIZE}} gridOptions=${{JTIME}} gridOptionsJobName={params.jobName} gridOptionsExecutive=${{ETIME}} executiveMemory=4 corMhapFilterThreshold=0.0000000002 gridOptionsCORMHAP="--time=72:00:00" corMhapOptions="--threshold 0.80 --num-hashes 512 --num-min-matches 3 --ordered-sketch-size 1000 --ordered-kmer-size 14 --min-olap-length 2000 --repeat-idf-scale 50" mhapMemory=60g mhapBlockSize=500 ovlMerThreshold=500 -pacbio-raw {input} stopOnLowCoverage=3 onSuccess=touch onFailure=./canuFailure.sh) 2> {log}
		else
			echo "Starting job {params.jobName}, {wildcards.ass_type} assembly."
			(canu -p {params.prefix} -d {params.dir} genomeSize=${{SIZE}} gridOptions=${{JTIME}} gridOptionsJobName={params.jobName} gridOptionsExecutive=${{ETIME}} executiveMemory=4 -pacbio-raw {input} stopOnLowCoverage=3 onSuccess=touch onFailure=./canuFailure.sh) 2> {log}
		fi
		
		while :
		do 
			if [ -f {params.dir}/{params.prefix} ]; then
				echo "Canu job {params.jobName} is finished!"
				mv {params.dir}/assembly.contigs.fasta {params.dir}/assembly.fasta
				break
			fi
			if [ -f {params.dir}/failure ]; then
				echo "Youre not a failure but Canu job {params.jobName} is"
				break
			fi
			sleep 5m 
		done
		
		"""

rule benchcanu:
	input:
		"3_{ass_type}_assembly/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta"
	output:
		"benchcanu/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{prefix}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}/job_finish_state.txt",
		"benchcanu/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{prefix}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}/submit_time.txt",
		"benchcanu/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{prefix}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}/start_time.txt",
		"benchcanu/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{prefix}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}/end_time.txt",
		"benchcanu/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{prefix}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}/walltime_reserved.txt",
		"benchcanu/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{prefix}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}/walltime_elapsed.txt",
		"benchcanu/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{prefix}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}/max_memory.txt",
		"benchcanu/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{prefix}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}/max_disk_write.txt",
		"benchcanu/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{prefix}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}/max_disk_read.txt",
		"benchcanu/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{prefix}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}/num_cores.txt",
	log:
		"logs/benchcanu/{ass_type}/{tool}/{prefix}_{read_select}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/benchcanu/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{prefix}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		1
	params:
		dir = "3_{ass_type}_assembly/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/",
		dir2 = "benchcanu/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{prefix}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}/",
	shell:
		"""
		grep -r 'State               :' {params.dir} > {params.dir2}job_finish_state.txt
		grep -r 'Submit              : 2020' {params.dir} > {params.dir2}submit_time.txt
		grep -r 'Start               : 2020' {params.dir} > {params.dir2}start_time.txt
		grep -r 'End                 : 2020' {params.dir} > {params.dir2}end_time.txt
		grep -r 'Walltime reserved   :' {params.dir} > {params.dir2}walltime_reserved.txt
		grep -r 'Walltime elapsed (%):' {params.dir} > {params.dir2}walltime_elapsed.txt
		grep -r '% Mem used (Max)    :' {params.dir} > {params.dir2}max_memory.txt
		grep -r 'Max Disk Write      :' {params.dir} > {params.dir2}max_disk_write.txt
		grep -r 'Max Disk Read       :' {params.dir} > {params.dir2}max_disk_read.txt
		grep -r 'Cores               :' {params.dir} > {params.dir2}num_cores.txt
		"""
rule contig_length:
	input:
		"3_genome_assembly/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
	output:
		"length/{ass_type}_{tool}_{prefix}_{read_select}_{kmer}_{cov}_{depth}.txt",	
	log:
		"logs/length/{ass_type}_{prefix}_{tool}_{read_select}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"logs/length/assemblytype_{ass_type}_prefix_{prefix}_assemblytool_{tool}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	wildcard_constraints:
		prefix = "m64015_90510_20042",
	shell:
		"""
		cat {input} | awk '$0 ~ ">" {{if (NR > 1) {{print c;}} c=0;printf substr($0,2,100) "\t"; }} $0 !~ ">" {{c+=length($0);}} END {{ print c; }}' > {output}
		"""

rule busco:
	input:
		fasta = "3_genome_assembly/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
	output:
		"busco/assemblytype_genome_assemblytool_{tool}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}/results/run_poales_odb10/full_table.tsv",
	log:
		"logs/busco/genome/{tool}_{read_select}_{prefix}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/busco/assemblytype_genome_readselect_{read_select}_assemblytool_{tool}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		MAX_THREADS
	shadow:
		"full",
	singularity:
		"docker://ezlabgva/busco:v4.0.5_cv1",
	params:
		outdir = "assemblytype_genome_assemblytool_{tool}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}",
		lineage_path = "./reference/poales_odb10",
	shell:
		"""
		(busco -f --in {input.fasta} --out results --lineage_dataset {params.lineage_path} --cpu {threads} --mode genome) 2> {log}
		mv results/short_summary.specific.poales_odb10.results.txt busco/{params.outdir}/results
		mv results/blast_db busco/{params.outdir}/results
#		mv results/logs busco/{params.outdir}/results
		mv results/run_poales_odb10/* busco/{params.outdir}/results/run_poales_odb10/
		"""

rule busco_polished:
	input:
		fasta = "4_genome_polished/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/{polishdepth}_{polishselect}/{polishround}_assembly.fasta",
	output:
		"busco/assemblytype_genome_assemblytool_{tool}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishselect_{polishselect}_polishdepth_{polishdepth}_polishround_{polishround}/results/run_poales_odb10/full_table.tsv",
	log:
		"logs/busco_polished/genome/{tool}_{read_select}_{prefix}_{kmer}_{cov}_{depth}_{polishselect}_{polishdepth}_{polishround}.log",
	benchmark:
		"benchmarks/buscopolished/assemblytype_genome_readselect_{read_select}_assemblytool_{tool}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishselect_{polishselect}_polishdepth_{polishdepth}_polishround_{polishround}.tsv",
	threads:
		MAX_THREADS
	shadow:
		"full",
	singularity:
		"docker://ezlabgva/busco:v4.0.5_cv1",
	params:
		outdir = "assemblytype_genome_assemblytool_{tool}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishselect_{polishselect}_polishdepth_{polishdepth}_polishround_{polishround}",
		lineage_path = "./reference/poales_odb10",
	shell:
		"""
		(busco -f --in {input.fasta} --out results --lineage_dataset {params.lineage_path} --cpu {threads} --mode genome) 2> {log}
		mv results/short_summary.specific.poales_odb10.results.txt busco/{params.outdir}/results
		mv results/blast_db busco/{params.outdir}/results
#		mv results/logs busco/{params.outdir}/results
		mv results/run_poales_odb10/* busco/{params.outdir}/results/run_poales_odb10/
		"""

rule gt_ltrharvest:
	input:
		"3_{ass_type}_assembly/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
	output:
		"ltr/harvest/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}/assembly.fa.harvest.scn",
	log:
		"logs/gt_ltrharvest/{ass_type}/{tool}_{read_select}_{prefix}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/gtltrharvest/assemblytype_{ass_type}_readselect_{read_select}_assemblytool_{tool}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		1
	shadow:
		"shallow",
	wildcard_constraints:
		ass_type = "genome",
	conda:
		"envs/genometools.yaml",
	params:
		dir = "ltr/harvest/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}/",
	shell:
		"""
			gt suffixerator -db {input} -indexname suffixer.fa -tis -suf -lcp -des -ssp -sds -dna
			gt ltrharvest -index suffixer.fa -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes > out.scn
			mv  out.scn {output}
			mv suffixer* {params.dir}
		"""

rule gt_ltrharvest_polished:
	input:
		"4_{ass_type}_polished/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/{polishdepth}_{polishselect}/{polishround}_assembly.fasta",
	output:
		"ltr/harvest/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishselect_{polishselect}_polishdepth_{polishdepth}_polishround_{polishround}/assembly.fa.harvest.scn",
	log:
		"logs/gt_ltrharvest_polished/{ass_type}/{tool}_{read_select}_{prefix}_{kmer}_{cov}_{depth}_{polishselect}_{polishdepth}_{polishround}.log",
	benchmark:
		"benchmarks/gtltrharvestpolished/assemblytype_{ass_type}_readselect_{read_select}_assemblytool_{tool}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishselect_{polishselect}_polishdepth_{polishdepth}_polishround_{polishround}.tsv",
	threads:
		1
	shadow:
		"shallow",
	wildcard_constraints:
		ass_type = "genome",
	conda:
		"envs/genometools.yaml",
	params:
		dir = "ltr/harvest/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishselect_{polishselect}_polishdepth_{polishdepth}_polishround_{polishround}/",
	shell:
		"""
			gt suffixerator -db {input} -indexname suffixer.fa -tis -suf -lcp -des -ssp -sds -dna
			gt ltrharvest -index suffixer.fa -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes > out.scn
			mv  out.scn {output}
			mv suffixer* {params.dir}
		"""

rule ltr_finder:
	input:
		"3_{ass_type}_assembly/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
	output:
		"ltr/finder/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}/assembly.fasta.finder.combine.scn",
	log:
		"logs/ltr_finder/{ass_type}/{tool}_{read_select}_{prefix}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/ltrfinder/assemblytype_{ass_type}_readselect_{read_select}_assemblytool_{tool}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		10
	shadow:
		"shallow",
	wildcard_constraints:
		ass_type = "genome",
	params:
		dir = "ltr/finder/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}/",
		finder = LTR_FINDER_PATH,
	shell:
		"""
		perl {params.finder} -seq {input} -threads {threads} -harvest_out -size 1000000 -time 300
		mv assembly* {params.dir}
		"""

rule ltr_finder_polished:
	input:
		"4_{ass_type}_polished/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/{polishdepth}_{polishselect}/{polishround}_assembly.fasta",
	output:
		"ltr/finder/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishselect_{polishselect}_polishdepth_{polishdepth}_polishround_{polishround}/{polishround}_assembly.fasta.finder.combine.scn",
	log:
		"logs/ltr_finder_polished/{ass_type}/{tool}_{read_select}_{prefix}_{kmer}_{cov}_{depth}_{polishselect}_{polishdepth}_{polishround}.log",
	benchmark:
		"benchmarks/ltrfinderpolished/assemblytype_{ass_type}_readselect_{read_select}_assemblytool_{tool}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishselect_{polishselect}_polishdepth_{polishdepth}_polishround_{polishround}.tsv",
	threads:
		10
	shadow:
		"shallow",
	wildcard_constraints:
		ass_type = "genome",
	params:
		dir = "ltr/finder/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishselect_{polishselect}_polishdepth_{polishdepth}_polishround_{polishround}/",
		finder = LTR_FINDER_PATH,
	shell:
		"""
		perl {params.finder} -seq {input} -threads {threads} -harvest_out -size 1000000 -time 300
		mv {wildcards.polishround}_assembly* {params.dir}
		"""

rule ltr_retriever:
	input:
		genome = "3_{ass_type}_assembly/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
		finder = "ltr/finder/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}/assembly.fasta.finder.combine.scn",
		harvest = "ltr/harvest/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}/assembly.fa.harvest.scn",
	output:
		scn = "ltr/retriever/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}/rawLTR.scn",
		lai =  "ltr/retriever/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}/LAI_score.txt",
		dummy = "ltr/retriever/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}/complete",
	log:
		"logs/ltr_retriever/{ass_type}/{tool}_{read_select}_{prefix}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/ltrretriever/assemblytype_{ass_type}_readselect_{read_select}_assemblytool_{tool}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		10
	shadow:
		"shallow",
	wildcard_constraints:
		ass_type = "genome",
	conda:
		"envs/ltr.yaml",
	params:
		dir = "ltr/retriever/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}/",
	shell:
		"""
		cat {input.finder} {input.harvest} > {output.scn}
		rm -f {params.dir}assembly.fasta
		cp {input.genome} {params.dir}
		LTR_retriever -genome {params.dir}assembly.fasta -inharvest {output.scn} -threads {threads}
		
		mv assembly.fasta.*.LAI {params.dir}LAI_score.txt
		touch {output.dummy}
		"""


rule ltr_retriever_polished:
	input:
		genome = "4_{ass_type}_polished/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/{polishdepth}_{polishselect}/{polishround}_assembly.fasta",
		finder = "ltr/finder/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishselect_{polishselect}_polishdepth_{polishdepth}_polishround_{polishround}/{polishround}_assembly.fasta.finder.combine.scn",
		harvest = "ltr/harvest/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishselect_{polishselect}_polishdepth_{polishdepth}_polishround_{polishround}/assembly.fa.harvest.scn",
	output:
		scn = "ltr/retriever/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishselect_{polishselect}_polishdepth_{polishdepth}_polishround_{polishround}/rawLTR.scn",
		lai =  "ltr/retriever/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishselect_{polishselect}_polishdepth_{polishdepth}_polishround_{polishround}/LAI_score.txt",
		dummy = "ltr/retriever/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishselect_{polishselect}_polishdepth_{polishdepth}_polishround_{polishround}/complete",
	log:
		"logs/ltr_retriever_polished/{ass_type}/{tool}_{read_select}_{prefix}_{kmer}_{cov}_{depth}_{polishselect}_{polishdepth}_{polishround}.log",
	benchmark:
		"benchmarks/ltrretrieverpolished/assemblytype_{ass_type}_readselect_{read_select}_assemblytool_{tool}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishselect_{polishselect}_polishdepth_{polishdepth}_polishround_{polishround}.tsv",
	threads:
		10
	shadow:
		"shallow",
	wildcard_constraints:
		ass_type = "genome",
		polishround = "(1|2)",
	conda:
		"envs/ltr.yaml",
	params:
		dir = "ltr/retriever/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishselect_{polishselect}_polishdepth_{polishdepth}_polishround_{polishround}/",
	shell:
		"""
		cat {input.finder} {input.harvest} > {output.scn}
		rm -f {params.dir}assembly.fasta
		cp {input.genome} {params.dir}
		LTR_retriever -genome {params.dir}{wildcards.polishround}_assembly.fasta -inharvest {output.scn} -threads {threads}
		
		mv {wildcards.polishround}_assembly.fasta.*.LAI {params.dir}LAI_score.txt
		touch {output.dummy}
		"""

rule quast:
	input:
		"3_{ass_type}_assembly/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta"
	output:
		"quast/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{prefix}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}/report.tsv",
	log:
		"logs/quast/{ass_type}/{tool}/{prefix}_{read_select}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/quast/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{prefix}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		MAX_THREADS
	params:
		out = "quast/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{prefix}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}",
#	conda:
#		"envs/quast.yaml",
	singularity:
		"docker://quay.io/biocontainers/quast:5.0.2--1"
	resources:
		time = lambda wildcards, input: (15 if wildcards.ass_type == "genome" else 1),
		mem_mb = lambda wildcards, input: (3000 if wildcards.ass_type == "genome" else 3000),
		cpu = lambda wildcards, input: (2 if wildcards.ass_type == "genome" else 1),
	shell:
		"""
		quast {input} --threads {threads} -o {params.out}
		"""


rule quastref:
	input:
		ass = "3_{ass_type}_assembly/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
		ref = "reference/{ass_type}.fasta",
	output:
		"quastref/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{prefix}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}/report.tsv",
	log:
		"logs/quastref/{ass_type}/{tool}/{prefix}_{read_select}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/quastref/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{prefix}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		MAX_THREADS
	params:
		out = "quastref/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{prefix}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}",
	singularity:
		"docker://quay.io/biocontainers/quast:5.0.2--1"
	resources:
		time = lambda wildcards, input: (60 if wildcards.ass_type == "genome" else 1),
		mem_mb = lambda wildcards, input: (15000 if wildcards.ass_type == "genome" else 3000),
		cpu = lambda wildcards, input: (5 if wildcards.ass_type == "genome" else 1),
	shell:
		"""
		if [ {wildcards.ass_type} == 'genome' ]; then
			quast {input.ass} --eukaryote --no-icarus --no-html --large -r {input.ref} --threads {threads} -o {params.out}
		else
			quast {input.ass} --no-icarus --no-html -r {input.ref} --threads {threads} -o {params.out}
		fi
		"""


rule quastref_polished:
	input:
		ass = "4_{ass_type}_polished/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/{polishdepth}_{polishselect}/{polishround}_assembly.fasta",
		ref = "reference/{ass_type}.fasta",
	output:
		"quastref/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{prefix}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishselect_{polishselect}_polishdepth_{polishdepth}_polishround_{polishround}/report.tsv",
	log:
		"logs/quastref_polished/{ass_type}/{tool}/{prefix}_{read_select}_{kmer}_{cov}_{depth}_{polishselect}_{polishdepth}_{polishround}.log",
	benchmark:
		"benchmarks/quastrefpolished/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{prefix}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishselect_{polishselect}_polishdepth_{polishdepth}_polishround_{polishround}.tsv",
	threads:
		MAX_THREADS
	params:
		out = "quastref/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{prefix}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishselect_{polishselect}_polishdepth_{polishdepth}_polishround_{polishround}",
	singularity:
		"docker://quay.io/biocontainers/quast:5.0.2--1"
	resources:
		time = lambda wildcards, input: (60 if wildcards.ass_type == "genome" else 1),
		mem_mb = lambda wildcards, input: (15000 if wildcards.ass_type == "genome" else 3000),
		cpu = lambda wildcards, input: (5 if wildcards.ass_type == "genome" else 1),
	shell:
		"""
		if [ {wildcards.ass_type} == 'genome' ]; then
			quast {input.ass} --eukaryote --no-icarus --no-html --large -r {input.ref} --threads {threads} -o {params.out}
		else
			quast {input.ass} --no-icarus --no-html -r {input.ref} --threads {threads} -o {params.out}
		fi
		"""

rule mummer:
	input:
		chlor = "reference/chloroplast.fasta",
		mito = "reference/mitochondria.fasta",
		genome = "reference/genome.fasta",
		query = "3_{ass_type}_assembly/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
	output:
		mums = "mummer/{ass_type}/prefix_{prefix}_assemblytool_{tool}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}.delta",
		gp = "mummer/{ass_type}/prefix_{prefix}_assemblytool_{tool}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}.gp",
	log:
		"logs/mummer/{prefix}_{ass_type}/{read_select}_assemblytool_{tool}_ref_query_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/mummer/prefix_{prefix}_assemblytool_{tool}_assemblytype_{ass_type}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		MAX_THREADS
	shadow:
		"shallow"
	resources:
		time = lambda wildcards, input: (100 if wildcards.ass_type == "genome" else 1),
		mem_mb = lambda wildcards, input: (30000 if wildcards.ass_type == "genome" else 200),
		cpu = lambda wildcards, input: (10 if wildcards.ass_type == "genome" else 1),
	params:
		pref = "prefix_{prefix}_assemblytool_{tool}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}",
		mummer = MUMMER_PATH,
	shell:
		"""
		if [ {wildcards.ass_type} == "chloroplast" ]; then
			REF={input.chlor}
		fi

		if [ {wildcards.ass_type} == "mitochondria" ]; then
			REF={input.mito}
		fi

		if [ {wildcards.ass_type} == "genome" ]; then
			REF={input.genome}
		fi

		({params.mummer}nucmer -t {threads} --prefix={params.pref} ${{REF}} {input.query}) 2> {log}
		{params.mummer}show-coords -r -c -H -d -o -T -l {params.pref}.delta > {params.pref}.coords
		{params.mummer}show-snps -C -l -r -T -H {params.pref}.delta > {params.pref}.snps
		{params.mummer}show-tiling {params.pref}.delta > {params.pref}.tiling
		{params.mummer}mummerplot --postscript --prefix={params.pref} {params.pref}.delta
		mv {params.pref}.* mummer/{wildcards.ass_type}

		"""


rule mummer_polished:
	input:
		chlor = "reference/chloroplast.fasta",
		mito = "reference/mitochondria.fasta",
		genome = "reference/genome.fasta",
		query = "4_{ass_type}_polished/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/{polishdepth}_{polishselect}/{polishround}_assembly.fasta",
	output:
		mums = "mummer/{ass_type}/prefix_{prefix}_assemblytool_{tool}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishselect_{polishselect}_polishdepth_{polishdepth}_polishround_{polishround}.delta",
		gp = "mummer/{ass_type}/prefix_{prefix}_assemblytool_{tool}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishselect_{polishselect}_polishdepth_{polishdepth}_polishround_{polishround}.gp",
	log:
		"logs/mummer_polished/{prefix}_{ass_type}/{read_select}_assemblytool_{tool}_ref_query_{kmer}_{cov}_{depth}_polishselect_{polishselect}_{polishdepth}_{polishround}.log",
	benchmark:
		"benchmarks/mummerpolished/prefix_{prefix}_assemblytool_{tool}_assemblytype_{ass_type}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishselect_{polishselect}_polishdepth_{polishdepth}_polishround_{polishround}.tsv",
	threads:
		MAX_THREADS
	shadow:
		"shallow"
	resources:
		time = lambda wildcards, input: (100 if wildcards.ass_type == "genome" else 1),
		mem_mb = lambda wildcards, input: (30000 if wildcards.ass_type == "genome" else 200),
		cpu = lambda wildcards, input: (10 if wildcards.ass_type == "genome" else 1),
	params:
		pref = "prefix_{prefix}_assemblytool_{tool}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishselect_{polishselect}_polishdepth_{polishdepth}_polishround_{polishround}",
		mummer = MUMMER_PATH,
	shell:
		"""
		if [ {wildcards.ass_type} == "chloroplast" ]; then
			REF={input.chlor}
		fi

		if [ {wildcards.ass_type} == "mitochondria" ]; then
			REF={input.mito}
		fi

		if [ {wildcards.ass_type} == "genome" ]; then
			REF={input.genome}
		fi

		({params.mummer}nucmer -t {threads} --prefix={params.pref} ${{REF}} {input.query}) 2> {log}
		{params.mummer}show-coords -r -c -H -d -o -T -l {params.pref}.delta > {params.pref}.coords
		{params.mummer}show-snps -C -l -r -T -H {params.pref}.delta > {params.pref}.snps
		{params.mummer}show-tiling {params.pref}.delta > {params.pref}.tiling
		{params.mummer}mummerplot --postscript --prefix={params.pref} {params.pref}.delta
		mv {params.pref}.* mummer/{wildcards.ass_type}

		"""

rule minimap2_1:
	input:
		draft = "3_{ass_type}_assembly/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
		reads = "2_{ass_type}_subset/{polishselect}/{prefix}_{kmer}_{cov}_{polishdepth}.fasta",
	output:
		temp("3_{ass_type}_assembly/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/clr_mapped_{polishselect}_{polishdepth}.sam"),
	log:
		"logs/minimap2_1/{ass_type}_{tool}_{prefix}_{read_select}_{kmer}_{cov}_{depth}_{polishselect}_{polishdepth}.log",
	benchmark:
		"benchmarks/minimap2run1/assemblytype_{ass_type}_readselect_{read_select}_assemblytool_{tool}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishdepth_{polishdepth}_polishselect_{polishselect}.tsv",
	conda:
		"envs/minimap2.yaml",
	params:
		dir = "3_{ass_type}_assembly/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}",
	resources:
		time = lambda wildcards, input: (600 if wildcards.ass_type == "genome" else 1),
		mem_mb = lambda wildcards, input: (30000 if wildcards.ass_type == "genome" else 3000),
		cpu = lambda wildcards, input: (10 if wildcards.ass_type == "genome" else 1),
	threads:
		MAX_THREADS
	shell:
		"""
		minimap2 -ax map-pb --eqx -m 5000 -t {threads} --secondary=no {input.draft} {input.reads} | samtools view -F 1796 - > {output}
		"""

rule racon_1:
	input:
		sam = "3_{ass_type}_assembly/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/clr_mapped_{polishselect}_{polishdepth}.sam",
		draft = "3_{ass_type}_assembly/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
		reads = "2_{ass_type}_subset/{polishselect}/{prefix}_{kmer}_{cov}_{polishdepth}.fasta",
	output:
		"4_{ass_type}_polished/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/{polishdepth}_{polishselect}/1_assembly.fasta",
	log:
		"logs/racon_1/{ass_type}_{tool}_{prefix}_{read_select}_{kmer}_{cov}_{depth}_{polishdepth}_{polishselect}.log",
	benchmark:
		"benchmarks/racon1/assemblytype_{ass_type}_readselect_{read_select}_assemblytool_{tool}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishdepth_{polishdepth}_polishselect_{polishselect}.tsv",
	resources:
		time = lambda wildcards, input: (400 if wildcards.ass_type == "genome" else 1),
		mem_mb = lambda wildcards, input: (60000 if wildcards.ass_type == "genome" else 3000),
		cpu = lambda wildcards, input: (20 if wildcards.ass_type == "genome" else 1),
	conda:
		"envs/racon.yaml",
	threads:
		MAX_THREADS
	shell:
		"""
		racon {input.reads} {input.sam} {input.draft} -u -t {threads} > {output}
		"""

rule minimap2_2:
	input:
		draft = "4_{ass_type}_polished/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/{polishdepth}_{polishselect}/1_assembly.fasta",
		reads = "2_{ass_type}_subset/{polishselect}/{prefix}_{kmer}_{cov}_{polishdepth}.fasta",
	output:
		temp("4_{ass_type}_polished/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/{polishdepth}_{polishselect}/clr_mapped.sam"),
	log:
		"logs/minimap2_2/{ass_type}_{tool}_{prefix}_{read_select}_{kmer}_{cov}_{depth}_{polishselect}_{polishdepth}.log",
	benchmark:
		"benchmarks/minimap2run2/assemblytype_{ass_type}_readselect_{read_select}_assemblytool_{tool}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishdepth_{polishdepth}_polishselect_{polishselect}.tsv",
	conda:
		"envs/minimap2.yaml",
	params:
		dir = "4_{ass_type}_polished/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/{polishdepth}_{polishselect}",
	resources:
		time = lambda wildcards, input: (600 if wildcards.ass_type == "genome" else 1),
		mem_mb = lambda wildcards, input: (30000 if wildcards.ass_type == "genome" else 3000),
		cpu = lambda wildcards, input: (10 if wildcards.ass_type == "genome" else 1),
	threads:
		MAX_THREADS
	shell:
		"""
		minimap2 -ax map-pb --eqx -m 5000 -t {threads} --secondary=no {input.draft} {input.reads} | samtools view -F 1796 - > {output}
		"""

rule racon_2:
	input:
		sam = "4_{ass_type}_polished/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/{polishdepth}_{polishselect}/clr_mapped.sam",
		draft = "4_{ass_type}_polished/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/{polishdepth}_{polishselect}/1_assembly.fasta",
		reads = "2_{ass_type}_subset/{polishselect}/{prefix}_{kmer}_{cov}_{polishdepth}.fasta",
	output:
		"4_{ass_type}_polished/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/{polishdepth}_{polishselect}/2_assembly.fasta",
	log:
		"logs/racon_2/{ass_type}_{tool}_{prefix}_{read_select}_{kmer}_{cov}_{depth}_{polishdepth}_{polishselect}.log",
	benchmark:
		"benchmarks/racon2/assemblytype_{ass_type}_readselect_{read_select}_assemblytool_{tool}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}_polishdepth_{polishdepth}_polishselect_{polishselect}.tsv",
	resources:
		time = lambda wildcards, input: (400 if wildcards.ass_type == "genome" else 1),
		mem_mb = lambda wildcards, input: (60000 if wildcards.ass_type == "genome" else 3000),
		cpu = lambda wildcards, input: (20 if wildcards.ass_type == "genome" else 1),
	conda:
		"envs/racon.yaml",
	threads:
		MAX_THREADS
	shell:
		"""
		racon {input.reads} {input.sam} {input.draft} -u -t {threads} > {output}
		"""

#=====================================================
# Chloroplast assembly
#=====================================================

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
