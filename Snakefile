# Instructions:

PREFIXES = ["m64015_90510_20042"]
MITOCONDRIA_REF = "reference/mitochondria.fasta"
CHLOROPLASTS_REF = "reference/chloroplast.fasta",
READ_COVERAGE = ["25"]
MAX_THREADS = 32
BBDUK_KMER_LENGTH = ["17"]
BBDUK_MIN_COVERAGE = ["0.7"]
LENGTH_CHLOROPLAST = ["134502"]
LENGTH_MITOCHONDRIA = ["415805"]
LENGTH_GENOME = ["5500000"]
ASSEMBLY_TOOLS = ["flye","raven","smartdenovo","canu"]
ASSEMBLY_TYPE = ["genome"]
READ_SELECTION_METHOD = ["longest"]

localrules: 
	all,
	canu,
	generate_coverage_list,

rule all:
	input:
		"SequelToolsResults/summaryTable.txt",
		expand("1_{ASS_TYPE}_reads/{PREFIX}_{KMER}_{COV}.fasta", ASS_TYPE = ASSEMBLY_TYPE, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE),
		expand("1_{ASS_TYPE}_reads/sorted/{PREFIX}_{KMER}_{COV}.fasta", ASS_TYPE = ASSEMBLY_TYPE, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE),
		expand("1_{ASS_TYPE}_reads/sorted/{PREFIX}_{KMER}_{COV}_coveragetable.txt", ASS_TYPE = ASSEMBLY_TYPE, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE),
		expand("2_{ASS_TYPE}_subset/{READ_SELECTION}/{PREFIX}_{KMER}_{COV}_{DEPTH}.fasta", ASS_TYPE = ASSEMBLY_TYPE, READ_SELECTION = READ_SELECTION_METHOD, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE, DEPTH = READ_COVERAGE),
		expand("3_{ASS_TYPE}_assembly/{TOOL}/{PREFIX}/{READ_SELECTION}_{KMER}_{COV}_{DEPTH}/assembly.fasta", ASS_TYPE = ASSEMBLY_TYPE, TOOL = ASSEMBLY_TOOLS, READ_SELECTION = READ_SELECTION_METHOD, PREFIX = PREFIXES, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE, DEPTH = READ_COVERAGE),
		expand("mummer/{ASS_TYPE}/prefix_{PREFIX}_assemblytool_{TOOL}_readselect_{READ_SELECTION}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}.mums", PREFIX = PREFIXES, ASS_TYPE = ASSEMBLY_TYPE, TOOL = ASSEMBLY_TOOLS, READ_SELECTION = READ_SELECTION_METHOD, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE, DEPTH = READ_COVERAGE), 
		expand("quast/assemblytype_{ASS_TYPE}_assemblytool_{TOOL}_prefix_{PREFIX}_readselect_{READ_SELECTION}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}/report.tsv", ASS_TYPE = ASSEMBLY_TYPE, TOOL = ASSEMBLY_TOOLS, PREFIX = PREFIXES, READ_SELECTION = READ_SELECTION_METHOD, KMER = BBDUK_KMER_LENGTH, COV = BBDUK_MIN_COVERAGE, DEPTH = READ_COVERAGE),

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
		"2_genome_subset/{read_select}/{prefix}_{kmer}_{cov}_{depth}.fasta",
	output:
		"3_genome_assembly/flye/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta"
	log:
		"logs/flye/genome/{read_select}_{prefix}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/flyegenome/assemblytype_genome_readselect_{read_select}_prefix_{prefix}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		MAX_THREADS
	params:
		out = "3_genome_assembly/flye/{prefix}/{read_select}_{kmer}_{cov}_{depth}",
		size = LENGTH_GENOME,
	conda:
		"envs/flye.yaml",
	shell:
		"""
		SIZE={params.size}
		echo "Starting genome assembly ..."
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
		jobName = "{depth}{kmer}{cov}",
#	conda:
#		"envs/canu.yaml",
	shell:
		"""
		cp canuFailure.sh {params.dir}

		if [ {wildcards.ass_type} == 'chloroplast' ]; then
			SIZE={params.chlor}
			JTIME="--time=00:40:00"
			echo "Starting job {params.jobName}, chloroplast assembly."
		fi

		if [ {wildcards.ass_type} == 'mitochondria' ]; then
			SIZE={params.mito}
			JTIME="--time=00:40:00"
			echo "Starting job {params.jobName}, mitochondria assembly."
		fi
		
		if [ {wildcards.ass_type} == 'genome' ]; then
			SIZE={params.genome}
			JTIME="--time=04:00:00"
			echo "Starting job {params.jobName}, genome assembly."
		fi
		
		(canu -p {params.prefix} -d {params.dir} genomeSize=${{SIZE}} gridOptions=${{JTIME}} gridOptionsJobName={params.jobName} gridOptionsExecutive="--time=00:10:00" executiveMemory=1 -pacbio-raw {input} onSuccess=touch onFailure=./canuFailure.sh) 2> {log}

		while :
		do 
			if [ -f {params.dir}/{params.prefix} ]; then
				echo "Canu job {params.jobName} is finished!"
				mv {params.dir}/assembly.contigs.fasta {params.dir}/assembly.fasta
				break
			fi
			if [ -f {params.dir}/failure ]; then
				echo "You're not a failure but Canu job {params.jobName} is"
				break
			fi
			sleep 5m 
		done
		
		"""

#rule arrow_polish:
#	input:
#		"3_{ass_type}_assembly/canu/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
#	output:
#		"4_{ass_type}_polished/canu/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
#	log:
#		
#	benchmark:
#
#	params:
#		
#	conda:
#		"envs/pbgcpp.yaml",
#	shell:
#		"""
#		quiver aligned_reads.bam 
#		"""

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
		1
	params:
		out = "quast/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{prefix}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}",
	conda:
		"envs/quast.yaml",
	shell:
		"""
		quast {input} --threads {threads} -o {params.out}
		"""

rule mummer:
	input:
		chlor = "reference/chloroplast.fasta",
		mito = "reference/mitochondria.fasta",
		query = "3_{ass_type}_assembly/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
	output:
		mums = "mummer/{ass_type}/prefix_{prefix}_assemblytool_{tool}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}.mums",
		gp = "mummer/{ass_type}/prefix_{prefix}_assemblytool_{tool}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}.gp",
	log:
		"logs/mummer/{prefix}_{ass_type}/{read_select}_assemblytool_{tool}_ref_query_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/mummer/prefix_{prefix}_assemblytool_{tool}_assemblytype_{ass_type}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		1
	shadow:
		"shallow"
	params:
		pref = "prefix_{prefix}_assemblytool_{tool}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}",
	conda:
		"envs/mummer.yaml",
	shell:
		"""

		if [ {wildcards.ass_type} == "chloroplast" ]; then

			mummer -mum -b -c {input.chlor} {input.query} > {params.pref}.mums
			mummerplot --postscript --prefix={params.pref} {params.pref}.mums
			nucmer -maxmatch -c 100 -p {params.pref} {input.chlor} {input.query}
			show-coords -r -c -H -d -o -T -l {params.pref}.delta > {params.pref}.coords
			show-snps -C -l -r -T -H {params.pref}.delta > {params.pref}.snps
			show-tiling {params.pref}.delta > {params.pref}.tiling
			mv {params.pref}.* mummer/{wildcards.ass_type}
		fi

		if [ {wildcards.ass_type} == "mitochondria" ]; then
			mummer -mum -b -c {input.mito} {input.query} > {params.pref}.mums
			mummerplot --postscript --prefix={params.pref} {params.pref}.mums
			nucmer -maxmatch -c 100 -p {params.pref} {input.chlor} {input.query}
			show-coords -r -c -H -d -o -T -l {params.pref}.delta > {params.pref}.coords
			show-snps -C -l -r -T -H {params.pref}.delta > {params.pref}.snps
			show-tiling {params.pref}.delta > {params.pref}.tiling
			mv {params.pref}.* mummer/{wildcards.ass_type}
		fi
		"""

rule run_mummer3:
	input:
		chlor = "reference/chloroplast.fasta",
		mito = "reference/mitochondria.fasta",
		query = "3_{ass_type}_assembly/{tool}/{prefix}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
	output:
		mums = "runmummer3/{ass_type}/prefix_{prefix}_assemblytool_{tool}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}.mums",
		gp = "runmummer3/{ass_type}/prefix_{prefix}_assemblytool_{tool}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}.gp",
	log:
		"logs/run_mummer3/{prefix}_{ass_type}/{read_select}_assemblytool_{tool}_ref_query_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/runmummer3/prefix_{prefix}_assemblytool_{tool}_assemblytype_{ass_type}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		1
	shadow:
		"shallow"
	params:
		pref = "prefix_{prefix}_assemblytool_{tool}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}",
	conda:
		"envs/mummer.yaml",
	shell:
		"""


		if [ {wildcards.ass_type} == "chloroplast" ]; then

			run-mummer3 {input.chlor} {input.query} {params.pref}.mums
			mummerplot --postscript --prefix={params.pref} {params.pref}.mums
			nucmer -maxmatch -c 100 -p {params.pref} {input.chlor} {input.query}
			show-coords -r -c -H -d -o -T -l {params.pref}.delta > {params.pref}.coords
			show-snps -C -l -r -T -H {params.pref}.delta > {params.pref}.snps
			show-tiling {params.pref}.delta > {params.pref}.tiling
			mv {params.pref}.* mummer/{wildcards.ass_type}
		fi

		if [ {wildcards.ass_type} == "mitochondria" ]; then
			mummer -mum -b -c {input.mito} {input.query} > {params.pref}.mums
			mummerplot --postscript --prefix={params.pref} {params.pref}.mums
			nucmer -maxmatch -c 100 -p {params.pref} {input.chlor} {input.query}
			show-coords -r -c -H -d -o -T -l {params.pref}.delta > {params.pref}.coords
			show-snps -C -l -r -T -H {params.pref}.delta > {params.pref}.snps
			show-tiling {params.pref}.delta > {params.pref}.tiling
			mv {params.pref}.* mummer/{wildcards.ass_type}
		fi

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
