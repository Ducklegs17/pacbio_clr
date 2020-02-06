

localrules:
	all

rule all:
	input:
		"SequelToolsResults/summaryTable.txt",
		"m64015_90510_20042.subreads.sam",

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
		"envs/sequelTools.yaml",
	shell:
		"""
		(bash SequelTools.sh -t Q -v -u {input}) 2> {log} 
		"""

rule sequelToolsSubset:
	input:
		bam = "0_raw/m64015_90510_20042.subreads.bam.1",
		file = "0_raw/bamList.txt",
	output:
		"m64015_90510_20042.subreads.sam",
	log:
		"logs/sequelToolsSubset/sequelToolsSubset.log",
	benchmark:
		"benchmarks/sequelToolsSubset.tsv",
	threads:
		1
	conda:
		"envs/sequelTools.yaml",
	shell:
		"""
		(bash SequelTools.sh -t S -v -u {input.file} -T r -f s -o 0_raw/subset) 2> {log}
		"""
	
