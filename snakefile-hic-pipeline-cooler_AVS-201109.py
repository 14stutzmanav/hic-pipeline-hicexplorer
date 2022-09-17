import pandas as pd
import shutil
import os
import re
import glob
import getpass
import sys
import argparse

configfile: './slurmConfig.json'

# Set Genome and Bin Size (edit this section if you're not working with dm6):
GenomeAssembly = 'dm6'
restrictionDigestPATH = str('/proj/mckaylab/users/astutzman/hic-analyses/dm6_arima_digest.bed')
chrsizePATH = str('/proj/mckaylab/users/astutzman/hic-analyses/juicer/AS1-21_juicer-output/dm6.chrom.sizes')
# Binsize is currently = 1. You can change this in line 121.

########################

# Set Module Versions:

bwaVer = str('bwa/0.7.17')
samtoolsVer = str('samtools/1.9')
pairtoolsVer = str('pairtools/0.3.0')
coolerVer = str('cooler/0.8.6')
hiceVer = str('hicexplorer/3.3.1')

#########################

def getTechReplicates(wildcards):
	"""
	Get the names of fastq files that will be concatenated together as technical replicates.
	Since data is paired-end, keep the read number in the file name.
	If there are no technical replicates, this simply cats the file to a new one w/ a different name
	"""
	readNumRegex = '_R{}'.format(wildcards.Num)

	# wildcards.tech is the prefix of the sample w/o the read number (1 or 2)
	techFilter = sampleDF [ sampleDF.techName == wildcards.tech ]
	fastqList = list(techFilter.fastqName)

	# we want to concatenate R1 technical reps and R2 technical reps seperately
	fastqList = [ fastqName for fastqName in fastqList 
			if re.search(readNumRegex, fastqName) ]

	# the fastq entry in the samplesheet does not contain the directory
	fastqList = [ 'Fastq/{}'.format(fastqName) for fastq in fastqList ]

	return(fastqList)

#########################

# Load Required Input Files:

BWAIndexPATH = str('/proj/seq/data/' + GenomeAssembly + '_UCSC/Sequence/BWAIndex/genome.fa')
sampleSheetPath = str('master-samplesheet-hic.csv')
sampleDF = pd.read_csv(sampleSheetPath, comment = '#')

techList = list(set(sampleDF.techName))
fastqList= list(set(sampleDF.fastqName))
readNumList = list(sampleDF.readNum)


#########################
localrules: all

rule all:
	input:
		expand("Fastq/{fastq}.fastq.gz", fastq=fastqList),
		expand("Fastq/{tech}_{readNum}.fastq.gz", tech=techList, readNum=readNumList),
		expand("Bam/{tech}_cooler.bam", tech=techList),
		expand("Pairs/{tech}_parsed.pairsam", tech=techList),
		expand("Pairs/{tech}_sorted.pairsam", tech=techList),
		expand("Pairs/{tech}_dedup.pairsam", tech=techList),
		expand("Pairs/{tech}.pairs", tech=techList),
		expand("Cool/{tech}.cool", tech=techList),
		expand("Cool/{tech}.mcool", tech=techList),
		expand("Bam/{tech}_mate{readNum}.bam", tech=techList, readNum=readNumList),
		expand("Mat/{tech}_hicMatrix.h5", tech =techList)


rule copyFiles:
	input:
		lambda x: list(sampleDF.htsfFile)
	output:
		expand('Fastq/{fastq}.fastq.gz', fastq = list(sampleDF.fastqName))
	message: "Copying files to Fastq directory with corrected file names."
	run:
		for htsf in list(sampleDF.htsfFile):
			#print('HERE', htsf)
			outFileFilt = sampleDF [ sampleDF.htsfFile == htsf ] 
			outFileBase = list(outFileFilt.fastqName)[0]
			outFile = 'Fastq/{fastq}.fastq.gz'.format(fastq = outFileBase)
			#print('THERE')
			#print(outFile)
			shutil.copyfile(htsf, outFile)
			print('copied file')

#rule combineTechReps:
#	input:
#		getTechReplicates
#	output:
#		'Fastq/{tech}_R{Num,[12]}.fastq.gz'
#	shell:
#		"""
#		cat {input} > {output}
#		"""
#
rule coolerMapping:
	input:
		fastq_R1 = 'Fastq/{tech}_R1.fastq.gz',
		fastq_R2 = 'Fastq/{tech}_R2.fastq.gz',
		BWAIndex = BWAIndexPATH
	
	output:
		bam = "Bam/{tech}_cooler.bam"
		
	params:
		tech = techList,
		module1 = bwaVer,
		module2 = samtoolsVer
	
	shell:
		"""
		module purge && module load {params.module1} && module load {params.module2}
		bwa mem -SP5M {input.BWAIndex} {input.fastq_R1} {input.fastq_R2} | samtools view -bhS - > {output.bam}
		"""	


rule coolerParse:
	input:
		bam = "Bam/{tech}_cooler.bam",
		chromSize = chrsizePATH
	output:
		parsed = "Pairs/{tech}_parsed.pairsam"
	params:
		tech=techList,
		moduleVer = pairtoolsVer
	shell:
		"""
		module purge && module load {params.moduleVer}
		
		samtools view -h {input.bam} |  pairtools parse -c {input.chromSize} -o {output.parsed}
		"""

rule coolerSort:
	input:
		parsed= "Pairs/{tech}_parsed.pairsam"
	output:
		sorte="Pairs/{tech}_sorted.pairsam"
	params:
		moduleVer=pairtoolsVer
	shell:
		"""
		module purge && module load {params.moduleVer}
		pairtools sort --nproc 8 --tempdir=./ -o {output.sorte} {input.parsed}
		"""

rule coolerDedup:
	input:
		sorte="Pairs/{tech}_sorted.pairsam"
	output:
		dedup="Pairs/{tech}_dedup.pairsam"
	params:
		moduleVer=pairtoolsVer
	shell:
		"""
		module purge && module load {params.moduleVer}
		pairtools dedup --mark-dups -o {output.dedup} {input.sorte}
		"""

rule coolerFilter:
	input:
		dedup="Pairs/{tech}_dedup.pairsam"
	output:
		filtered="Pairs/{tech}_filtered.pairsam"
	params:
		moduleVer=pairtoolsVer
	shell:
		"""
		module purge && module load {params.moduleVer}
		pairtools select \ '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU)' \ -o {output.filtered} {input.dedup}
		"""

rule coolerSplit:
	input:
		filtered="Pairs/{tech}_filtered.pairsam"
	output:
		pairs="Pairs/{tech}.pairs"
	params:
		moduleVer=pairtoolsVer
	shell:
		"""
		module purge && module load {params.moduleVer}
		pairtools split --output-pairs {output.pairs} {input.filtered}
		"""

rule coolerBinning:
	input:
		chrsize = chrsizePATH,
		pairs = "Pairs/{tech}.pairs"
	output:
		cool = "Cool/{tech}.cool",
		mcool = "Cool/{tech}.mcool"
	params:
		moduleVer = coolerVer
	shell:
		"""
		module purge && module load {params.moduleVer}
		cooler cload pairix {input.chrsize}:1 {input.pairs} {output.cool}
		cooler balance {output.cool}
		cooler zoomify {output.cool}
		"""

rule hiceMapping:
	input:
		fastq = 'Fastq/{tech}_{readNum}.fastq.gz',
                BWAIndex = BWAIndexPATH

	output:
		bam='Bam/{tech}_mate{readNum}.bam',
		log='Logs/{tech}_mate{readNum}.log'
	params:
                module1 = bwaVer,
                module2 = samtoolsVer
	shell:
		"""
		module purge && module load {params.module1} {params.module2}
		bwa mem -A1 -B4  -E50 -L0  {input.BWAIndex} \
		    {input.fastq} 2>>{output.log} | samtools view -Shb - > {output.bam}
		"""

rule hiceBuildingMatrix:
	input:
		mate1 = 'Bam/{tech}_mateR1.bam',
		mate2 = 'Bam/{tech}_mateR2.bam',
		restriction = restrictionDigestPATH
	output:
		matrix='Mat/{tech}_hicMatrix.h5',
		qc=directory('QCfiles/{tech}_hicQC')
	params:
		moduleVer=hiceVer
	shell:
		"""
		module purge && module load {params.moduleVer}
		hicBuildMatrix -s {input.mate1} {input.mate2} -rs {input.restriction} -o {output.matrix} --QCfolder {output.qc}
		"""


