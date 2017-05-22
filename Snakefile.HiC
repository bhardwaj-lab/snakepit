# bwa mem mapping options:
#       -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [1]
#       -B INT        penalty for a mismatch [4]
#       -O INT[,INT]  gap open penalties for deletions and insertions [6,6]
#       -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [1,1] # this is set very high to avoid gaps 
#                                  at restriction sites. Setting the gap extension penalty high, produces better results as
#                                  the sequences left and right of a restriction site are mapped independently.
#       -L INT[,INT]  penalty for 5'- and 3'-end clipping [5,5] # this is to to no penalty.

MAPPER_CMD = "/package/bwa-0.7.12/bin/bwa mem -A1 -B4  -E50 -L0 -t 10"
SAMTOOLS = "/package/samtools-1.2/bin/samtools"

HICEXPLORER_PATH = "/data/manke/group/ramirez/tools/HiCExplorer_v0.2/bin"

import os

os.environ['PYTHONPATH'] = '/data/manke/group/ramirez/tools/HiCExplorer_v0.2/'

rule all:
	input: "hic/merge.h5"

rule map_fastq_single_end:
	input: "fastq/{sample}.fastq.gz"
	output: "bam/{sample}.bam"
	log:    "bam/{sample}.log"
	shell: 
		"""
		echo "mapping {input}" > {log}
		{MAPPER_CMD} {INDEX} {input}  2>> {log} |
		{SAMTOOLS} view -Shb - > {output}
		"""

rule get_HindIII_bed:
	output: "HindIII.bed"
	shell: "{HICEXPLORER_PATH}/findRestSite -f {GENOME} --searchPattern AAGCTT -o {output}"

rule get_DpnII_bed:
	output: "DpnII.bed"
	shell: "{HICEXPLORER_PATH}/findRestSite -f {GENOME} --searchPattern GATC -o {output}"

rule make_hic_matrix:
	input: "bam/{sample}_1.bam", "bam/{sample}_2.bam", RESTRICTION_FILE
	output: "hic/{sample}.h5", bam="bam/{sample}_R12.bam"
	log: "hic/{sample}.log"
	shell:
		"""
	        {HICEXPLORER_PATH}/hicBuildMatrix -s {input[0]} {input[1]} -rs {input[2]} \
	        --restrictionSequence {RESTRICTION_SEQUENCE} \
	        --minDistance {MIN_RS_DISTANCE} \
	        --maxDistance {MAX_RS_DISTANCE} \
	        -b {output.bam} -o {output[0]} &> {log} 
		"""
rule sort_bam:
	input: "bam/{sample}.bam"
	output: "bam/{sample}_sorted.bam"
	params: prefix="bam/_{sample}"
	shell:
		"""
		{SAMTOOLS} sort {input} -@20 -T {params.prefix} -O bam -o {output[0]}
		"""

rule index_bam:
	input: "bam/{sample}.bam"
	output: "bam/{sample}.bam.bai"
	shell: "{SAMTOOLS} index {input}"


#rule download:
#	output: "fastq/{sample}_1.fastq.gz", "fastq/{sample}_2.fastq.gz" 
#	params: sra="{sample}"
#	shell: "/package/sratoolkit.2.5.0/bin/fastq-dump.2.5.0 --split-files --gzip --outdir fastq {params.sra}"

def ftp_from_sra(sra):
	if sra.endswith("_1") or sra.endswith("_2"):
	   mate = sra[-2:]
	   sra = sra[:-2]
	else:
	  mate = ''
	print("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{0}/00{1}/{2}/{2}{3}.fastq.gz".format(sra[0:6], sra[-1], sra, mate))
	return "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{0}/00{1}/{2}/{2}{3}.fastq.gz".format(sra[0:6], sra[-1], sra, mate)
        # example:  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR165/001/SRR1658531/SRR1658531_1.fastq.gz
	#           ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR165/009/SRR1658529/SRR1658529_1.fastq.gz

rule download_ena:
	params: sra ="{sample}"
	output: "fastq/{sample}.fastq.gz"
	run:
	  ftp = ftp_from_sra(params.sra)
	  shell("wget -q -P fastq {}".format(ftp))

rule merge:
	input: expand("hic/{sample}.h5", sample=SAMPLES)
	output: "hic/merge.h5"
	shell: "{HICEXPLORER_PATH}/hicSumMatrices -m {input} -o {output}"


