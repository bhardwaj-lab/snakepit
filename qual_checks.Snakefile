# @author: Vivek Bhardwaj (@vivekbhr)
# @date: Feb 15, 2017
# @desc: Snakemake pipeline for RNA editing detection following variant calling
#
# Usage: snakemake --snakefile qual_checks.Snakefile --jobs 2 --config [config params] -c "SlurmEasy -t {threads} -n {rule}"
## needs: fastqc, multiqc

from os.path import join
# Globals ---------------------------------------------------------------------

# Full path to output folder.
OUTPUT_DIR = config['outdir']


# A snakemake regular expression matching the forward mate FASTQ files.
SAMPLES, = glob_wildcards(join(OUTPUT_DIR, 'raw/{sample}_R1.fastq.gz'))
READS = ['R1','R2']
DIRS = ['raw','trimmed']
print(SAMPLES)

rule all:
    input:
        expand(join(OUTPUT_DIR, '{dir}' ,'fastqc', '{sample}_{read}_fastqc.zip'), dir = DIRS , sample = SAMPLES, read = READS),
        expand(join(OUTPUT_DIR, '{dir}' , 'fastqc' , 'multiqc_report.html'), dir = DIRS)

print(expand(join(OUTPUT_DIR, '{dir}' ,'fastqc', '{sample}_{read}_fastqc.zip'), dir = DIRS , sample = SAMPLES, read = READS))

rule fastqc:
    output: "{folder}/fastqc/{file_basename}_fastqc.zip"
    params: out="{folder}/fastqc"
    input: "{folder}/{file_basename}.fastq.gz"
    shell:
        'fastqc -o {params.out} {input}'

rule multiqc:
    output: "{folder}/fastqc/multiqc_report.html"
    input: "{folder}"
    params: out="{folder}/fastqc"
    shell:
        'multiqc -o {params.out} {input}'
