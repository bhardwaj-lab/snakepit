# @author: Vivek Bhardwaj (@vivekbhr)
# @date: Feb 15, 2017
# @desc: Snakemake pipeline for RNA editing detection following variant calling
#
# Usage: snakemake --snakefile qual_checks.Snakefile --jobs 2 -c "SlurmEasy -t {threads} -n {rule}"

## needs: trimgalore, fastqc, multiqc
from os.path import join
# Globals ---------------------------------------------------------------------

# Full path to output folder.
OUTPUT_DIR = config['workdir']#"01_fastq"


# A snakemake regular expression matching the forward mate FASTQ files.
SAMPLES, = glob_wildcards(join(OUTPUT_DIR, 'fastq/{sample}_R1.fastq.gz'))
READS = ['R1','R2']
print(SAMPLES)

rule all:
    input:
        expand(join(OUTPUT_DIR, 'qual_trimmed', '{sample}_{read}.fastq.gz'),sample = SAMPLES, read = READS),
        expand(join(OUTPUT_DIR, 'qual_trimmed/fastqc', '{sample}_{read}_fastqc.zip'),sample = SAMPLES, read = READS),
        expand(join(OUTPUT_DIR, 'qual_trimmed/fastqc', 'multiqc_report.html'),sample = SAMPLES)

print(expand(join(OUTPUT_DIR, 'qual_trimmed', '{sample}_{read}.fastq.gz'),sample = SAMPLES, read = READS))
#print(expand(join(OUTPUT_DIR, 'qual_trimmed/fastqc', '{sample}_{read}_fastqc.zip'),sample = SAMPLES, read = READS))

rule trim:
    output:
        R1 = "{file_basename}_R1_val_1.fq.gz",
        R2 = "{file_basename}_R2_val_2.fq.gz"
    params:
        out="qual_trimmed"
    input: "{file_basename}_R1.fastq.gz","{file_basename}_R2.fastq.gz"
    shell:
        'trim_galore --paired -q 20 -o {params.out} {input}'

rule rename:
    output:
        R1 = "{file_basename}_R1.fastq.gz",
        R2 = "{file_basename}_R2.fastq.gz"
    input:
        R1 = "{file_basename}_R1_val_1.fq.gz",
        R2 = "{file_basename}_R2_val_2.fq.gz"
    shell:
        'mv {input.R1} {output.R1}; mv {input.R2} {output.R2}'

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
