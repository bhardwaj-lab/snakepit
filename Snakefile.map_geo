# @author: Vivek Bhardwaj (@vivekbhr)
# @date: October 31, 2016
# @desc: Snakemake pipeline for https://www.broadinstitute.org/gatk/guide/article?id=3891
#        Based off of code by https://github.com/slowkow/snakefiles/blob/master/star/star.snakefile
#
# Usage: snakemake --snakefile map_geo.Snakefile --cores 20 --jobs 4 -c "SlurmEasy -t {threads} -n {rule}"


from os.path import join, dirname
from subprocess import check_output


# Globals ---------------------------------------------------------------------
THREADS=16

# Full path to genome fasta.
STARIDX = <starindex_dir>
# Full path to gene model annotations for splice aware alignment.
GTF = <annotation.gtf>
# Full path to a folder that holds all of your FASTQ files.
FASTQ_DIR = "01_fastq"
# Full path to output folder.
OUTPUT_DIR = "02_mapping"

# A snakemake regular expression matching the forward mate FASTQ files.
SAMPLES, = glob_wildcards(join(FASTQ_DIR, '{sample,[^/]+}_1.fastq'))

print(SAMPLES)
# Patterns for the 1st mate and the 2nd mate using the 'sample' wildcard.
PATTERN_R1 = '{sample}_1.fastq'
PATTERN_R2 = '{sample}_2.fastq'

#rules -----------------------------------------------------------------------

rule all:
    input:
        expand(join(OUTPUT_DIR, '{sample}', 'Aligned.out.bam'), sample = SAMPLES)

# 1. Map paired-end RNA-seq reads to the genome.
# 2. Make a coordinate sorted BAM with genomic coordinates.
# 3. Count the number of reads mapped to each gene.
# 4. Count the number of reads supporting each splice junction.
rule star:
    input:
        r1 = join(FASTQ_DIR, PATTERN_R1),
        r2 = join(FASTQ_DIR, PATTERN_R2),
        gtf = GTF,
        idx = STARIDX
    output:
        sam = join(OUTPUT_DIR, '{sample}', 'Aligned.out.sam'),
        counts = join(OUTPUT_DIR, '{sample}', 'ReadsPerGene.out.tab'),
        sj = join(OUTPUT_DIR, '{sample}', 'SJ.out.tab')
    log:
        join(OUTPUT_DIR, '{sample}', 'star.map.log')
    threads:
        THREADS
    run:
        # Map reads with STAR.
        shell('/package/STAR-2.5.2b/bin/STAR'
              ' --runThreadN {threads}'
              ' --genomeDir {input.idx}'
              ' --sjdbGTFfile {input.gtf}'
              ' --readFilesIn {input.r1} {input.r2}'
              # BAM file in transcript coords, in addition to genomic BAM file.
              ' --quantMode GeneCounts'
              # Basic 2-pass mapping, with all 1st pass junctions inserted
              # into the genome indices on the fly.
              ' --twopassMode Basic'
              # By default, this prefix is "./".
              ' --outFileNamePrefix ' + join(OUTPUT_DIR, '{wildcards.sample}') + '/'
              # If exceeded, the read is considered unmapped.
              ' --outFilterMultimapNmax 20'
              # Minimum overhang for unannotated junctions.
              ' --alignSJoverhangMin 8'
              # Minimum overhang for annotated junctions.
              ' --alignSJDBoverhangMin 1'
              # Minimum intron length.
              ' --alignIntronMin 1'
              # Maximum intron length.
              ' --alignIntronMax 1000000'
              # Maximum genomic distance between mates.
              ' --alignMatesGapMax 1000000'
              ' > {log} 2>&1')

# add read groups, sort
rule picard_cleanstar:
    input:
        sam =rules.star.output.sam
    output:
        bam = join(OUTPUT_DIR, '{sample}', 'Aligned.out.bam'),
    log:
        join(OUTPUT_DIR, '{sample}', 'picard.clean.log')
    run:
        shell('module load picard-tools/2.3.0;'
              'java -jar /package/picard-tools-2.3.0/picard.jar'
              ' AddOrReplaceReadGroups'
              ' I={input.sam}'
              ' O={output}'
              ' SORT_ORDER=coordinate CREATE_INDEX=false RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample'
              ' > {log} 2>&1')
