# @author: Vivek Bhardwaj (@vivekbhr)
# @desc: Snakemake workflow for creating pseudo-bulk coverage bigwigs from bam files and a tsv file with cell->cluster information
#
# Usage: snakemake -s pseudoBulk_coverage.Snakefile --cores 20 --jobs 4 --config [config params] -c "SlurmEasy -t {threads} -n {rule}"
## needs: deeptools, samtools
import csv
import os
infolder = config['bamfolder']
samples = config['samples']
clusters_tsv = config['clusters']
splitbam_path = config['splitbam']
bctag = config['barcodetag']
deeptoolsArgs = config['deeptoolsArgs']
clusterRegex = config['clusterRegex']
# -bl {params.blacklist} -ignore X Y M "
# "--effectiveGenomeSize 2652783500 "

print(config)
def getClusters(tsv):
    cluster_dict = {}
    with open(tsv) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            cluster_dict[row[0]] = row[1]
    clusters = set(x for x in cluster_dict.values())
    return(list(clusters))

clusters = getClusters(clusters_tsv[0])

rule all:
    input: expand("clusters_split/{sample}_{cluster}.bam", sample = samples, cluster = clusters),
           expand("clusters_split/{cl}_merged.bam", cl = clusters),
           expand("clusters_split/{cl}_merged.bam.bai", cl = clusters),
           expand("clusters_split/{cl}_merged.bw", cl = clusters)
print(samples)
print(clusters)
## split bam by cluster
rule splitBAM:
    input:
        bam = infolder+"/{sample}.bam",
        bai = infolder+"/{sample}.bam.bai",
        clusters = '{sample}.txt'
    output:
        bams = expand("clusters_split/{{sample}}_{cluster}.bam", cluster = clusters)
    wildcard_constraints:
        cluster=clusterRegex
    params:
        outdir = "clusters_split/{sample}_",
        splitbam = splitbam_path,
        bctag = bctag
    threads: 1
    shell:
        "{params.splitbam} -i {input.bam} -t {params.bctag} -c {input.clusters} -o {params.outdir}"

## index bam
rule index:
    input: "clusters_split/{sample}_{cluster}.bam"
    output: temp("clusters_split/{sample}_{cluster}.bam.bai")
    wildcard_constraints:
        cluster=clusterRegex
    threads: 1
    shell:
        "samtools index {input}"

## make bigwig
rule bamcov:
    input:
        bam = "clusters_split/{sample}_{cluster}.bam",
        bai = "clusters_split/{sample}_{cluster}.bam.bai"
    output: "clusters_split/{sample}_{cluster}.bw"
    wildcard_constraints:
        cluster=clusterRegex
    params:
        args = deeptoolsArgs
    threads: 20
    shell:
        "bamCoverage -p {threads} --skipNAs "
        "{params.args} -b {input.bam} -o {output}"

## optionally, merge bw across samples
rule mergeBam:
    input:
         expand("clusters_split/{sample}_{{cl}}.bam", sample = samples)
    output: "clusters_split/{cl}_merged.bam"
    wildcard_constraints:
        cl=clusterRegex
    threads: 10
    shell:
        "samtools merge -@ {threads} {output} {input}"

rule idx_mergeBam:
    input: "clusters_split/{cl}_merged.bam"
    output: "clusters_split/{cl}_merged.bam.bai"
    wildcard_constraints:
        cl=clusterRegex
    threads: 1
    shell:
        "samtools index {input}"

rule bamcov_mergeBam:
    input:
        bam = "clusters_split/{cl}_merged.bam",
        bai = "clusters_split/{cl}_merged.bam.bai"
    output: "clusters_split/{cl}_merged.bw"
    wildcard_constraints:
        cl=clusterRegex
    params:
        args = deeptoolsArgs
    threads: 20
    log: "{cl}_bamCoverage.log"
    shell:
        "bamCoverage -p {threads} --skipNAs "
        "{params.args} -b {input.bam} -o {output} > {log} 2>&1"
