# conda create -n meme -c bioconda meme ucsc-bedgraphtobigwig snakemake
import re
import glob

genome_fasta='/hpc/hub_oudenaarden/vbhardwaj/annotations/mm10_gencode23/genome_and_annotation/GRCm38.p6.genome.fa'
chrsizes='/hpc/hub_oudenaarden/vbhardwaj/annotations/mm10_gencode23/STARindex/no_junctions/chrNameLength.txt'
motif_meme=glob.glob('*.meme')

motif_names=[re.sub("\.meme", "", x) for x in motif_meme]

print(motif_names)
rule all:
    input:
        expand("{motif}/fimo.bw", motif = motif_names),
        expand("{motif}/fimo.bed", motif = motif_names)

rule fimo:
    input:
        genome = genome_fasta,
        meme = "{motif}.meme"
    output: "{motif}/fimo.tsv"
    params:
        outdir = "{motif}"
    conda: "meme.yaml"
    shell:
        "fimo --max-stored-scores 10000000 -oc {params.outdir} {input.meme} {input.genome}"

rule fimo_bed:
    input: "{motif}/fimo.tsv"
    output:
        bed = "{motif}/fimo.bed",
        bg = "{motif}/fimo.bg"
    conda: "meme.yaml"
    shell:
        """
        awk 'OFS="\\t" {{ print $3, $4, $5, $2, $7, $6 }}' {input} > {output.bed}; \
        awk 'OFS="\\t" {{ print $3, $4, $5, $7 }}' {input} > {output.bg}
        """

rule fimo_bw:
    input:
        bg = "{motif}/fimo.bg",
        sizes = chrsizes
    output: "{motif}/fimo.bw"
    conda: "meme.yaml"
    shell:
        "bedGraphToBigwig {input.bg} {input.sizes} {output}"
