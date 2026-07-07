
## needs: meme, ucsc-bedGraphToBigWig, bedtools
# conda create -n meme -c bioconda meme ucsc-bedgraphtobigwig bedtools snakemake snakemake-executor-plugin-slurm
import re
import glob

genome_fasta=config['fasta']#/hpc/hub_oudenaarden/vbhardwaj/annotations/mm10_gencode23/genome_and_annotation/GRCm38.p6.genome.fa'
chrsizes=config['chrsizes']#/hpc/hub_oudenaarden/vbhardwaj/annotations/mm10_gencode23/STARindex/no_junctions/chrNameLength.txt'
jaspar_memefile=config['jaspar']# motif file obtained from jaspar; with all motifs in .MEME format (https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_non-redundant_pfms_meme.txt)
motif_names=config['motifs'].split(",") # motif names separated by comma

# optional: given a bed file, only use those regions from genome_fasta
try:
    bedfile=config['bedfile']
except KeyError:
    bedfile=None

#motif_meme=glob.glob('*.meme')
#motif_names=[re.sub("\.meme", "", x) for x in motif_meme]

print(motif_names)
rule all:
    input:
        expand("{motif}/fimo.bw", motif = motif_names),
        expand("{motif}/fimo.bed", motif = motif_names)

rule get_motif:
    input: jaspar_memefile
    output: "{motif}/{motif}.meme"
    params:
        motif_name = "{motif}"
    #conda: "meme.yaml"
    shell:
        """
        sed -n '1,/^MOTIF/p' {input} | sed '$d' > {output} && \
        sed -En "/(^|[[:space:]]){params.motif_name}([[:space:]]|$)/,/^$/p" {input} >> {output}
        """

rule fimo:
    input:
        genome = genome_fasta,
        meme = "{motif}/{motif}.meme"
    output: "{motif}/fimo.tsv"
    params:
        outdir = "{motif}"
    #conda: "meme.yaml"
    shell:
        "fimo --max-stored-scores 10000000 -oc {params.outdir} {input.meme} {input.genome}"

rule fimo_bed:
    input: "{motif}/fimo.tsv"
    output:
        bed = "{motif}/fimo.bed",
        bg = temp("{motif}/fimo.bg")
    #conda: "meme.yaml"
    shell:
        """
        awk 'OFS="\\t" {{ if(NR>1) {{print $3, $4, $5, $2, $7, $6}} }}' {input} | head -n -4 | \
        sort -k1,1 -k2,2n - > {output.bed} && \
        awk 'OFS="\\t" {{ if(NR>1) {{print $3, $4, $5, $7}} }}' {input} | head -n -4 | \
        sort -k1,1 -k2,2n | bedtools merge -i - -c 4 -o sum > {output.bg}
        """

rule fimo_bw:
    input:
        bg = "{motif}/fimo.bg",
        sizes = chrsizes
    output: "{motif}/fimo.bw"
    #conda: "meme.yaml"
    shell:
        " bedGraphToBigWig {input.bg} {input.sizes} {output}"
