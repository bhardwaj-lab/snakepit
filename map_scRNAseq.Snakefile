# @author: Vivek Bhardwaj (@vivekbhr)
# @desc: Snakemake workflow for creating pseudo-bulk coverage bigwigs from bam files and a tsv file with cell->cluster information
#
# Usage: snakemake -s map_scRNAseq.Snakefile --cores 20 --jobs 4 --config [config params] -c "SlurmEasy -t {threads} -n {rule}"
## needs: STAR > 1.7, deeptools, samtools
samples = config['samplenames'].split(,)
BCwhiteList = config['barcodes']
tempDir = config['tempdir']#/hpc/hub_oudenaarden/vbhardwaj/tempdir"
GTF = config['gtf']#"/hpc/hub_oudenaarden/vbhardwaj/annotations/mm10_gencode23/genome_and_annotation/gencode.vM23.annotation.gtf"
star_index = config['star_index']#"/hpc/hub_oudenaarden/vbhardwaj/annotations/mm10_gencode23/STARindex/no_junctions/"
STARsoloCoords = [9, 6, 1, 8] if config['protocol'] == "vasa" else [17, 12, 1, 16] ## UMI start, UMI length, Cell BC start, cell BC length (1-based)

trim_dir = "FASTQ_trimmed"
rule all:
    input:
        fileList = [expand("STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/barcodes.tsv", sample = samples),
                    expand("bigwigs/{sample}.bw", sample = samples)]

# mapped using STARsolo so possible to split SAM files per cell type later on
rule mapReads:
    input:
        read1 = trim_dir + "/{sample}_R1.fastq.gz",
        read2 = trim_dir + "/{sample}_R2.fastq.gz"
    output:
        bam = "STARsolo/{sample}.sorted.bam",
        raw_counts = "STARsolo/{sample}/{sample}.Solo.out/Gene/raw/matrix.mtx",
        filtered_counts = "STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/matrix.mtx",
        filtered_bc = "STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/barcodes.tsv"
    params:
        gtf = GTF,
        index = star_index,
        prefix = "STARsolo/{sample}/{sample}.",
        samsort_memory = '2G',
        sample_dir = "STARsolo/{sample}",
        bclist = BCwhiteList,
        UMIstart = STARsoloCoords[0],
        UMIlen = STARsoloCoords[1],
        CBstart = STARsoloCoords[2],
        CBlen = STARsoloCoords[3],
        outdir = "STARsolo/",
        tempDir = tempDir,
        sample = "{sample}"
    log: "logs/mapReads_{sample}.err"
    threads: 20
    shell:
        """
        TMPDIR={params.tempDir}
        MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
        ( [ -d {params.sample_dir} ] || mkdir -p {params.sample_dir} )
        STAR --runThreadN {threads} \
            --sjdbOverhang 100 \
            --outSAMunmapped Within \
            --outSAMtype BAM SortedByCoordinate \
            --outBAMsortingBinsN 20 \
            --outSAMattributes NH HI AS nM CB UB \
            --sjdbGTFfile {params.gtf} \
            --genomeDir {params.index} \
            --readFilesIn {input.read2} {input.read1}\
            --readFilesCommand zcat \
            --outFileNamePrefix {params.prefix} \
        --soloFeatures Gene Velocyto \
        --soloUMIstart {params.UMIstart} \
        --soloUMIlen {params.UMIlen} \
        --soloCBstart {params.CBstart} \
        --soloCBlen {params.CBlen} \
        --soloCBwhitelist {params.bclist} \
        --soloBarcodeReadLength 0 \
        --soloStrand Forward\
        --soloCBmatchWLtype Exact \
        --soloType CB_UMI_Simple \
        --soloUMIdedup Exact 2> {log}
        ln -r -s {params.outdir}/{params.sample}/{params.sample}.Aligned.sortedByCoord.out.bam {output.bam} 2>> {log}
        rm -rf $MYTEMP
        """

rule indexBAM:
    input: "STARsolo/{sample}.sorted.bam"
    output: "STARsolo/{sample}.bam.bai"
    log: "logs/indexBAM.{sample}.log"
    threads: 1
    shell: 'samtools index {input} > {output} 2>> {log}'

# Get bigwig files of the trimmed fastq files check for A/T stretches
rule getBW:
    input:
        bam = "STARsolo/{sample}.sorted.bam",
        idx = "STARsolo/{sample}.bam.bai"
    output: "bigwigs/{sample}.bw"
    log: "logs/getBW_{sample}.err"
    threads: 20
    shell:
        "bamCoverage --minMappingQuality 10 --samFlagExclude 256 -p {threads} --normalizeUsing CPM -b {input.bam} -o {output} 2> {log}"
