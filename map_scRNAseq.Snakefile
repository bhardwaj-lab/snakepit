# @author: Vivek Bhardwaj (@vivekbhr)
# @desc: Snakemake workflow for mapping VASA or 10x samples using STAR
# setup: mamba create -n rnaseq_tchic -c bioconda -c conda-forge fastqc cutadapt star multiqc deeptools samtools snakemake
# Usage: snakemake -s map_scRNAseq.Snakefile --cores 20 --jobs 4 --config [config params] -c "SlurmEasy -t {threads} -n {rule}"
## needs: STAR > 2.7, deeptools, samtools
## for both VASA and 10x data, R2 has the sequence and R1 has the barcode, therefore can be mapped with the same parameters
import os
import glob

BCwhiteList = config['barcodes']
tempDir = config['tempdir']#/hpc/hub_oudenaarden/vbhardwaj/tempdir"
trim = config['trim']
GTF = config['gtf']#"/hpc/hub_oudenaarden/vbhardwaj/annotations/mm10_gencode23/genome_and_annotation/gencode.vM23.annotation.gtf"
star_index = config['star_index']#"/hpc/hub_oudenaarden/vbhardwaj/annotations/mm10_gencode23/STARindex/no_junctions/"
STARsoloCoords = [1, 6, 7, 8] if config['protocol'] == "vasa" else [17, 12, 1, 16] ## UMI start, UMI length, Cell BC start, cell BC length (1-based)

# optional args
try:
    indir=config['indir']
except KeyError:
    indir='FASTQ'
# assuming files are named _R1/R2.fastq.gz, get sample names
try:
    samples = config['samples'].strip().split(',')
except KeyError:
    infiles = sorted(glob.glob(os.path.join(indir, '*_R1.fastq.gz')))
    samples = [os.path.basename(x)[:-len('_R1.fastq.gz')] for x in infiles]

# organism type for max intron length
try:
    organism = config['organism']
except KeyError:
    organism = None

if trim:
    trim_dir = "FASTQ_trimmed"
    other=expand("FASTQ_trimmed/FastQC/{sample}{read}_fastqc.html", sample = samples, read=['_R1', '_R2'])
else:
    trim_dir = indir
    other=[]

rule all:
    input:
        fileList = [expand("STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/matrix.mtx", sample = samples),
                    expand("bigwigs/{sample}.bw", sample = samples),
                    "QC/multiqc_report.html",
                    other]

if trim:
    rule cutadapt:
        input:
            R1 = os.path.join(indir,"{sample}_R1.fastq.gz"),
            R2 = os.path.join(indir,"{sample}_R2.fastq.gz")
        output:
            R1 = temp("FASTQ_trimmed/{sample}_R1.fastq.gz"),
            R2 = temp("FASTQ_trimmed/{sample}_R2.fastq.gz"),
            QC = "QC/cutadapt/{sample}.out"
        params:
            opts = "-A W{'10'}"
        threads: 10
        resources:
            mem_mb=50000
        shell:
            """
            cutadapt -j {threads} {params.opts} -e 0.1 -q 16 -O 3 --trim-n --minimum-length 10 \
            -a AGATCGGAAGAGC -A AGATCGGAAGAGC --nextseq-trim=16 \
            -b TGGAATTCTCGGGTGCCAAGG -B TGGAATTCTCGGGTGCCAAGG \
            -b ATCTCGTATGCCGTCTTCTGCTTG -B ATCTCGTATGCCGTCTTCTGCTTG \
            -b GTTCAGAGTTCTACAGTCCGACGATC -B GTTCAGAGTTCTACAGTCCGACGATC \
            -o "{output.R1}" -p "{output.R2}" "{input.R1}" "{input.R2}" > {output.QC}
            """

    rule FastQC_trimmed:
        input:
            "FASTQ_trimmed/{sample}{read}.fastq.gz"
        output:
            "FASTQ_trimmed/FastQC/{sample}{read}_fastqc.html"
        params:
            outdir = "FASTQ_trimmed/FastQC"
        log:
            out = "logs/FastQC_trimmed.{sample}{read}.out",
            err = "logs/FastQC_trimmed.{sample}{read}.err"
        threads: 2
        shell: "fastqc -o {params.outdir} {input} > {log.out} 2> {log.err}"

# mapped using STARsolo so possible to split SAM files per cell type later on
rule mapReads:
    input:
        read1 = trim_dir + "/{sample}_R1.fastq.gz",
        read2 = trim_dir + "/{sample}_R2.fastq.gz"
    output:
        bam = "STARsolo/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        raw_counts = "STARsolo/{sample}/{sample}.Solo.out/Gene/raw/matrix.mtx",
        filtered_counts = "STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/matrix.mtx",
        filtered_bc = "STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/barcodes.tsv",
        tmpdir=temp(directory("STARsolo/{sample}/{sample}._STARtmp")),
        tmpgenome=temp(directory("STARsolo/{sample}/{sample}._STARgenome"))
    params:
        star_bugfix='/hpc/hub_oudenaarden/vbhardwaj/programs/STAR-2.7.10a_alpha_220506/source',# use this unil version 2.7.11 is released, to avoid segfault
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
    resources:
        mem_mb=100000
    shell:
        """
    TMPDIR={params.tempDir}
    MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
    ( [ -d {params.sample_dir} ] || mkdir -p {params.sample_dir} )
    maxIntronLen=`expr $(awk '{{ print $3-$2 }}' {params.index}/sjdbList.fromGTF.out.tab | sort -n -r | head -1) + 1`
    {params.star_bugfix}/STAR --runThreadN {threads} \
    --sjdbOverhang 100 \
    --outSAMunmapped Within \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMsortingBinsN 20 \
    --genomeDir {params.index} \
    --readFilesIn {input.read2} {input.read1}\
    --readFilesCommand zcat \
    --outFileNamePrefix {params.prefix} \
    --outSAMattributes NH HI AS nM MD jM jI MC ch CB UB GX GN \
    --outSAMstrandField intronMotif \
    --sjdbGTFfile {params.gtf} \
    --outFilterType BySJout \
    --outFilterIntronMotifs RemoveNoncanonical \
    --alignIntronMax ${{maxIntronLen}} \
    --alignSJoverhangMin 8 \
    --soloFeatures Gene Velocyto \
    --soloCBstart {params.CBstart} \
    --soloCBlen {params.CBlen} \
    --soloUMIstart {params.UMIstart} \
    --soloUMIlen {params.UMIlen} \
    --soloCBwhitelist {params.bclist} \
    --soloBarcodeReadLength 0 \
    --soloStrand Forward \
    --soloCBmatchWLtype Exact \
    --soloType CB_UMI_Simple \
    --soloUMIdedup NoDedup 2> {log}
    rm -rf $MYTEMP
        """

#    --soloUMIdedup Exact
rule filterBAMunique:
    input: "STARsolo/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
    output: temp("STARsolo/{sample}.uniqueReads.sam")
    threads: 10
    resources:
        mem_mb=50000
    shell:
        "samtools view -@ {threads} -h -d CB -F 4 -F 256 -F 2048 -q 255 {input} | grep -w -v 'CB:Z:-' > {output}"

rule dedupBAMunique:
    input: "STARsolo/{sample}.uniqueReads.sam",
    output:
        bam="STARsolo/{sample}.uniqueReads.bam",
        qc="QC/umi_dedup/{sample}_per_umi_per_position.tsv"
    params:
        sample = "{sample}",
        tempdir= tempDir
    log: "logs/filterBAM.{sample}.log"
    threads: 1
    resources:
        mem_mb=80000
    shell:
        """
        umi_tools dedup --per-cell --cell-tag CB --umi-tag UB --extract-umi-method tag \
        --method unique --spliced-is-unique \
        --output-stats=QC/umi_dedup/{params.sample} \
        --temp-dir {params.tempdir} -L {log} -i -I {input} > {output.bam}
        """

rule indexBAM:
    input: "STARsolo/{sample}.uniqueReads.bam"
    output: "STARsolo/{sample}.uniqueReads.bam.bai"
    log: "logs/indexBAM.{sample}.log"
    threads: 1
    shell: 'samtools index {input}'

# Get bigwig files of the trimmed fastq files check for A/T stretches
rule getBW:
    input:
        bam = "STARsolo/{sample}.uniqueReads.bam",
        idx = "STARsolo/{sample}.uniqueReads.bam.bai"
    output: "bigwigs/{sample}.bw"
    log: "logs/getBW_{sample}.err"
    threads: 20
    shell:
        "bamCoverage -p {threads} --normalizeUsing CPM -b {input.bam} -o {output} > {log} 2>&1"

rule multiQC:
    input: expand("STARsolo/{sample}.uniqueReads.bam", sample = samples)
    output: "QC/multiqc_report.html"
    params:
        outdir = "QC"
    log: "logs/multiqc.log"
    threads: 1
    shell:
        "multiqc -f -o {params.outdir} . > {log} 2>&1"
