# @author: Vivek Bhardwaj (@vivekbhr)
# @date: May 22 2017
#
# Usage: snakemake --snakefile Snakefile.fwd-rev --cores 20 --jobs 3 -c "SlurmEasy -t {threads} -n {rule}"



from os.path import join, dirname

# step 1 . bamcoverage +/- strand for each file
# 2. bigwigcompare (add) the replicates
# 3. bigwigcompare (subrtact) the + file with - file for each sample



OUTDIR = "fwd-rev"
INDIR = <indir>
THREADS = 10

SAMPLES,REPS, = glob_wildcards(join(INDIR,'{sample,[^/]+}{rep,[ab]}_dupFilt.bam'))

STRAND = ['forward', 'reverse']

if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

rule all :
    input : expand(join(OUTDIR, '{sample}{rep}.{strand}.bw'), sample = set(SAMPLES), rep = set(REPS), strand = STRAND),
            expand(join(OUTDIR, '{sample}.{strand}_merge.bw'), sample = SAMPLES, strand = STRAND),
            expand(join(OUTDIR, '{sample}.fwd-rev.bw'), sample = SAMPLES)

rule bamcoverage:
    input : join(INDIR, '{basename}{rep}_dupFilt.bam')
    output :
        join(OUTDIR, '{basename}{rep}.{strand}.bw')
    threads : THREADS
    params :
        strand = '{strand}'
    shell:
        '/package/deeptools-2.5.1/bin/bamCoverage -bs 1 -p {threads} --Offset 1 5 --normalizeUsingRPKM '
        '--filterRNAstrand {params.strand} -b {input} -o {output}'

rule mergereps:
    input : '{sample}a.{strand}.bw', '{sample}b.{strand}.bw'
    output : '{sample}.{strand}_merge.bw'
    threads : THREADS
    shell:
        '/package/deeptools-2.5.1/bin/bigwigCompare -bs 1 -p {threads}  --ratio mean -b1 {input[0]} -b2 {input[1]} -o {output}'

rule mergestrands:
    input : '{sample}.forward_merge.bw', '{sample}.reverse_merge.bw'
    output : '{sample}.fwd-rev.bw'
    threads : THREADS
    shell:
        '/package/deeptools-2.5.1/bin/bigwigCompare -bs 1 -p {threads}  --ratio subtract -b1 {input[0]} -b2 {input[1]} -o {output}'
