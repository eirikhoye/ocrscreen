"""
Module bwa.

Per sequencing run: bwa-mem2 alignment in Hi-C mode (-SP5M) streamed directly into
pairtools parse (read-name ordered SAM, as pairtools requires) and pairtools sort.
Follows the 4DN Hi-C processing pipeline / Open2C distiller layout. No BAM is kept.
"""


rule bwa_mem2_index:
    input:
        GENOME_FASTA,
    output:
        BWA_INDEX,
    params:
        prefix=lambda wc, output: output[0].removesuffix(".0123"),
    resources:
        mem_mb=96000,
    log:
        os.path.join(LOGS, "resources", "bwa_mem2_index.log"),
    conda:
        "../envs/hic.yaml"
    shell:
        "bwa-mem2 index -p {params.prefix} {input} > {log} 2>&1"


rule bwa_parse_sort:
    input:
        fastq1=FASTQ_1,
        fastq2=FASTQ_2,
        index=BWA_INDEX,
        chrom_sizes=CHROM_SIZES,
    output:
        pairs=RUN_PAIRS,
        stats=RUN_PARSE_STATS,
    params:
        prefix=lambda wc, input: input.index[0].removesuffix(".0123"),
        parse=config["hic"]["parse_options"],
        assembly=BUILD,
        sort_threads=config["hic"]["threads"]["sort"],
        tmpdir=lambda wc, output: os.path.join(os.path.dirname(output.pairs), f"tmp_{wc.run}"),
    threads: config["hic"]["threads"]["bwa"]
    log:
        bwa=os.path.join(LOGS, "{cell_line}", "hic", "bwa_{run}.log"),
        parse=os.path.join(LOGS, "{cell_line}", "hic", "parse_{run}.log"),
        sort=os.path.join(LOGS, "{cell_line}", "hic", "sort_{run}.log"),
    conda:
        "../envs/hic.yaml"
    shell:
        """
        mkdir -p {params.tmpdir}
        bwa-mem2 mem -SP5M -t {threads} {params.prefix} {input.fastq1} {input.fastq2} 2> {log.bwa} | \
            pairtools parse {params.parse} -c {input.chrom_sizes} --assembly {params.assembly} \
                --output-stats {output.stats} 2> {log.parse} | \
            pairtools sort --nproc {params.sort_threads} --tmpdir {params.tmpdir} \
                -o {output.pairs} 2> {log.sort}
        rm -rf {params.tmpdir}
        """
