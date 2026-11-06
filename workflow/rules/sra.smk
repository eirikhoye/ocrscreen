"""
Module sra.

Download each Hi-C sequencing run listed in config/hic_samples.tsv from SRA and
convert to gzipped paired FASTQ. Set hic.max_spots to subsample (testing only).
"""


rule download_sra:
    output:
        temp(SRA_FILE),
    params:
        max_size=config["hic"].get("sra_max_size", "200g"),
    log:
        os.path.join(LOGS, "{cell_line}", "hic", "prefetch_{run}.log"),
    conda:
        "../envs/sra.yaml"
    shell:
        "prefetch {wildcards.run} --output-file {output} --max-size {params.max_size} > {log} 2>&1"


if config["hic"].get("max_spots"):

    rule fastq_dump_subsample:
        """Subsample the first hic.max_spots read pairs (was hic_abc sra_subsample.smk)."""
        output:
            fastq1=FASTQ_1,
            fastq2=FASTQ_2,
        params:
            outdir=lambda wc, output: os.path.dirname(output.fastq1),
            max_spots=config["hic"]["max_spots"],
        log:
            os.path.join(LOGS, "{cell_line}", "hic", "fastq_dump_{run}.log"),
        conda:
            "../envs/sra.yaml"
        shell:
            """
            fastq-dump {wildcards.run} --gzip --split-3 --skip-technical \
                --minSpotId 1 --maxSpotId {params.max_spots} -O {params.outdir} > {log} 2>&1
            """

else:

    rule fasterq_dump:
        input:
            sra=SRA_FILE,
        output:
            fastq1=FASTQ_1,
            fastq2=FASTQ_2,
        params:
            outdir=lambda wc, output: os.path.dirname(output.fastq1),
            tmp=lambda wc, output: os.path.join(os.path.dirname(output.fastq1), f"tmp_{wc.run}"),
        threads: 16
        log:
            os.path.join(LOGS, "{cell_line}", "hic", "fasterq_dump_{run}.log"),
        conda:
            "../envs/sra.yaml"
        shell:
            """
            fasterq-dump {input.sra} --split-3 --skip-technical --threads {threads} \
                -O {params.outdir} --temp {params.tmp} > {log} 2>&1
            pigz -p {threads} {params.outdir}/{wildcards.run}_1.fastq {params.outdir}/{wildcards.run}_2.fastq
            rm -rf {params.tmp} {params.outdir}/{wildcards.run}.fastq
            """
