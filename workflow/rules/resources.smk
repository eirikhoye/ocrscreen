"""
Download / prepare reference files and tools. Always included; rules only run when
an enabled module needs their outputs. Files that cannot be downloaded (screen
library, CrisprCounter.jar, UniBind databases) are described in resources/README.md.
"""


localrules:
    download_genome,
    download_gtf,
    download_chain,
    download_psem,
    download_dnase_peaks,
    get_juicer_tools,
    get_run_jacks,
    download_precomputed_hic,


rule download_genome:
    output:
        GENOME_FASTA,
    params:
        url=config["reference"]["fasta_url"],
    log:
        os.path.join(LOGS, "resources", "download_genome.log"),
    conda:
        "../envs/hic.yaml"
    shell:
        "(wget -q -O - {params.url} | gunzip -c > {output}) 2> {log}"


rule chrom_sizes:
    input:
        GENOME_FASTA,
    output:
        fai=GENOME_FAI,
        sizes=CHROM_SIZES,
    log:
        os.path.join(LOGS, "resources", "chrom_sizes.log"),
    conda:
        "../envs/hic.yaml"
    shell:
        """
        samtools faidx {input} 2> {log}
        cut -f1,2 {output.fai} > {output.sizes}
        """


rule download_gtf:
    output:
        GTF,
    params:
        url=config["reference"]["gtf_url"],
    log:
        os.path.join(LOGS, "resources", "download_gtf.log"),
    conda:
        "../envs/download.yaml"
    shell:
        "(wget -q -O - {params.url} | gunzip -c > {output}) 2> {log}"


rule filter_gtf_protein:
    """Gene records of protein-coding genes only (used by ABC runs with gtf: protein_coding)."""
    input:
        GTF,
    output:
        GTF_PROTEIN,
    log:
        os.path.join(LOGS, "resources", "filter_gtf_protein.log"),
    conda:
        "../envs/download.yaml"
    shell:
        """
        awk '$3 == "gene" && $0 ~ /gene_type "protein_coding"/' {input} > {output} 2> {log}
        """


rule download_chain:
    output:
        CHAIN_HG19_TO_HG38,
    params:
        url=config["reference"]["chain_hg19_to_hg38_url"],
    log:
        os.path.join(LOGS, "resources", "download_chain.log"),
    conda:
        "../envs/download.yaml"
    shell:
        "wget -q -O {output} {params.url} 2> {log}"


rule download_psem:
    """STARE PSEM (JASPAR/HOCOMOCO/Kellis human PWMs, STARE release 2.2)."""
    output:
        PSEM,
    params:
        url=config["abc_stare"]["psem_url"],
    log:
        os.path.join(LOGS, "resources", "download_psem.log"),
    conda:
        "../envs/download.yaml"
    shell:
        "wget -q -O {output} {params.url} 2> {log}"


rule download_dnase_peaks:
    """ENCODE/Duke DNase-seq peaks (hg19) from the UCSC wgEncodeOpenChromDnase track."""
    output:
        os.path.join(RESOURCES, "encode", "{name}.narrowPeak"),
    wildcard_constraints:
        name="wgEncodeOpenChromDnase[A-Za-z0-9]+Pk",
    params:
        url=lambda wc: config["reference"]["dnase_url_template"].format(name=wc.name),
    log:
        os.path.join(LOGS, "resources", "download_dnase_{name}.log"),
    conda:
        "../envs/download.yaml"
    shell:
        "(wget -q -O - {params.url} | gunzip -c > {output}) 2> {log}"


rule get_juicer_tools:
    output:
        JUICER_JAR,
    params:
        url=config["hic"]["juicer"]["url"],
    log:
        os.path.join(LOGS, "resources", "get_juicer_tools.log"),
    conda:
        "../envs/download.yaml"
    shell:
        "wget -q -O {output} {params.url} 2> {log}"


rule get_run_jacks:
    """run_JACKS.py from the JACKS repository (pinned commit; not installed by pip)."""
    output:
        RUN_JACKS,
    params:
        url=config["jacks"]["run_jacks_url"],
    log:
        os.path.join(LOGS, "resources", "get_run_jacks.log"),
    conda:
        "../envs/download.yaml"
    shell:
        "wget -q -O {output} {params.url} 2> {log}"


rule download_precomputed_hic:
    """Published .hic maps listed in config hic.precomputed_hic_urls (used with set_HIC)."""
    output:
        os.path.join(RESOURCES, "hic", "{name}.hic"),
    params:
        url=lambda wc: config["hic"]["precomputed_hic_urls"][
            os.path.join(RESOURCES, "hic", f"{wc.name}.hic")
        ],
    log:
        os.path.join(LOGS, "resources", "download_hic_{name}.log"),
    conda:
        "../envs/download.yaml"
    shell:
        "wget -q -O {output} {params.url} 2> {log}"
