"""
Module pairs2hic.

Merged, deduplicated pairs of a cell line -> .hic with juicer_tools pre (needed for the
KR normalisation used by STARE).
"""


rule pairs2hic:
    input:
        pairs=MERGED_PAIRS,
        chrom_sizes=CHROM_SIZES,
        jar=JUICER_JAR,
    output:
        MERGED_HIC,
    params:
        memory=config["hic"]["juicer"]["memory_mb"],
    threads: config["hic"]["juicer"]["threads"]
    resources:
        mem_mb=config["hic"]["juicer"]["memory_mb"] + 4000,
    log:
        os.path.join(LOGS, "{cell_line}", "hic", "pairs2hic.log"),
    conda:
        "../envs/juicer.yaml"
    shell:
        """
        java -Xmx{params.memory}m -Djava.awt.headless=true -jar {input.jar} pre \
            --threads {threads} {input.pairs} {output} {input.chrom_sizes} > {log} 2>&1
        """
