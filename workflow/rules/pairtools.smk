"""
Module pairtools.

Hi-C pairs -> contact maps, following the 4DN / Open2C recommendations:

1. runs (SRR) of the same library are merged *before* deduplication,
2. each library (replicate) is deduplicated separately (PCR/optical duplicates are
   library specific) and gets its own statistics and a QC .cool,
3. deduplicated pairs of all libraries of a cell line are merged (pairtools merge),
4. the merged pairs are binned (cooler cload pairs) and zoomified/balanced (.mcool).

The merged pairs also feed module pairs2hic (.hic for juicer/STARE).
"""


rule pairtools_merge_runs:
    input:
        lambda wc: expand(RUN_PAIRS, cell_line=wc.cell_line, run=hic_runs_of(wc.cell_line, wc.replicate)),
    output:
        LIB_PAIRS,
    threads: config["hic"]["threads"]["pairtools"]
    log:
        os.path.join(LOGS, "{cell_line}", "hic", "merge_runs_{replicate}.log"),
    conda:
        "../envs/hic.yaml"
    shell:
        "pairtools merge --nproc {threads} -o {output} {input} > {log} 2>&1"


rule pairtools_dedup:
    input:
        LIB_PAIRS,
    output:
        nodups=LIB_NODUPS,
        dups=LIB_DUPS,
        unmapped=LIB_UNMAPPED,
        stats=LIB_DEDUP_STATS,
    params:
        max_mismatch=config["hic"]["dedup_max_mismatch"],
    threads: config["hic"]["threads"]["pairtools"]
    log:
        os.path.join(LOGS, "{cell_line}", "hic", "dedup_{replicate}.log"),
    conda:
        "../envs/hic.yaml"
    shell:
        """
        pairtools dedup \
            --max-mismatch {params.max_mismatch} \
            --mark-dups \
            --nproc-in {threads} --nproc-out {threads} \
            --output {output.nodups} \
            --output-dups {output.dups} \
            --output-unmapped {output.unmapped} \
            --output-stats {output.stats} \
            {input} > {log} 2>&1
        """


rule cooler_library:
    """Per-library contact matrix at hic.qc_resolution (replicate QC / HiCRep)."""
    input:
        pairs=LIB_NODUPS,
        chrom_sizes=CHROM_SIZES,
    output:
        LIB_COOL,
    params:
        res=QC_RES,
        assembly=BUILD,
    log:
        os.path.join(LOGS, "{cell_line}", "hic", "cooler_{replicate}.log"),
    conda:
        "../envs/hic.yaml"
    shell:
        """
        cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 --assembly {params.assembly} \
            {input.chrom_sizes}:{params.res} {input.pairs} {output} > {log} 2>&1
        """


rule pairtools_merge_libraries:
    input:
        lambda wc: expand(LIB_NODUPS, cell_line=wc.cell_line, replicate=hic_replicates_of(wc.cell_line)),
    output:
        MERGED_PAIRS,
    threads: config["hic"]["threads"]["pairtools"]
    log:
        os.path.join(LOGS, "{cell_line}", "hic", "merge_libraries.log"),
    conda:
        "../envs/hic.yaml"
    shell:
        "pairtools merge --nproc {threads} -o {output} {input} > {log} 2>&1"


rule pairtools_stats_merged:
    """Combined statistics of all libraries of a cell line (pairtools stats --merge)."""
    input:
        lambda wc: expand(LIB_DEDUP_STATS, cell_line=wc.cell_line, replicate=hic_replicates_of(wc.cell_line)),
    output:
        MERGED_STATS,
    log:
        os.path.join(LOGS, "{cell_line}", "hic", "stats_merged.log"),
    conda:
        "../envs/hic.yaml"
    shell:
        "pairtools stats --merge {input} -o {output} > {log} 2>&1"


rule cooler_merged:
    input:
        pairs=MERGED_PAIRS,
        chrom_sizes=CHROM_SIZES,
    output:
        MERGED_COOL,
    params:
        res=BASE_RES,
        assembly=BUILD,
    log:
        os.path.join(LOGS, "{cell_line}", "hic", "cooler_merged.log"),
    conda:
        "../envs/hic.yaml"
    shell:
        """
        cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 --assembly {params.assembly} \
            {input.chrom_sizes}:{params.res} {input.pairs} {output} > {log} 2>&1
        """


rule cooler_zoomify:
    input:
        MERGED_COOL,
    output:
        MERGED_MCOOL,
    params:
        resolutions=config["hic"]["zoomify_resolutions"],
    threads: config["hic"]["threads"]["pairtools"]
    log:
        os.path.join(LOGS, "{cell_line}", "hic", "cooler_zoomify.log"),
    conda:
        "../envs/hic.yaml"
    shell:
        """
        cooler zoomify --balance --nproc {threads} -r {params.resolutions} \
            -o {output} {input} > {log} 2>&1
        """
