"""
Module hic2KR.

Dump KR-normalised intra-chromosomal contacts from a .hic file at the STARE bin size
(STARE's Juicebox_normalization.sh). The .hic is the merged map produced by
pairs2hic, or the precomputed file of a cell line when module set_HIC is enabled.
"""


rule hic2KR:
    input:
        hic=lambda wc: hic_file(wc.cell_line),
        jar=JUICER_JAR,
    output:
        contacts_dir=directory(CONTACTS_DIR),
    params:
        norm=KR["norm"],
        bin_size=KR["bin_size"],
        chromosomes=",".join(KR["chromosomes"]),
    log:
        os.path.join(LOGS, "{cell_line}", "hic", "hic2KR.log"),
    conda:
        "../envs/stare.yaml"
    shell:
        """
        bash "$(command -v Juicebox_normalization.sh)" -j {input.jar} -h {input.hic} \
            -n {params.norm} -d {output.contacts_dir} -b {params.bin_size} \
            -c {params.chromosomes} > {log} 2>&1
        """
