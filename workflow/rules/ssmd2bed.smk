"""
Module ssmd2bed.

SSMD scores of each screen replicate as BED files at sgRNA level and OCR level
(OCR score = minimum guide Score). Requires RSL_counter and crispr_library outputs.
"""


rule ssmd2bed:
    input:
        scores=SSMD_SCORES,
        coords=LIBRARY_COORDS,
    output:
        bed=SSMD_BED,
    params:
        script=script("ssmd2bed.py"),
    log:
        os.path.join(LOGS, "{cell_line}", "ssmd2bed_{sample}_{mode}.log"),
    conda:
        "../envs/crispr.yaml"
    shell:
        """
        python {params.script} \
            -r {input.scores} \
            -c {input.coords} \
            -n {wildcards.sample} \
            -m {wildcards.mode} \
            -o {output.bed} > {log} 2>&1
        """
