"""
Module crispr_library (was crispr_pipeline add_coords.smk).

Adds OCR coordinates to the sgRNA library of each cell line and writes a copy of the
library with a header row, which CrisprCounter and JACKS use downstream.
"""


rule merge_library_coords:
    input:
        ref=lambda wc: CELL_CFG[wc.cell_line]["library"],
        coord=lambda wc: CELL_CFG[wc.cell_line]["coords"],
    output:
        coords=LIBRARY_COORDS,
        library=LIBRARY,
    params:
        coord_cols=lambda wc: ["N20", f"Coord_{wc.cell_line}"],
    log:
        os.path.join(LOGS, "{cell_line}", "merge_library_coords.log"),
    conda:
        "../envs/crispr.yaml"
    script:
        "../scripts/merge_library_coords.py"
