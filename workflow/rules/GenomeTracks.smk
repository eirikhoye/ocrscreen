"""
Module GenomeTracks.

pyGenomeTracks plot of the FDR-filtered CRISPRi enhancer activity and gene annotation
for each region in config genome_tracks.regions. Requires JACKS outputs.
"""


rule make_ini:
    input:
        bedfile=ENHANCER_BED_FILTERED,
        gtf=GTF,
    output:
        ini=TRACKS_INI,
    log:
        os.path.join(LOGS, "{cell_line}", "make_tracks_ini.log"),
    conda:
        "../envs/pygenometracks.yaml"
    shell:
        "make_tracks_file --trackFiles {input.bedfile} {input.gtf} -o {output.ini} > {log} 2>&1"


rule plot_track:
    input:
        ini=TRACKS_INI,
    output:
        TRACK_PLOT,
    params:
        region=lambda wc: config["genome_tracks"]["regions"][wc.region],
    log:
        os.path.join(LOGS, "{cell_line}", "plot_track_{region}.log"),
    conda:
        "../envs/pygenometracks.yaml"
    shell:
        """
        pyGenomeTracks --tracks {input.ini} --region {params.region} -o {output} \
            --plotWidth 600 --height 300 > {log} 2>&1
        """
