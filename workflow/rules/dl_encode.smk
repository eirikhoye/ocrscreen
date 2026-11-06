"""
Module dl_encode.

Download ENCODE tracks listed in config encode (file name -> URL) into
resources/encode/tracks/ (e.g. RWPE-1 H3K27ac fold-change bigWig used for
genome-browser figures).
"""


localrules:
    download_encode_track,


rule download_encode_track:
    output:
        ENCODE_TRACK,
    params:
        url=lambda wc: config["encode"][wc.name],
    log:
        os.path.join(LOGS, "resources", "encode_{name}.log"),
    conda:
        "../envs/download.yaml"
    shell:
        "wget -c -q -O {output} {params.url} 2> {log}"
