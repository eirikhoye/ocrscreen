"""Pairwise HiCRep stratum-adjusted correlation (SCC) between Hi-C libraries.

Snakemake script (rule ``hicrep_scc``). Calls the ``hicrep`` command line tool
(https://github.com/dejunlin/hicrep) for every pair of single-resolution .cool files
and writes a long table: library1, library2, chrom, scc, plus a per-pair mean.
"""
import itertools
import os
import subprocess
import tempfile

import cooler

cools = list(snakemake.input)
h = str(snakemake.params.h)
dbpmax = str(snakemake.params.dbpmax)


def lib_name(path):
    return os.path.basename(path).split(".")[0]


rows = []
with open(snakemake.log[0], "w") as log, tempfile.TemporaryDirectory() as tmp:
    for a, b in itertools.combinations(cools, 2):
        out = os.path.join(tmp, "scc.txt")
        cmd = ["hicrep", a, b, out, "--h", h, "--dBPMax", dbpmax]
        log.write(" ".join(cmd) + "\n")
        log.flush()
        subprocess.run(cmd, check=True, stdout=log, stderr=log)
        with open(out) as fh:
            sccs = [line.strip() for line in fh if line.strip() and not line.startswith("#")]
        chroms = cooler.Cooler(a).chromnames
        if len(chroms) != len(sccs):
            chroms = [f"chrom_{i + 1}" for i in range(len(sccs))]
        values = []
        for chrom, scc in zip(chroms, sccs):
            rows.append((lib_name(a), lib_name(b), chrom, scc))
            try:
                values.append(float(scc))
            except ValueError:
                pass
        finite = [v for v in values if v == v]
        mean = sum(finite) / len(finite) if finite else float("nan")
        rows.append((lib_name(a), lib_name(b), "mean", f"{mean:.6f}"))

with open(snakemake.output[0], "w") as fh:
    fh.write("library1\tlibrary2\tchrom\tscc\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")
