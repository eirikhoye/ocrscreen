"""
Shared definitions: sample sheets, output path patterns, set_HIC logic and the
collection of final targets for ``rule all`` from the enabled modules.
"""

import os
import re
from itertools import combinations

import pandas as pd


# -------------------------------------------------------------------------------------
# Modules
# -------------------------------------------------------------------------------------
KNOWN_MODULES = [
    # CRISPRi screen
    "crispr_library",
    "RSL_counter",
    "pathway",
    "ssmd2bed",
    "JACKS",
    "GenomeTracks",
    # Hi-C
    "sra",
    "bwa",
    "pairtools",
    "pairs2hic",
    "hic2KR",
    "set_HIC",
    "hic_qc",
    # Integration
    "ABC_STARE",
    "crispri_abc_significant_hits",
    "GSEA",
    "dl_encode",
    "TFenrich",
    "TFenrich_OCR_universe",
]

# Modules that share one rule file (rule names must stay unique).
MODULE_RULE_FILE = {"TFenrich_OCR_universe": "TFenrich"}

MODULES = list(config["modules"])
_unknown = [m for m in MODULES if m not in KNOWN_MODULES]
if _unknown:
    raise ValueError(f"Unknown module(s) in config['modules']: {_unknown}. Known: {KNOWN_MODULES}")


def module_rule_files(modules):
    """Rule files to include for the enabled modules (deduplicated, in order)."""
    files = []
    for m in modules:
        f = MODULE_RULE_FILE.get(m, m)
        if f not in files:
            files.append(f)
    return files


# -------------------------------------------------------------------------------------
# General paths
# -------------------------------------------------------------------------------------
RESULTS = config["results"].rstrip("/")
RESOURCES = config["resources"].rstrip("/")
LOGS = os.path.join(RESULTS, "logs")
SCRIPTS = os.path.join(workflow.basedir, "scripts")


def script(name):
    """Absolute path of a script in workflow/scripts (for shell: rules)."""
    return os.path.join(SCRIPTS, name)


# -------------------------------------------------------------------------------------
# Cell lines and sample sheets
# -------------------------------------------------------------------------------------
CELL_LINES = list(config["cell_lines"])
CELL_CFG = config["cell_lines"]

CRISPR_SHEET = pd.read_csv(config["samples"]["crispr"], sep="\t", dtype=str, comment="#")
validate(CRISPR_SHEET, "../schemas/crispr_samples.schema.yaml")
HIC_SHEET = pd.read_csv(config["samples"]["hic"], sep="\t", dtype=str, comment="#")
validate(HIC_SHEET, "../schemas/hic_samples.schema.yaml")

for _sheet, _name in [(CRISPR_SHEET, "crispr"), (HIC_SHEET, "hic")]:
    _missing = set(_sheet["cell_line"]) - set(CELL_LINES)
    if _missing:
        raise ValueError(f"{_name} sample sheet has cell lines not in config['cell_lines']: {_missing}")

STAGES = list(config["crispr"]["stages"])
if len(STAGES) != 2:
    raise ValueError("crispr.stages must list exactly two time points: [reference, late]")

# CRISPR screen replicates ("samples" in the original pipeline), e.g. RWPE_R1
CRISPR_SAMPLES = sorted(CRISPR_SHEET["sample"].unique())
CRISPR_CELL_LINES = [c for c in CELL_LINES if c in set(CRISPR_SHEET["cell_line"])]


def crispr_samples_of(cell_line):
    return sorted(CRISPR_SHEET.loc[CRISPR_SHEET["cell_line"] == cell_line, "sample"].unique())


def crispr_rows(sample):
    """Sample-sheet rows (as dicts) of one screen replicate."""
    rows = CRISPR_SHEET[CRISPR_SHEET["sample"] == sample].to_dict("records")
    if not rows:
        raise ValueError(f"CRISPR sample {sample} not in sample sheet")
    return rows


def crispr_rows_of(cell_line):
    return CRISPR_SHEET[CRISPR_SHEET["cell_line"] == cell_line].to_dict("records")


def crispr_cell_of(sample):
    return crispr_rows(sample)[0]["cell_line"]


def crispr_replicate_pair(cell_line):
    """The two screen replicates of a cell line (used for R1/R2 comparisons)."""
    samples = crispr_samples_of(cell_line)
    if len(samples) != 2:
        raise ValueError(
            f"{cell_line}: replicate-comparison steps expect exactly two screen replicates, got {samples}"
        )
    return samples


def replicate_pair_files(pattern):
    """Input function: ``pattern`` for the two screen replicates of {cell_line}."""

    def _f(wc):
        return [
            pattern.format(cell_line=wc.cell_line, sample=s)
            for s in crispr_replicate_pair(wc.cell_line)
        ]

    return _f


# Hi-C: run (SRR) -> replicate/library -> cell line
HIC_CELL_LINES_ALL = [c for c in CELL_LINES if c in set(HIC_SHEET["cell_line"])]


def hic_runs_of(cell_line, replicate=None):
    sel = HIC_SHEET["cell_line"] == cell_line
    if replicate is not None:
        sel &= HIC_SHEET["replicate"] == replicate
    return sorted(HIC_SHEET.loc[sel, "run"].unique())


def hic_replicates_of(cell_line):
    return sorted(HIC_SHEET.loc[HIC_SHEET["cell_line"] == cell_line, "replicate"].unique())


# -------------------------------------------------------------------------------------
# set_HIC: use a precomputed .hic file or contact directory instead of processing reads
# -------------------------------------------------------------------------------------
USE_SET_HIC = "set_HIC" in MODULES


def uses_precomputed_hic(cell_line):
    if not USE_SET_HIC:
        return False
    c = CELL_CFG[cell_line]
    return bool(c.get("precomputed_hic") or c.get("precomputed_contacts"))


if USE_SET_HIC and not any(uses_precomputed_hic(c) for c in CELL_LINES):
    raise ValueError(
        "Module set_HIC is enabled but no cell line sets precomputed_hic or "
        "precomputed_contacts. Set one of them or remove set_HIC from modules."
    )

# Cell lines whose Hi-C is processed from raw reads by this workflow
HIC_CELL_LINES = [c for c in HIC_CELL_LINES_ALL if not uses_precomputed_hic(c)]


# -------------------------------------------------------------------------------------
# Wildcard constraints
# -------------------------------------------------------------------------------------
def _alternatives(values):
    return "|".join(re.escape(v) for v in values) if values else "(?!x)x"


wildcard_constraints:
    cell_line=_alternatives(CELL_LINES),
    sample=_alternatives(CRISPR_SAMPLES),
    run="[^/.]+",
    replicate="[^/.]+",
    abc_run=_alternatives(list(config["abc_stare"]["runs"])),
    region="[^/.]+",
    gene_set="[^/]+",
    mode="sgrna|ocr",


# -------------------------------------------------------------------------------------
# References and tools (rules in resources.smk)
# -------------------------------------------------------------------------------------
BUILD = config["reference"]["build"]
GENOME_DIR = os.path.join(RESOURCES, "genome")
GENOME_FASTA = os.path.join(GENOME_DIR, f"{BUILD}.fa")
GENOME_FAI = f"{GENOME_FASTA}.fai"
CHROM_SIZES = os.path.join(GENOME_DIR, f"{BUILD}.chrom.sizes")
BWA_INDEX_PREFIX = os.path.join(GENOME_DIR, "bwa-mem2", BUILD)
BWA_INDEX_EXT = ["0123", "amb", "ann", "bwt.2bit.64", "pac"]
BWA_INDEX = [f"{BWA_INDEX_PREFIX}.{ext}" for ext in BWA_INDEX_EXT]
GTF = os.path.join(GENOME_DIR, os.path.basename(config["reference"]["gtf_url"]).removesuffix(".gz"))
GTF_PROTEIN = f"{GTF.removesuffix('.gtf')}.protein_coding.genes.gtf"
CHAIN_HG19_TO_HG38 = os.path.join(
    GENOME_DIR, os.path.basename(config["reference"]["chain_hg19_to_hg38_url"])
)
PSEM = os.path.join(RESOURCES, "stare", os.path.basename(config["abc_stare"]["psem_url"]))
JUICER_JAR = config["hic"]["juicer"]["jar"]
RUN_JACKS = os.path.join(RESOURCES, "tools", "run_JACKS.py")


# -------------------------------------------------------------------------------------
# Output patterns: CRISPRi screen  (results/<CELL>/crispr/...)
# -------------------------------------------------------------------------------------
CR = os.path.join(RESULTS, "{cell_line}", "crispr")
LIBRARY = os.path.join(CR, "library", "{cell_line}_library.csv")
LIBRARY_COORDS = os.path.join(CR, "library", "{cell_line}_library_coords.csv")
SAMPLE_TABLE = os.path.join(CR, "sample_files", "{sample}.tsv")
COUNTER_PROPERTIES = os.path.join(CR, "config_files", "{sample}.properties")
RSL_COUNTS = os.path.join(CR, "RSL_counts", "{sample}.csv")
PROCESSED_COUNTS = os.path.join(CR, "RSL_counts_processed", "{sample}_processed.csv")
RSL_COUNTS_BINNED = os.path.join(CR, "RSL_binned_counts", "{sample}.csv")
SSMD_PREFIX = os.path.join(CR, "SSMD_score", "{sample}_RSL")
SSMD_SCORES = f"{SSMD_PREFIX}_SSMD_Score.tsv"
SSMD_NORMALIZED = f"{SSMD_PREFIX}_normalized.csv"
SSMD_BED = os.path.join(CR, "SSMD_score", "{sample}_SSMD_score_{mode}.bed")
GMT_FILE = os.path.join(CR, "RSL_counts_GMT", "{sample}_processed.gmt")
PATHWAY_PREFIX = os.path.join(CR, "mageck_pathway", "{sample}")
PATHWAY_OUT = f"{PATHWAY_PREFIX}.pathway_summary.txt"
SSMD_FIGURE = os.path.join(CR, "SSMD_figures", "{cell_line}_RSL_SSMD_score.pdf")

# JACKS (results/<CELL>/jacks/...)
JK = os.path.join(RESULTS, "{cell_line}", "jacks")
AGGREGATE_RSL = os.path.join(JK, "aggregated_counts.csv")
RSL_SYNTHETIC_REPLICATES = os.path.join(JK, "RSL_synthetic_replicates.csv")
REPMAP = os.path.join(JK, "repmap.tab")
JACKS_PREFIX = os.path.join(JK, "{cell_line}_result")
JACKS_SUFFIXES = [
    "_gene_JACKS_results.txt",
    "_gene_pval_JACKS_results.txt",
    "_gene_std_JACKS_results.txt",
    "_grna_JACKS_results.txt",
    "_JACKS_results_full.pickle",
    "_logfoldchange_means.txt",
    "_logfoldchange_std.txt",
]
JACKS_OUTPUTS = [f"{JACKS_PREFIX}{s}" for s in JACKS_SUFFIXES]
JACKS_GENE = f"{JACKS_PREFIX}_gene_JACKS_results.txt"
JACKS_PVAL = f"{JACKS_PREFIX}_gene_pval_JACKS_results.txt"
ENHANCER_BED_UNFILTERED = os.path.join(JK, "enhancer_activity_per_cell.bed")
ENHANCER_BED_FILTERED = os.path.join(JK, "enhancer_activity_per_cell_filtered.bed")
ENHANCER_NARROWPEAK = os.path.join(JK, "enhancer_activity_smallpeaks.narrowPeak")

# Genome tracks
TRACKS_INI = os.path.join(RESULTS, "{cell_line}", "genome_tracks", "tracks.ini")
TRACK_PLOT = os.path.join(RESULTS, "{cell_line}", "genome_tracks", "{region}.png")

# -------------------------------------------------------------------------------------
# Output patterns: Hi-C  (results/<CELL>/hic/...)
# -------------------------------------------------------------------------------------
HC = os.path.join(RESULTS, "{cell_line}", "hic")
SRA_FILE = os.path.join(HC, "fastq", "{run}.sra")
FASTQ_1 = os.path.join(HC, "fastq", "{run}_1.fastq.gz")
FASTQ_2 = os.path.join(HC, "fastq", "{run}_2.fastq.gz")
RUN_PAIRS = os.path.join(HC, "runs", "{run}.sorted.pairs.gz")
RUN_PARSE_STATS = os.path.join(HC, "runs", "{run}.parse.stats")
LIB_PAIRS = os.path.join(HC, "libraries", "{replicate}.sorted.pairs.gz")
LIB_NODUPS = os.path.join(HC, "libraries", "{replicate}.nodups.pairs.gz")
LIB_DUPS = os.path.join(HC, "libraries", "{replicate}.dups.pairs.gz")
LIB_UNMAPPED = os.path.join(HC, "libraries", "{replicate}.unmapped.pairs.gz")
LIB_DEDUP_STATS = os.path.join(HC, "libraries", "{replicate}.dedup.stats")
QC_RES = config["hic"]["qc_resolution"]
LIB_COOL = os.path.join(HC, "libraries", f"{{replicate}}.{QC_RES}.cool")
MERGED_PAIRS = os.path.join(HC, "{cell_line}.merged.nodups.pairs.gz")
MERGED_STATS = os.path.join(HC, "{cell_line}.merged.stats")
BASE_RES = config["hic"]["cooler_base_resolution"]
MERGED_COOL = os.path.join(HC, f"{{cell_line}}.{BASE_RES}.cool")
MERGED_MCOOL = os.path.join(HC, "{cell_line}.mcool")
MERGED_HIC = os.path.join(HC, "{cell_line}.hic")
KR = config["hic"]["kr"]
CONTACTS_DIR = os.path.join(HC, f"contacts_{KR['norm']}_{KR['bin_size']}")
HICREP = os.path.join(HC, "qc", "{cell_line}.hicrep_scc.tsv")


def hic_file(cell_line):
    """The .hic file used for one cell line (precomputed with set_HIC, else merged)."""
    pre = CELL_CFG[cell_line].get("precomputed_hic")
    if uses_precomputed_hic(cell_line) and pre:
        return pre
    return MERGED_HIC.format(cell_line=cell_line)


def contacts_dir(cell_line):
    """KR-normalised contact directory used by STARE for one cell line."""
    pre = CELL_CFG[cell_line].get("precomputed_contacts")
    if uses_precomputed_hic(cell_line) and pre:
        return pre
    return CONTACTS_DIR.format(cell_line=cell_line)


# -------------------------------------------------------------------------------------
# Output patterns: integration
# -------------------------------------------------------------------------------------
ABC_DIR = os.path.join(RESULTS, "{cell_line}", "abc_stare", "{abc_run}")
ABC_METADATA = os.path.join(ABC_DIR, "{abc_run}_metadata.amd.tsv")
ABC_SCORES = os.path.join(ABC_DIR, "{cell_line}_{abc_run}_ABC_scoredInteractions.txt.gz")
ABC_BEDPE = os.path.join(ABC_DIR, "{cell_line}_{abc_run}_abc_scores.bedpe")

HITS_DIR = os.path.join(RESULTS, "{cell_line}", "crispri_abc")
SIG_CONTACTS = os.path.join(HITS_DIR, "{cell_line}_sig_contacts.tsv")
SIG_CRE_BED = os.path.join(HITS_DIR, "{cell_line}_sig_CRE_hg19.bed")

GSEA_DIR = os.path.join(RESULTS, "{cell_line}", "gsea")
ABC_FILTERED_FDR = os.path.join(GSEA_DIR, "{cell_line}_abc_scores_filtered_jacks_fdr.txt")
ABC_RANK = os.path.join(GSEA_DIR, "{cell_line}_abc_scores.rnk")
GSEA_RES = os.path.join(GSEA_DIR, "{gene_set}", "gsea_results.tsv")

TF_DIR = os.path.join(RESULTS, "{cell_line}", "tfenrich")
SIG_CRE_HG38 = os.path.join(TF_DIR, "{cell_line}_sig_CRE_hg38.bed")
SIG_CRE_UNMAPPED = os.path.join(TF_DIR, "{cell_line}_sig_CRE_unmapped.bed")
OCR_UNIVERSE_HG19 = os.path.join(TF_DIR, "{cell_line}_ocr_universe_hg19.bed")
OCR_UNIVERSE_HG38 = os.path.join(TF_DIR, "{cell_line}_ocr_universe_hg38.bed")
OCR_UNIVERSE_UNMAPPED = os.path.join(TF_DIR, "{cell_line}_ocr_universe_unmapped.bed")
UNIBIND_OUT = os.path.join(TF_DIR, "unibind_universe", "allEnrichments.tsv")
UNIBIND_OUT_OCR = os.path.join(TF_DIR, "ocr_universe", "allEnrichments.tsv")

ENCODE_TRACK = os.path.join(RESOURCES, "encode", "tracks", "{name}")


# -------------------------------------------------------------------------------------
# Input functions
# -------------------------------------------------------------------------------------
def counter_fastqs(wc):
    """Guide (R1) and RSL (R2) FASTQs of one screen replicate."""
    fq_dir = CELL_CFG[wc.cell_line]["fastq_dir"]
    fastqs = []
    for row in crispr_rows(wc.sample):
        r1 = row["fastq"]
        fastqs += [os.path.join(fq_dir, r1), os.path.join(fq_dir, r1.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz"))]
    return fastqs


def binned_counts_of_cell(wc):
    """Binned RSL counts of all screen replicates of {cell_line}."""
    return [
        RSL_COUNTS_BINNED.format(cell_line=wc.cell_line, sample=s)
        for s in crispr_samples_of(wc.cell_line)
    ]


ABC_RUNS = config["abc_stare"]["runs"]


def abc_activity_bed(wc):
    """Enhancer-activity BED of an ABC run (config abc_stare.runs.<run>.activity)."""
    activity = ABC_RUNS[wc.abc_run]["activity"]
    if activity == "dnase":
        return CELL_CFG[wc.cell_line]["dnase_peaks"]
    if activity == "jacks_fdr":
        return ENHANCER_BED_FILTERED.format(cell_line=wc.cell_line)
    if activity == "jacks_all":
        return ENHANCER_BED_UNFILTERED.format(cell_line=wc.cell_line)
    raise ValueError(f"Unknown abc_stare activity '{activity}' (run {wc.abc_run})")


def abc_gtf(wc):
    """Gene annotation of an ABC run (all genes or protein-coding genes only)."""
    return GTF_PROTEIN if ABC_RUNS[wc.abc_run]["gtf"] == "protein_coding" else GTF


# -------------------------------------------------------------------------------------
# Final targets for rule all
# -------------------------------------------------------------------------------------
def get_final_output(modules):
    out = []
    xs = CRISPR_SAMPLES
    xc = CRISPR_CELL_LINES

    def per_sample(pattern, **kw):
        return [pattern.format(cell_line=crispr_cell_of(s), sample=s, **kw) for s in xs]

    if "crispr_library" in modules:
        out += expand(LIBRARY_COORDS, cell_line=xc)

    if "RSL_counter" in modules:
        for p in [SAMPLE_TABLE, COUNTER_PROPERTIES, RSL_COUNTS, PROCESSED_COUNTS, RSL_COUNTS_BINNED, SSMD_SCORES]:
            out += per_sample(p)

    if "pathway" in modules:
        out += per_sample(GMT_FILE) + per_sample(PATHWAY_OUT)
        out += expand(SSMD_FIGURE, cell_line=xc)

    if "ssmd2bed" in modules:
        for mode in ["sgrna", "ocr"]:
            out += per_sample(SSMD_BED, mode=mode)

    if "JACKS" in modules:
        out += expand(
            [AGGREGATE_RSL, REPMAP, RSL_SYNTHETIC_REPLICATES, ENHANCER_BED_FILTERED,
             ENHANCER_BED_UNFILTERED, ENHANCER_NARROWPEAK],
            cell_line=xc,
        )
        out += expand(JACKS_OUTPUTS, cell_line=xc)

    if "GenomeTracks" in modules:
        out += expand(TRACK_PLOT, cell_line=xc, region=list(config["genome_tracks"]["regions"]))

    # Hi-C (only cell lines processed from reads; set_HIC cell lines are skipped)
    for c in HIC_CELL_LINES:
        reps = hic_replicates_of(c)
        runs = hic_runs_of(c)
        if "sra" in modules:
            out += expand([FASTQ_1, FASTQ_2], cell_line=c, run=runs)
        if "bwa" in modules:
            out += expand([RUN_PAIRS, RUN_PARSE_STATS], cell_line=c, run=runs)
        if "pairtools" in modules:
            out += expand([LIB_NODUPS, LIB_DEDUP_STATS, LIB_COOL], cell_line=c, replicate=reps)
            out += expand([MERGED_PAIRS, MERGED_STATS, MERGED_MCOOL], cell_line=c)
        if "pairs2hic" in modules:
            out.append(MERGED_HIC.format(cell_line=c))
        if "hic_qc" in modules and len(reps) > 1:
            out.append(HICREP.format(cell_line=c))

    if "hic2KR" in modules:
        # contacts computed here: from the merged .hic or from a precomputed .hic
        for c in CELL_LINES:
            if (c in HIC_CELL_LINES or uses_precomputed_hic(c)) and contacts_dir(c) == CONTACTS_DIR.format(cell_line=c):
                out.append(contacts_dir(c))

    if "ABC_STARE" in modules:
        for c in CELL_LINES:
            for r in config["abc_stare"]["runs"]:
                out += [p.format(cell_line=c, abc_run=r) for p in (ABC_METADATA, ABC_SCORES, ABC_BEDPE)]

    if "crispri_abc_significant_hits" in modules:
        out += expand([SIG_CONTACTS, SIG_CRE_BED], cell_line=xc)

    if "GSEA" in modules:
        out += expand(GSEA_RES, cell_line=xc, gene_set=config["gsea"]["gene_sets"])

    if "dl_encode" in modules:
        out += expand(ENCODE_TRACK, name=list(config.get("encode") or {}))

    if "TFenrich" in modules:
        out += expand(UNIBIND_OUT, cell_line=xc)

    if "TFenrich_OCR_universe" in modules:
        out += expand(UNIBIND_OUT_OCR, cell_line=xc)

    return list(dict.fromkeys(out))
