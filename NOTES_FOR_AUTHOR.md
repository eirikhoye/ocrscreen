# Notes for the author

This file documents how `ocrscreen` was assembled from `crispr_pipeline` and `hic_abc`
(both left untouched). It is meant for you, not for readers of the paper repository —
**delete it before publishing**.

Built unattended on 2026-11-06 from crispr_pipeline HEAD `3956e19` (2026-07-10; the
working tree had one uncommitted notebook change, ignored) and hic_abc HEAD `78bd8e2`
(2026-09-29).

---

## 1. Validation results

Snakemake 9.14.3 (`crispr_pipeline/.venv/bin/snakemake`, run with cwd = `ocrscreen`,
`PYTHONDONTWRITEBYTECODE=1`).

| check | result |
|---|---|
| `snakemake --lint --configfile config/config.yaml` | **"Congratulations, your workflow is in a good condition!"** |
| `snakemake -n --configfile config/config.yaml` | parses; stops with `MissingInputException` for `resources/crispr/LibraryFile_RWPE.csv` / `Final_List_RWPE.csv` (private inputs, expected — see §4) |
| dry-run, all 18 modules except `set_HIC`, dummy inputs (`/storage/scratch/tmp/ocrscreen_test/test_config.yaml`) | DAG builds, **191 jobs**: 11 SRA runs → 9 libraries (RWPE rep1/rep2; HMEC HIC058–HIC064, HIC064 = 3 runs merged before dedup) → 2 merged maps; 8 STARE runs (4 per cell line) |
| same + `set_HIC`, HMEC from a precomputed `.hic` | 137 jobs (only RWPE reads processed; HMEC `hic2KR` uses the given `.hic`) |
| same + RWPE from a precomputed contact dir | 116 jobs (no `hic2KR` for RWPE) |
| `set_HIC` with `resources/hic/GSE63525_HMEC_combined.hic` | 138 jobs incl. GEO download |
| `--config 'modules=[crispr_library,RSL_counter,pathway]'` | 37 jobs |
| Hi-C only + `hic.max_spots` (subsampling) | 68 jobs, uses `fastq_dump_subsample` |
| `snakemake --lint` with the test override | only lint: UniBind rules have no conda/container, because the test sets `unibind.container: null` |
| `python -m py_compile workflow/scripts/*.py` | OK (18 files) |
| `perl -c bin_truncated_rsl_counts.pl` | OK |
| R script | not checked (no `Rscript` on the host) |

Functional spot check (no pipeline jobs were run): `crispri_ocr_abc_hits.py` (modified)
was run on the existing `results_rwpe_RSLnt3` files + `ABC_SRR7760386` ABC output
(read-only, output in scratch). The written `--sig-bed` is **identical** to
`hic_abc/workflow/notebooks/rwpe_sig_CRE_hg19.bed` (445 OCRs). The contacts table has the
same 3033 rows and values as `hic_abc/test_sig_contacts.tsv` except `FDR_R1/FDR_R2/
concordance_tier`, which is explained by the pathway summaries having been regenerated
with 10⁶ permutations (2026-06-22) after that test file was written (2026-06-03); p-values
are identical.

**Not validated**: conda environments were not built (no conda/mamba on the host) —
see §3 for pins that are unverified.

## 2. Design decisions

### Structure
- Standard layout (`workflow/Snakefile`, `workflow/rules/*.smk`, `envs/`, `scripts/`,
  `schemas/`, `config/`). One rule file per module; `common.smk` holds sample-sheet
  parsing, path patterns, input functions and `get_final_output()`; `resources.smk`
  (always included) holds download rules for references/tools.
- One `modules:` list (merges `enabled_modules` and `modules`). As in both originals, a
  module contributes its rules *and* its targets; inputs from a disabled module must
  already exist. `TFenrich` and `TFenrich_OCR_universe` share one rule file
  (`TFenrich.smk`) because they share the liftOver step (avoids the original rule-name
  collision `bed_liftover`/`unibind_enrich`).
- Everything is namespaced as `results/<CELL>/...` and wildcards `{cell_line}`,
  `{sample}` (screen replicate, e.g. `RWPE_R1`), `{run}`, `{replicate}`, `{abc_run}`.
  RWPE and HMEC run in one invocation; global variables such as `CELL = CELLS[0]`,
  `SAMPLES`, `HIC`, `UNMAPPED` (defined twice in hic_abc) are gone.
- CRISPR samples moved from nested YAML to `config/crispr_samples.tsv`; Hi-C runs to
  `config/hic_samples.tsv` (cell_line, replicate, run). Both validated with JSON schemas.
- Inline `run:` blocks were moved to `workflow/scripts/*.py` (`script:` directive) with
  unchanged logic: merge_library_coords, make_crisprcounter_sample_table, aggregate_rsl,
  rsl_synthetic_replicates (incl. `pivot_counts_by_umi`), abc_to_bedpe,
  filter_abc_by_jacks_fdr, rank_genes_for_gsea.
- No absolute paths in rules/scripts. `module load …` lines were replaced by conda envs.

### Parameters taken as "the paper run"
- CRISPR: the latest per-cell-line configs `config_rwpe_RSLnt3.yaml` /
  `config_hmec_RSLnt3.yaml` (used by hic_abc notebooks via `results_*_RSLnt3`):
  `nr_nt_RSL: 3`, `bin_min_count: 0`, `SSMD_trunclen: 1`, `SSMD_mincount: 0`,
  `mageck_pathway_permutations: 1000000`, `fdr_threshold: 0.05`, neg. controls NCR/NCH,
  positive controls `ctrl_genes_single.txt`. `config_rwpe_hmec.yaml` (combined library
  `HMEC_RWPE_result_NEW`, old `merge_refs`) is an older design and was not ported.
- Hi-C/ABC: `hic_abc/config/config_{rwpe,hmec}.yaml` (bin 5000, KR, `all_tss`, DNase
  narrowPeak column 7, juicer 48 GB).
- The SSMD keys were misnamed in the original (`SSMD_trunclen` is the R script's
  *minimum control count*, `SSMD_mincount` its *filter-on-sum flag*); renamed to
  `crispr.ssmd.min_control_count` / `filter_on_sum` with identical values.

### ABC/STARE merge
- One `ABC_STARE` module; each entry of `abc_stare.runs` = activity source × STARE options,
  always using the **cell line's own Hi-C contacts** from `hic2KR` (or `set_HIC`):
  - `dnase` = hic_abc design (DNase narrowPeak, `-n 7`, `all_tss`) — used by
    `crispri_abc_significant_hits` (as in the notebooks);
  - `jacks_fdr`, `jacks_all`, `jacks_fdr_protein` = crispr_pipeline designs
    (`-n 5+`, `5_tss`; filtered/unfiltered JACKS activity; protein-coding GTF) — the last
    one feeds `GSEA` as in `RUN_GSEA.smk`.
- STARE output names depend on the column name (`…_scoredInteractions_c7.txt.gz`,
  `…_RWPE.txt.gz`); the rule copies the file to a stable
  `<CELL>_<run>_ABC_scoredInteractions.txt.gz`.
- `stare_rmtree.py` is kept (identical in both repos).

### Hi-C best practice (new)
- Per run: `bwa-mem2 mem -SP5M | pairtools parse --min-mapq 30 --walks-policy 5unique
  --max-inter-align-gap 30 --add-columns mapq --drop-sam --drop-seq | pairtools sort`.
- Runs of one library merged (`pairtools merge`) **before** `pairtools dedup
  --max-mismatch 3 --mark-dups` per library; per-library `dedup.stats` (pairtools stats
  format), per-library QC `.cool` at 100 kb; libraries merged with `pairtools merge`;
  merged `.cool` at 1 kb → `cooler zoomify --balance` (`.mcool`); merged pairs →
  `juicer_tools pre` → `.hic` → `Juicebox_normalization.sh` KR contacts at 5 kb.
- Combined stats with `pairtools stats --merge`. Optional `hic_qc` (HiCRep SCC between
  libraries) — not enabled in `config.yaml`.
- Decision: MAPQ ≥ 30 and `5unique` follow common 4DN/Open2C practice; hic_abc used
  pairtools defaults. Change `hic.parse_options` if you prefer the old behaviour.
- No BAMs are written (bwa output is streamed into pairtools).
- HMEC: each GEO sample (HIC058–HIC064, all "in situ Hi-C" per GEO; checked for HIC058)
  is treated as one library/replicate; HIC064 has 3 runs (SRR1658686–88) of one library.
  All 7 libraries are merged, like Rao et al.'s "combined" map.

### Accessions found / looked up
- RWPE-1: SRR7760386 / SRR7760387 (GSM3358193 rep1, GSM3358194 rep2, GSE118514,
  PRJNA485892) — from hic_abc configs, confirmed via ENA.
- HMEC: not in any config (the original used `HMEC_combined.hic`, 7.6 GB, matching GEO
  `GSE63525_HMEC_combined.hic`). Runs looked up via NCBI E-utilities: GSM1551607–
  GSM1551613 → SRR1658680–SRR1658688.
- DNase peaks: the files used (`/storage/mathelierarea/processed/rozabl/.../wgEncodeOpenChromDnase{Rwpe1,Hmec}Pk.narrowPeak`)
  are md5-identical to the gunzipped UCSC downloads; the STARE PSEM is md5-identical to
  the STARE GitHub file; `run_JACKS.py` is identical to JACKS commit dd5c4be.

## 3. Behaviour changes vs. the original code (please review)

1. **RWPE ABC used one replicate only, without deduplication**: per hic_abc README the
   `.hic` was built from `SRR7760386.sorted.pairs.gz` (sorted, *not* deduped), and
   notebooks read `ABC_SRR7760386`. Now: both replicates, deduplicated per library, merged.
   RWPE ABC scores and everything downstream (sig contacts, Fig. 3) will change.
2. **hic_abc parsed a coordinate-sorted BAM** (`samtools sort` after bwa, then
   `samtools view | pairtools parse`). pairtools parse needs read-name–grouped input, so
   the original pairs were likely mis-parsed (mates not adjacent). Fixed by streaming.
3. bwa flags: original `-5SP`; now `-SP5M` (4DN).
4. **HMEC now reprocessed from reads by default.** To reproduce the original exactly,
   enable `set_HIC` and set `cell_lines.HMEC.precomputed_hic:
   resources/hic/GSE63525_HMEC_combined.hic` (downloaded automatically).
5. **crispr_pipeline ABC_STARE used the average hg19 Hi-C** (`average.hic.v2.191020`,
   ABC paper; rule `copy_KRobserved`) with `-k 5`. STARE's `-k` is the Hi-C bin size in
   bp, so `-k 5` looks like a unit bug. All runs now use cell-line Hi-C with
   `-k 5000`. To use the average Hi-C again, point `precomputed_contacts` at a contact
   directory with `set_HIC` (and check the file layout STARE expects; `-k` comes from
   `hic.kr.bin_size`). `copy_KRobserved` was dropped.
6. `plot_SSMD` (RSL_counter.smk) had a stray quote (`--output {output} '`) and
   hard-coded `"NCR"` for negative controls → fixed (HMEC uses NCH).
7. `RSL_counter.smk` and `pathway.smk` both defined `make_gmt`/`mageck_pathway`; the
   RSL_counter version called `convert_RSL_counts_to_gmt.py` with the old 2-argument
   signature (would fail). Only the `pathway` version is kept; `plot_SSMD` moved to
   `pathway` (it needs the pathway output).
8. hic2KR/ABC_STARE in hic_abc used non-wildcard outputs (`contacts`, `ABC_output`) that
   only worked for one sample; now per cell line.
9. JACKS, ABC `jacks_*` runs and GSEA were only ever run for RWPE (HMEC result dirs have
   no JACKS outputs); they now run for both cell lines.
10. SRA: `--split-files` → `--split-3 --skip-technical` + pigz; `.sra` is `temp()`.
11. GSEA: six near-identical rules → one rule with a `{gene_set}` wildcard; outputs move
    from `gsea*/gsea_results.tsv` to `gsea/<GENE_SET>/gsea_results.tsv`.
12. `crispri_abc_significant_hits.smk` was empty in hic_abc; implemented with
    `crispri_ocr_abc_hits.py`, which now also writes `<CELL>_sig_CRE_hg19.bed`
    (previously written by hand in the `analyze_abc_score_(old_for_reference)*` notebooks)
    and takes `--coord-col` (was hard-coded `Coord_RWPE`, so HMEC would have failed).
13. TFenrich used the binary `./resources/liftOver` → conda `ucsc-liftover`.
14. Results layout and some file names changed (`results/<CELL>/crispr/...` instead of
    `results_rwpe_RSLnt3/...`; sample tables are `.tsv`). Notebooks will need new paths.
15. `make_enhancer_activity_bed.py` + `_unfiltered.py` merged (`--no_fdr_filter`).
16. GenomeTracks: the "TEST REMOVE" rules (`make_config`, `plot_bigwig_track`, using
    pyGenomeTracks example files) and the H3K27ac bigBed conversion (`ConvertFiles`) were
    dropped; region moved to config. The odd plot size (`--plotWidth 600 --heigh 300`)
    was kept as `--plotWidth 600 --height 300`.

## 4. Inputs you must provide (original locations, for your local config)

Use a local override config (not committed) such as `config/local.yaml` and run with
`--configfile config/config.yaml config/local.yaml`, or symlink into `resources/`/`data/`:

| key | original path |
|---|---|
| `cell_lines.RWPE.fastq_dir` | `/storage/mathelierarea/processed/eirikhoy/cisRegVar/data/cisCRISPRi_data_from_KI/screening_data/rwpe1/Anthony_RWPE/Anthony_RWPE_FASTQ/Read1_guide/` |
| `cell_lines.HMEC.fastq_dir` | `/storage/mathelierarea/processed/eirikhoy/cisRegVar/data/cisCRISPRi_data_from_KI/screening_data/hmec/Anthony_FASTQ_Me16C/` |
| `cell_lines.RWPE.library` | `…/screening_data/rwpe1/Anthony_RWPE/Anthony_RWPE_Readcounting/LibraryFile_RWPE.csv` |
| `cell_lines.HMEC.library` | `…/screening_data/hmec/Anthony_Me16C_Results/LibraryFile_Me16C.csv` |
| `cell_lines.<CELL>.coords` | `/storage/mathelierarea/processed/eirikhoy/cisRegVar/data/metadata/Final_List_{RWPE,HMEC}.csv` |
| `crispr.crisprcounter.jar` | `/storage/mathelierarea/processed/eirikhoy/CRISPR-pooled-RSL/containers/crisprcounter-perl/CrisprCounter.jar` |
| `unibind.enrich_script` | `/storage/mathelierarea/processed/anthoma/Projects/unibind_enrichment/bin/UniBind_enrich.sh` |
| `unibind.container` | `hic_abc/resources/unibind_enrichment_0.2.sif` |
| `unibind.lola_db`, `lola_universe` | `hic_abc/resources/hg38_robust_UniBind_LOLA{,_universe}.RDS` |
| juicer jar (optional, else downloaded) | `hic_abc/resources/juicer_tools_1.22.01.jar` |

The library and coordinate tables are 3–12 MB each (above the 1 MB limit for
`resources/`) and are study data → release them with the paper (supplement/Zenodo) and
fill in `resources/README.md`.

UniBind with Singularity: as in the hic_abc README you probably need
`--singularity-args "--bind <unibind_enrichment dir>,<repo>,<tmp>"`.

## 5. Unverified assumptions / open questions

- **Conda pins not solved here**: `sra-tools=3.1.1`, `pigz=2.8`, `r-reshape2=1.4.4`,
  `r-plyr=1.8.9`, `r-base=4.3.2`, `wget=1.21.4`, `perl=5.32`, `openjdk=17`, `hicrep==0.2.6`
  were chosen without a solver. Pins for crispr/gseapy/pygenometracks/jacks/stare/hic come
  from built envs (`crispr_pipeline/.snakemake/conda/*`, `hic_abc/.venv`). Run
  `snakemake --use-conda --conda-create-envs-only` once and adjust.
- CrisprCounter originally ran with the system/module Java; now `openjdk=17` from conda.
- `juicer_tools pre` on pairtools ≥1.0 pairs with extra `mapq1/mapq2` columns: the
  original ran it on pairs without extra columns. I expect juicer to ignore extra
  columns, but please check on a subsample (`hic.max_spots`).
- `Juicebox_normalization.sh` is called from the `stare-abc` conda package (present in
  `bin/` of the built env); it needs Java → `openjdk` added to `stare.yaml`.
- `hicrep` CLI call/format in `hicrep_pairwise.py` (`hicrep a.cool b.cool out --h --dBPMax`)
  follows the hicrep README; not run.
- UniBind: whether `UniBind_enrich.sh` must live inside the container or be bind-mounted
  is as in the original (bind-mounted).
- GSEApy with Enrichr library names needs internet at run time (as in the original).
- `plot_SSMD.py`/`ssmd2bed.py`/`make_enhancer_activity_bed.py` assume the late time point
  is `D14` and two replicates named `<CELL>_R1`, `<CELL>_R2`.
- The HMEC library is called `Me16C` (file name); I assumed this is the HMEC (HMEC/Me16C)
  library as configured in `config_hmec*.yaml`.
- Study summary in README is inferred from code/configs (D4 vs D14 CRISPRi fitness screen
  tiling OCRs); please rewrite with the paper abstract.
- Paper figure mapping in README is partial (Fig. 2/suppl. unknown; Fig. 3 notebook also
  contains its own GSEA on ABC contacts, which is *not* the `GSEA` module).

## 6. Files deliberately excluded

**crispr_pipeline** (tracked files not copied):
- Results/logs/figures/binaries: `results/`, `results_rwpe/`, `results_rwpe_hmec/`, `logs/`,
  `figures/*.svg`, `abc_pairplot.png`, `hmec_DAG.png`, `output.png`, `nohup.out`, `temp.csv`,
  `coord_search.txt`, `.complete`, `scripts` (symlink to juicer CPU scripts).
- Notebooks and their outputs: `workflow/notebooks/**` (all `.ipynb`, `.png/.svg/.pdf`,
  `rwpe_df`, `your_file.bg.bgz*`).
- Resources: `JASPAR2024_CORE_non-redundant_pfms_transfac.txt` (unused),
  `ctrl_genes.txt`, `ctrl_genes_pipe.txt`, `ctrl_genes_single_no_comma.txt` (only used by
  the old combined config / `ctrl_genes_pos`, which no rule reads), `engreitz_abc/…`,
  `gencode.v48lift37.annotation.gff3`, `hg19.fa.fai`, `hg38ToHg19.over.chain.gz`,
  `juicer_tools_1.22.01.jar`, `md5sum.txt`, `ncbi_dataset*`, `README.md` (NCBI boilerplate).
- Configs: `config_files/metadata_*.txt`, all `config/config_*.yaml` (merged into
  `config/config.yaml`).
- Rules: `RSL_count_TMP.smk` (draft), `merge_refs.smk` (old combined-library design),
  `merge_references.smk`, `Enformer.smk`, `enhancer_activity_score.smk` (empty/broken),
  `format_coordinates.smk` (own configfile, unused), `MAGECK.smk` (standalone with own
  `rule all`/configfile, not a module), `sra_fasterq_dump.smk`, `ABC_STARE_enhancer.smk`,
  `ABC_STARE_hic_matrix.smk` (Hi-C drafts superseded by hic_abc), `ConvertFiles.smk` (only
  fed the test track rule).
- Scripts not called by any rule: `average_SSMD_scores.py` (empty),
  `convert_input_to_JACKS.py`, `convert_to_bed.py`, `filter_aggregate_counts.py`,
  `get_sgRNA_OCR_coordinates.py`, `make_enhancer_activity_bed_SSMD_score.py`,
  `make_enhancer_activity_bed_unfiltered.py` (merged), `makeSampleInfo.py`,
  `merge_mageck_ref.py`, `merge_mageck_refs_alternative.py`, `parse_CrisprCounter.py`
  (only in the TMP draft), `workflow/utils/` (inlined into `plot_SSMD.py`).
- `workflow/envs/cisregvar.yaml` (→ `crispr.yaml` without pytorch), `python_env.yaml`
  (only used by `merge_refs`), `UCSC_tools.yaml` (→ `ucsc.yaml`).
- `workflow/Snakefile` `report: "report/workflow.rst"` (file does not exist).

**hic_abc**:
- `allEnrichments.tsv`, `dag.svg`, `output.png`, all `*.log`, `test_merged_nodups.pairs.gz`,
  `test_process_filter_CRISPRi_res.sh`, `test_sig_contacts.tsv`, `hicexplorer/`, `tmp/`,
  `.cache/`, `.config/`, `.pytest_cache/`, `.venv/`, `results/`.
- `workflow/notebooks/**` (all notebooks, CSVs, figures, gmt files).
- `resources/` (bwa-mem2 binaries, liftOver binary, juicer jar, LOLA RDS, `.sif`/`.tar`,
  chain file, gmt files) — downloaded or documented instead.
- `config/config_rwpe_5mln.yaml` (subsampling test → `hic.max_spots`),
  `config_rwpe.yaml`, `config_hmec.yaml` (merged).
- Rules merged/replaced: `sra_subsample.smk` (→ `sra` + `max_spots`), `TFenrich_OCR_universe.smk`
  (→ `TFenrich.smk`), `pairs2hic.smk` (was an incomplete `java` call), `set_HIC.smk`
  (logic → `common.smk`), `crispri_abc_significant_hits.smk` (was empty).

## 7. TODOs for you

- [ ] Fill in citation/DOI (README, `CITATION.cff`), repository URL, GEO accession of the
      screen FASTQs, and release of library/coordinate tables.
- [ ] Replace the inferred study summary with the abstract; complete the figure table.
- [ ] Build the envs once (`--conda-create-envs-only`) and fix any pin that does not solve.
- [ ] Test the Hi-C branch end-to-end on a subsample (`hic.max_spots: 5000000`), incl.
      `juicer_tools pre` on the new pairs and `Juicebox_normalization.sh`.
- [ ] Decide whether HMEC should be reprocessed (default) or taken from GEO (`set_HIC`)
      for the paper, and whether to keep MAPQ ≥ 30.
- [ ] Re-run ABC/STARE + downstream with the merged RWPE replicates and update figures.
- [ ] Consider adding cleaned figure notebooks (e.g. `workflow/notebooks/`) that read from
      the new `results/<CELL>/…` layout.
- [ ] Optionally enable `hic_qc` and report replicate SCC.
- [ ] Delete this file before publishing.
