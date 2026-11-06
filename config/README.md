# Configuration

`config/config.yaml` is the full analysis of the paper (RWPE-1 + HMEC, all modules).
It is validated against `workflow/schemas/config.schema.yaml` at start-up. Relative
paths are relative to the repository root (the directory you run `snakemake` from).
Several config files can be combined (`--configfile a.yaml b.yaml`); later files
override earlier ones, lists are replaced.

## Sample sheets

### `crispr_samples.tsv` — one row per screen FASTQ

| column | meaning |
|---|---|
| `sample` | screen replicate, `<cell_line>_<replicate>` (e.g. `RWPE_R1`); the unit of RSL counting/SSMD/pathway |
| `cell_line` | key in `cell_lines` |
| `replicate` | replicate label (`R1`, `R2`) |
| `stage` | time point, one of `crispr.stages` |
| `label` | CrisprCounter sample name = count-table column, `<sample>_<stage>` |
| `group` | CrisprCounter group (`rep1`, `rep2`) |
| `control` | JACKS control sample, `<cell_line>_<reference stage>` |
| `fastq` | guide-read FASTQ name (`*_R1_001.fastq.gz`) in `cell_lines.<CELL>.fastq_dir`; the RSL read is the matching `*_R2_001.fastq.gz` |

Steps that compare replicates (`plot_SSMD`, `crispri_abc_significant_hits`) expect exactly
two screen replicates per cell line. `plot_SSMD.py` and `ssmd2bed.py` use the column
names `<sample>_D14.*`, i.e. the late time point must be called `D14`.

### `hic_samples.tsv` — one row per Hi-C sequencing run

| column | meaning |
|---|---|
| `cell_line` | key in `cell_lines` |
| `replicate` | library / biological replicate. Runs with the same `replicate` are merged **before** deduplication; replicates are deduplicated separately and merged afterwards |
| `run` | SRA run accession (downloaded by module `sra`) |
| `geo_sample`, `description` | documentation only |

## Keys

### Top level

| key | default | meaning |
|---|---|---|
| `modules` | all except `set_HIC`, `hic_qc` | modules to run, see the module table in the top-level README |
| `results` | `results` | output directory |
| `resources` | `resources` | directory for downloaded references/tools |
| `samples.crispr` | `config/crispr_samples.tsv` | CRISPRi sample sheet |
| `samples.hic` | `config/hic_samples.tsv` | Hi-C sample sheet |

### `cell_lines.<CELL>`

| key | meaning |
|---|---|
| `fastq_dir` | directory with the screen FASTQs of this cell line |
| `library` | sgRNA library CSV without header: `sgRNA_id,protospacer(N20),OCR_id` |
| `coords` | OCR coordinate table (tab-separated) with columns `N20`, `Coord_<CELL>` (`chr:start-end`, hg19), `Guide ID`, `Peak ID` |
| `neg_con` | prefix of negative-control guide sets (`NCR` for RWPE, `NCH` for HMEC); used by SSMD, GMT construction, JACKS (`--ctrl_genes`) and plotting |
| `dnase_peaks` | DNase-seq narrowPeak (hg19). Activity for ABC run `dnase` and universe for `TFenrich_OCR_universe`. Files named `resources/encode/wgEncodeOpenChromDnase*Pk.narrowPeak` are downloaded from UCSC |
| `precomputed_hic` | (module `set_HIC`) `.hic` file used instead of processing reads; KR contacts are dumped from it by `hic2KR` |
| `precomputed_contacts` | (module `set_HIC`) directory of KR contacts (`chr<N>_KR_Contacts.txt.gz`) used directly by STARE |

### `reference`

| key | meaning |
|---|---|
| `build` | genome build label (`hg19`); used for file names and `--assembly` |
| `fasta_url` | genome FASTA (gzip) → `resources/genome/<build>.fa` |
| `gtf_url` | gene annotation (gzip) → `resources/genome/<basename>.gtf`; a protein-coding-genes subset is derived from it |
| `chain_hg19_to_hg38_url` | liftOver chain for UniBind (hg38) |
| `dnase_url_template` | URL template for `download_dnase_peaks` (`{name}` = file stem) |

### `crispr`

| key | default | meaning |
|---|---|---|
| `stages` | `[D4, D14]` | reference and late time point |
| `crisprcounter.jar` | `resources/tools/CrisprCounter.jar` | CrisprCounter (not redistributable) |
| `crisprcounter.template` | `resources/crisprcounter_template.properties` | properties template (`truncumilen = 10`, i.e. full RSL) |
| `crisprcounter.java_mem_gb` | 32 | Java heap |
| `rsl_bin_nt` | 3 | RSLs are binned by their first `n` bases (4ⁿ bins); bins are the JACKS synthetic replicates. The paper used 3 (older runs: 1 and 2) |
| `rsl_bin_min_count` | 0 | keep an RSL if control or treatment count ≥ this |
| `ssmd.min_control_count` | 1 | IRA-SSMD argument 3: minimum control count (old key `SSMD_trunclen`) |
| `ssmd.filter_on_sum` | 0 | IRA-SSMD argument 4: 1 = filter on control+treatment sum instead of control (old key `SSMD_mincount`) |
| `positive_controls` | `resources/positive_control_genes.txt` | positive-control genes highlighted in the SSMD figure |
| `pathway.permutations` | 1000000 | `mageck pathway --permutation` |
| `pathway.ranking_column` | 7 | `mageck pathway --ranking-column` (column of the SSMD table) |

### `jacks`

| key | default | meaning |
|---|---|---|
| `run_jacks_url` | JACKS commit dd5c4be | `run_JACKS.py` download URL |
| `ctrl_sample` | `Control` | `--ctrl_sample` |
| `sgrna_hdr`, `gene_hdr` | `sgRNA`, `open_reg` | library column names |
| `fdr_threshold` | 0.05 | BH-FDR threshold for the filtered enhancer-activity BED and for module `GSEA` |

### `hic`

| key | default | meaning |
|---|---|---|
| `max_spots` | null | if set, only the first N read pairs per run are downloaded (`fastq-dump`; testing) |
| `sra_max_size` | `200g` | `prefetch --max-size` |
| `threads.bwa`, `threads.sort`, `threads.pairtools` | 32, 16, 8 | threads per step |
| `parse_options` | see file | `pairtools parse` options (MAPQ ≥ 30, `5unique` walks, mapq columns) |
| `dedup_max_mismatch` | 3 | `pairtools dedup --max-mismatch` (bp) |
| `cooler_base_resolution` | 1000 | base resolution of the merged `.cool` |
| `zoomify_resolutions` | 1 kb … 1 Mb | resolutions of the balanced `.mcool` |
| `qc_resolution` | 100000 | resolution of the per-library QC `.cool` (and HiCRep) |
| `hicrep.h`, `hicrep.dbpmax` | 1, 5000000 | HiCRep smoothing and max. distance (module `hic_qc`) |
| `juicer.jar`, `juicer.url` | juicer_tools 1.22.01 | jar path and download URL |
| `juicer.memory_mb`, `juicer.threads` | 48000, 16 | `juicer_tools pre` resources |
| `kr.norm`, `kr.bin_size`, `kr.chromosomes` | KR, 5000, 1–22,X | contact dump for STARE (`Juicebox_normalization.sh`); `bin_size` is also STARE `-k` |
| `precomputed_hic_urls` | HMEC combined map | path → URL for precomputed `.hic` files under `resources/hic/` |

### `abc_stare`

| key | meaning |
|---|---|
| `psem_url` | STARE PSEM (JASPAR/HOCOMOCO/Kellis, STARE PWMs 2.2) |
| `threads` | STARE cores |
| `runs.<name>.activity` | `dnase` (cell-line DNase peaks), `jacks_fdr` (JACKS-significant OCRs), `jacks_all` (all OCRs) |
| `runs.<name>.column` | STARE `-n` (activity column(s); `7` = narrowPeak signalValue, `5+` = JACKS score column) |
| `runs.<name>.tss` | STARE `-i` (`all_tss` or `5_tss`) |
| `runs.<name>.gtf` | `all` or `protein_coding` |

### Downstream modules

| key | meaning |
|---|---|
| `crispri_hits.abc_run` | ABC run intersected with the significant OCRs (paper: `dnase`) |
| `gsea.abc_run` | ABC run used for gene ranking (original: JACKS-filtered activity, protein-coding genes) |
| `gsea.permutations`, `gsea.gene_sets` | GSEApy prerank permutations and Enrichr libraries |
| `genome_tracks.regions` | name → region (`chr:start-end`) plotted by `GenomeTracks` |
| `encode` | file name → URL downloaded by `dl_encode` into `resources/encode/tracks/` |
| `unibind.enrich_script` | `bin/UniBind_enrich.sh` of a unibind_enrichment checkout |
| `unibind.container` | Singularity image (`--use-singularity`); `null` to run without container |
| `unibind.lola_db`, `unibind.lola_universe` | UniBind LOLA databases (Zenodo 4161027) |

## Mapping from the original pipelines

| original (crispr_pipeline / hic_abc) | here |
|---|---|
| `enabled_modules` / `modules` | `modules` |
| `add_coords` | `crispr_library` |
| `RUN_GSEA` | `GSEA` |
| `sra_subsample` (`nr_reads`) | `sra` with `hic.max_spots` |
| `samples` (nested YAML) | `config/crispr_samples.tsv` |
| `samples` (SRR list) | `config/hic_samples.tsv` |
| `nr_nt_RSL`, `bin_min_count` | `crispr.rsl_bin_nt`, `crispr.rsl_bin_min_count` |
| `SSMD_trunclen`, `SSMD_mincount` | `crispr.ssmd.min_control_count`, `crispr.ssmd.filter_on_sum` |
| `Neg_Con`, `ctrl_genes_neg` | `cell_lines.<CELL>.neg_con` |
| `mageck_pathway_permutations` | `crispr.pathway.permutations` |
| `hic_processed` | `cell_lines.<CELL>.precomputed_hic` |
| `enhancer_activity`, `all_ocr_narrowpeak` | `cell_lines.<CELL>.dnase_peaks` |
| `juicebox_norm`, `stare_path`, `rmtree_path`, `*_env` | provided by `workflow/envs/*.yaml` / `workflow/scripts` |
