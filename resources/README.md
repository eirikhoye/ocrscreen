# resources/

Reference data and external tools. Only two small text files are tracked in git;
everything else is downloaded by the workflow or must be provided by the user.

## Tracked files

| file | origin |
|---|---|
| `crisprcounter_template.properties` | CrisprCounter properties template from [CRISPR-pooled-RSL](https://github.com/agata-sm/CRISPR-pooled-RSL) (`assets/template.properties`, MIT, © 2022 Agata Smialowska) |
| `positive_control_genes.txt` | positive-control gene targets of the libraries (ribosomal protein genes and tumour-suppressor / Hippo-pathway genes such as TP53, CDKN1A, RB1, NF2, LATS2); used to colour the SSMD figure (`crispr.positive_controls`). Copied from crispr_pipeline `resources/ctrl_genes_single.txt` |

## Downloaded automatically (when a module needs them)

| path | rule | source |
|---|---|---|
| `genome/hg19.fa`, `.fai`, `hg19.chrom.sizes` | `download_genome`, `chrom_sizes` | UCSC hg19 |
| `genome/bwa-mem2/hg19.*` | `bwa_mem2_index` | built locally |
| `genome/gencode.v48lift37.annotation.gtf` (+ `.protein_coding.genes.gtf`) | `download_gtf`, `filter_gtf_protein` | GENCODE release 48, GRCh37 mapping |
| `genome/hg19ToHg38.over.chain.gz` | `download_chain` | UCSC liftOver |
| `stare/Jaspar_Hocomoco_Kellis_human_transfac.txt` | `download_psem` | STARE repository, PWMs/2.2 |
| `encode/wgEncodeOpenChromDnase{Rwpe1,Hmec}Pk.narrowPeak` | `download_dnase_peaks` | UCSC `wgEncodeOpenChromDnase` (ENCODE/Duke) |
| `encode/tracks/*` | `download_encode_track` (module `dl_encode`) | ENCODE portal |
| `tools/juicer_tools_1.22.01.jar` | `get_juicer_tools` | Aiden lab |
| `tools/run_JACKS.py` | `get_run_jacks` | JACKS GitHub, commit dd5c4be |
| `hic/GSE63525_HMEC_combined.hic` | `download_precomputed_hic` (only with `set_HIC`) | GEO GSE63525 |

## Must be provided

| path (default in `config/config.yaml`) | what |
|---|---|
| `crispr/LibraryFile_RWPE.csv`, `crispr/LibraryFile_Me16C.csv` | sgRNA libraries (no header: `sgRNA_id,N20,OCR_id`); 3–4 MB each, `<to be released with the paper: supplementary table / Zenodo TBD>` |
| `crispr/Final_List_RWPE.csv`, `crispr/Final_List_HMEC.csv` | OCR coordinate tables (tab-separated; columns incl. `N20`, `Coord_RWPE`/`Coord_HMEC`, `Guide ID`, `Peak ID`), ~10–12 MB each, same release |
| `data/crispr/<CELL>/*.fastq.gz` (`cell_lines.<CELL>.fastq_dir`) | demultiplexed screen FASTQs (`*_R1_001` guide read, `*_R2_001` RSL read), `<GEO accession TBD>` |
| `tools/CrisprCounter.jar` | CrisprCounter (T. Kivioja, University of Helsinki). Not redistributable; contact the CRISPR-pooled-RSL authors |
| `unibind_enrichment/` | `git clone https://bitbucket.org/CBGR/unibind_enrichment.git resources/unibind_enrichment` |
| `unibind_enrichment_0.2.sif` | Apptainer/Singularity image with the unibind_enrichment dependencies (`unibind.container`); alternatively deploy the tool locally (`./deploy.sh`) and set `unibind.container: null` |
| `unibind/hg38_robust_UniBind_LOLA.RDS`, `unibind/hg38_robust_UniBind_LOLA_universe.RDS` | UniBind LOLA databases, Zenodo doi:10.5281/zenodo.4161027 (e.g. `resources/unibind_enrichment/bin/zenodo_fetch 4161027 -f hg38_robust_UniBind_LOLA.RDS -f hg38_robust_UniBind_LOLA_universe.RDS`) |
