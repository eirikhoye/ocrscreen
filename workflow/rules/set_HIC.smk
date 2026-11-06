"""
Module set_HIC (alternative to sra/bwa/pairtools/pairs2hic).

Use a precomputed Hi-C map for selected cell lines, configured per cell line:

  cell_lines.<CELL>.precomputed_hic:      .hic file -> KR contacts made by hic2KR
  cell_lines.<CELL>.precomputed_contacts: directory of KR contacts used directly by STARE

Cell lines with a precomputed map are skipped by the Hi-C processing modules; the
selection logic lives in common.smk (uses_precomputed_hic, hic_file, contacts_dir).
Published maps listed in hic.precomputed_hic_urls are downloaded by
rule download_precomputed_hic (resources.smk).
"""
