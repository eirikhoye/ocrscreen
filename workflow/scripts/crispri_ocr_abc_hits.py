import pandas as pd
import numpy as np
import scipy.stats as stats
import argparse
import gzip
from dataclasses import dataclass, field


# ─────────────────────────────────────────────────────────────────────────────
#   Data container
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class CRISPRiData:
    """Holds all loaded DataFrames for one CRISPRi analysis run."""
    pathway_R1: pd.DataFrame = field(default_factory=pd.DataFrame)
    pathway_R2: pd.DataFrame = field(default_factory=pd.DataFrame)
    pathway: pd.DataFrame = field(default_factory=pd.DataFrame)

    ssmd_R1: pd.DataFrame = field(default_factory=pd.DataFrame)
    ssmd_R2: pd.DataFrame = field(default_factory=pd.DataFrame)
    ssmd: pd.DataFrame = field(default_factory=pd.DataFrame)

    coords: pd.DataFrame = field(default_factory=pd.DataFrame)
    ssmd_coords: pd.DataFrame = field(default_factory=pd.DataFrame)
    coords_unique: pd.DataFrame = field(default_factory=pd.DataFrame)

    pathway_coords: pd.DataFrame = field(default_factory=pd.DataFrame)
    pathway_coords_ocr: pd.DataFrame = field(default_factory=pd.DataFrame)
    pathway_coords_other: pd.DataFrame = field(default_factory=pd.DataFrame)

    abc: pd.DataFrame = field(default_factory=pd.DataFrame)


# ─────────────────────────────────────────────────────────────────────────────
#   Helper functions
# ─────────────────────────────────────────────────────────────────────────────

def split_ocr_other(df, column="Pathway", MATCH=r"^R\d+$"):
    """
    Split a DataFrame into OCR rows and non-OCR rows based on a regex match.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe.
    column : str
        Column to match against.
    MATCH : str
        Regex pattern; rows where *column* matches are treated as OCRs.

    Returns
    -------
    df_ocr : pd.DataFrame
        Rows where *column* matches the regex.
    df_other : pd.DataFrame
        Remaining rows.
    """
    mask = df[column].astype(str).str.match(MATCH)
    return df[mask].copy(), df[~mask].copy()


def filter_by_overlap(df, query_coord, chrom_col='#chr', start_col='start', end_col='end'):
    """
    Filter a DataFrame to rows whose genomic interval overlaps a query coordinate.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with genomic interval columns (e.g. ABC-STARE output).
    query_coord : str
        Coordinate string in "chr:start-end" format.
    chrom_col, start_col, end_col : str
        Column names for chromosome, start, and end positions.

    Returns
    -------
    pd.DataFrame
        Rows whose interval overlaps the query.
    """
    chrom, coords = query_coord.split(':')
    q_start, q_end = map(int, coords.split('-'))
    chrom = chrom.replace('chr', '')

    return df[
        (df[chrom_col].astype(str) == chrom) &
        (df[start_col] <= q_end) &
        (df[end_col] >= q_start)
    ]


def filter_by_overlap_list(df, query_coords, chrom_col='#chr', start_col='start', end_col='end'):
    """
    Filter a DataFrame to rows overlapping any of a list of query coordinates.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with genomic interval columns.
    query_coords : str or list[str]
        One or more coordinate strings in "chr:start-end" format.
    chrom_col, start_col, end_col : str
        Column names for chromosome, start, and end positions.

    Returns
    -------
    pd.DataFrame
        Deduplicated rows overlapping any query coordinate.
    """
    if isinstance(query_coords, str):
        query_coords = [query_coords]

    filtered_list = [
        filter_by_overlap(df, qc, chrom_col, start_col, end_col)
        for qc in query_coords
    ]
    return pd.concat(filtered_list).drop_duplicates().reset_index(drop=True)


def split_promoter_enhancer(df, tss_col='TSS-dist', promoter_window=2000):
    """
    Split contacts into promoter and enhancer sets based on TSS distance.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with a TSS-distance column.
    tss_col : str
        Column containing signed distance to the nearest TSS (bp).
    promoter_window : int
        Rows with |TSS-dist| <= this value are classified as promoters.

    Returns
    -------
    promoters : pd.DataFrame
        Rows within the promoter window.
    enhancers : pd.DataFrame
        Rows beyond the promoter window.
    """
    promoters = df[df[tss_col].abs() <= promoter_window].copy()
    enhancers = df[df[tss_col].abs() > promoter_window].copy()
    return promoters, enhancers


def _split_coord_string(coord):
    """
    Parse a "chr:start-end" or "N:start-end" coordinate string.

    Parameters
    ----------
    coord : str
        Coordinate string, e.g. "chr1:1000-2000" or "1:1000-2000".

    Returns
    -------
    chrom : str
        Chromosome without the 'chr' prefix.
    start : int
    end : int
    """
    chrom, positions = coord.split(':')
    start, end = map(int, positions.split('-'))
    chrom = chrom.replace('chr', '')
    return chrom, start, end


# ─────────────────────────────────────────────────────────────────────────────
#   Load data
# ─────────────────────────────────────────────────────────────────────────────

def load_data(
    pathway_R1_file,
    pathway_R2_file,
    ssmd_R1_file,
    ssmd_R2_file,
    coords_file,
    abc_file,
) -> CRISPRiData:
    """
    Load all input files and return a populated CRISPRiData container.

    Parameters
    ----------
    pathway_R1_file : str
        Path to Replicate 1 MAGeCK pathway TSV.
    pathway_R2_file : str
        Path to Replicate 2 MAGeCK pathway TSV.
    ssmd_R1_file : str
        Path to Replicate 1 SSMD scores TSV.
    ssmd_R2_file : str
        Path to Replicate 2 SSMD scores TSV.
    coords_file : str
        Path to sgRNA / OCR coordinates CSV (columns: sgRNA, open_reg, Coord_RWPE, …).
    abc_file : str
        Path to gzipped ABC model scores TSV.

    Returns
    -------
    CRISPRiData
    """
    data = CRISPRiData()

    data.pathway_R1 = pd.read_csv(pathway_R1_file, sep='\t')
    data.pathway_R2 = pd.read_csv(pathway_R2_file, sep='\t')

    data.pathway = pd.merge(
        data.pathway_R1,
        data.pathway_R2,
        on='Pathway',
        suffixes=("_R1", "_R2"),
    )

    data.ssmd_R1 = pd.read_csv(ssmd_R1_file, sep='\t')
    data.ssmd_R2 = pd.read_csv(ssmd_R2_file, sep='\t')

    data.ssmd = pd.merge(
        data.ssmd_R1,
        data.ssmd_R2,
        on='guide.set',
        suffixes=("_R1", "_R2"),
    )

    data.coords = pd.read_csv(coords_file)

    data.ssmd_coords = pd.merge(
        data.ssmd,
        data.coords,
        left_on='guide.set',
        right_on='sgRNA',
        how='inner',
    )

    # One row per OCR (first occurrence kept; all guides for the same OCR share
    # the same open_reg coordinates, so the choice is arbitrary).
    data.coords_unique = data.coords.drop_duplicates(subset='open_reg')

    data.pathway_coords = pd.merge(
        data.pathway,
        data.coords_unique,
        left_on='Pathway',
        right_on='open_reg',
        how='inner',
    )

    data.pathway_coords_ocr, data.pathway_coords_other = split_ocr_other(
        data.pathway_coords
    )

    with gzip.open(abc_file) as f:
        data.abc = pd.read_csv(f, sep='\t')

    return data


# ─────────────────────────────────────────────────────────────────────────────
#   Fisher's combined p-value across replicates
# ─────────────────────────────────────────────────────────────────────────────

def _fisher_row(row):
    """Apply Fisher's method to combine p-values from two replicates (one row)."""
    ps = [max(row['p_R1'], 1e-300), max(row['p_R2'], 1e-300)]
    chi2 = -2 * sum(np.log(ps))
    return stats.chi2.sf(chi2, df=2 * len(ps))


def compute_concordance(data: CRISPRiData, coord_col='Coord_RWPE'):
    """
    Compute replicate concordance and Fisher's combined p-values for all OCRs.

    Prints QC summaries to stdout (captured by Snakemake logs).

    Parameters
    ----------
    data : CRISPRiData
        Populated data container from load_data().

    Returns
    -------
    concordance : pd.DataFrame
        One row per OCR with ES, FDR, p-values, concordance tier, mean ES,
        Fisher combined p, BH-adjusted Fisher p, and -log10(padj).
    concordance_coord : pd.DataFrame
        concordance joined with genomic coordinates and a PeakID column.
    """
    print(f"  Pathway R1:  {data.pathway_R1.shape}")
    print(f"  Pathway R2:  {data.pathway_R2.shape}")
    print(f"  Library:     {data.coords.shape}")
    print(f"  ABC:         {data.abc.shape}")
    print(f"  SSMD R1:     {data.ssmd_R1.shape}")
    print(f"  SSMD R2:     {data.ssmd_R2.shape}")

    # ── QC ────────────────────────────────────────────────────────────────────
    print("\n── QC ──")
    print(f"ABC duplicates: {data.abc.duplicated().sum()}")
    print(f"Library duplicates (sgRNA): {data.coords.duplicated(subset='sgRNA').sum()}")
    print(f"ABC missing values:\n{data.abc.isnull().sum()[data.abc.isnull().sum() > 0]}")
    print(f"Library missing values:\n{data.coords.isnull().sum()[data.coords.isnull().sum() > 0]}")
    print(f"ABC score range: {data.abc['ABC-Score'].min():.4f} – {data.abc['ABC-Score'].max():.4f}")
    print(f"All ABC scores >= 0.02 (pre-filtered): {(data.abc['ABC-Score'] >= 0.02).all()}")

    OCR_PATTERN = r'^R\d+$'

    p1_ocr = data.pathway_R1[
        data.pathway_R1['Pathway'].astype(str).str.match(OCR_PATTERN)
    ].copy()
    p2_ocr = data.pathway_R2[
        data.pathway_R2['Pathway'].astype(str).str.match(OCR_PATTERN)
    ].copy()

    p1_sig = p1_ocr[p1_ocr['FDR'] < 0.05]
    p2_sig = p2_ocr[p2_ocr['FDR'] < 0.05]

    r1_sig_set = set(p1_sig['Pathway'])
    r2_sig_set = set(p2_sig['Pathway'])

    print(f"\nSignificant OCRs (FDR<0.05):")
    print(f"  R1 only:      {len(r1_sig_set - r2_sig_set)}")
    print(f"  R2 only:      {len(r2_sig_set - r1_sig_set)}")
    print(f"  Both (inter): {len(r1_sig_set & r2_sig_set)}")
    print(f"  Union:        {len(r1_sig_set | r2_sig_set)}")

    # ── Build concordance table ────────────────────────────────────────────────
    p1_ocr_all = p1_ocr.rename(columns={
        'ES': 'ES_R1', 'FDR': 'FDR_R1', 'p': 'p_R1',
        'p_permutation': 'p_perm_R1', 'Size': 'Size_R1',
        'Ranking': 'Ranking_R1', 'Hits': 'Hits_R1', 'LFC': 'LFC_R1',
    })
    p2_ocr_all = p2_ocr.rename(columns={
        'ES': 'ES_R2', 'FDR': 'FDR_R2', 'p': 'p_R2',
        'p_permutation': 'p_perm_R2', 'Size': 'Size_R2',
        'Ranking': 'Ranking_R2', 'Hits': 'Hits_R2', 'LFC': 'LFC_R2',
    })

    concordance = pd.merge(
        p1_ocr_all[['Pathway', 'ES_R1', 'FDR_R1', 'p_R1']],
        p2_ocr_all[['Pathway', 'ES_R2', 'FDR_R2', 'p_R2']],
        on='Pathway',
        how='outer',
    )

    def _concordance_tier(row):
        r1_sig = row['FDR_R1'] < 0.05
        r2_sig = row['FDR_R2'] < 0.05
        if r1_sig and r2_sig:
            return 'Both replicates'
        elif r1_sig:
            return 'R1 only'
        elif r2_sig:
            return 'R2 only'
        else:
            return 'none'

    concordance['concordance_tier'] = concordance.apply(_concordance_tier, axis=1)
    concordance['mean_ES'] = concordance[['ES_R1', 'ES_R2']].mean(axis=1)
    concordance['fisher_p'] = concordance.apply(_fisher_row, axis=1)
    concordance['fisher_p_adj'] = stats.false_discovery_control(
        concordance['fisher_p'], method='bh'
    )
    concordance = concordance.sort_values(by='fisher_p_adj')
    concordance['-log10(padj)'] = -np.log10(concordance['fisher_p_adj'])

    print(f"\nConcordance table: {concordance.shape}")
    print(concordance['concordance_tier'].value_counts())
    print(f"\nFisher combined p-adj < 0.05: {(concordance['fisher_p_adj'] < 0.05).sum()}")

    # ── Add genomic coordinates ────────────────────────────────────────────────
    concordance_coord = pd.merge(
        data.pathway_coords[['Pathway', coord_col]],
        concordance,
        on='Pathway',
    )

    concordance_coord[['#chr', 'start', 'end']] = pd.DataFrame(
        concordance_coord[coord_col].apply(_split_coord_string).tolist(),
        index=concordance_coord.index,
    )

    concordance_coord['PeakID'] = (
        concordance_coord['#chr'].astype(str) + ':' +
        concordance_coord['start'].astype(str) + '-' +
        concordance_coord['end'].astype(str)
    )

    return concordance, concordance_coord


# ─────────────────────────────────────────────────────────────────────────────
#   Merge significant OCRs with ABC contacts
# ─────────────────────────────────────────────────────────────────────────────

def merge_filtered_ocrs_with_abc(
    data: CRISPRiData,
    concordance: pd.DataFrame,
    concordance_coord: pd.DataFrame,
    coord_col: str = 'Coord_RWPE',
):
    """
    Filter OCRs to those significant in both replicates and the Fisher test,
    then join with ABC model contacts and annotate as promoter or enhancer.

    Parameters
    ----------
    data : CRISPRiData
        Populated data container.
    concordance : pd.DataFrame
        Output of compute_concordance().
    concordance_coord : pd.DataFrame
        Output of compute_concordance().

    Returns
    -------
    sig_df : pd.DataFrame
        Significant OCR–gene contacts with a 'type' column ('promoter'/'enhancer').
    pathway_coords_ocr_sig : pd.DataFrame
        The significant OCRs (one row per OCR) with their coordinates.
    """
    ocr_sig_fisher_filter = concordance[
        (concordance['fisher_p_adj'] < 0.05) &
        (concordance['p_R1'] < 0.05) &
        (concordance['p_R2'] < 0.05)
    ]['Pathway']

    pathway_coords_ocr_sig = data.pathway_coords_ocr[
        data.pathway_coords_ocr['Pathway'].isin(ocr_sig_fisher_filter)
    ]

    sig_contacts = filter_by_overlap_list(data.abc, pathway_coords_ocr_sig[coord_col])
    sig_contacts = pd.merge(sig_contacts, concordance_coord, on='PeakID')
    sig_contacts = sig_contacts.drop(
        columns=['#chr_y', 'start_y', 'end_y'], errors='ignore'
    ).rename(
        columns={'#chr_x': '#chr', 'start_x': 'start', 'end_x': 'end'}
    )

    sig_promoters, sig_enhancers = split_promoter_enhancer(sig_contacts)
    sig_promoters = sig_promoters.copy()
    sig_enhancers = sig_enhancers.copy()
    sig_promoters['type'] = 'promoter'
    sig_enhancers['type'] = 'enhancer'

    sig_df = pd.concat([sig_promoters, sig_enhancers], axis=0, ignore_index=True)
    return sig_df, pathway_coords_ocr_sig


def write_sig_ocr_bed(pathway_coords_ocr_sig: pd.DataFrame, coord_col: str, output_file: str):
    """
    Write the significant OCRs (CREs) as a 4-column BED file (chrom, start, end, name).

    This reproduces the "*_sig_CRE_hg19.bed" files used as input for the UniBind
    TF-enrichment analysis.
    """
    bed = pathway_coords_ocr_sig[coord_col].str.extract(r"(chr[\w]+):(\d+)-(\d+)")
    bed.columns = ['chrom', 'start', 'end']
    bed['start'] = bed['start'].astype(int)
    bed['end'] = bed['end'].astype(int)
    bed['name'] = pathway_coords_ocr_sig['Pathway']
    bed.to_csv(output_file, sep='\t', header=False, index=False)


# ─────────────────────────────────────────────────────────────────────────────
#   Entry point
# ─────────────────────────────────────────────────────────────────────────────

def main():
    """
    Process and filter MAGeCK pathway enrichment and SSMD scores using
    Fisher's combined p-values from the pathway output, apply FDR correction,
    and remove OCRs with an original p-value >= 0.05 in either replicate.

    Input files (positional arguments):
        pathway_R1   -- Replicate 1 MAGeCK pathway TSV
        pathway_R2   -- Replicate 2 MAGeCK pathway TSV
        ssmd_R1      -- Replicate 1 SSMD scores TSV
        ssmd_R2      -- Replicate 2 SSMD scores TSV
        coords       -- sgRNA / OCR coordinates CSV (columns: sgRNA, open_reg, Coord_<CELL>)
        abc          -- Gzipped ABC model scores TSV
        output_file  -- Output path for significant OCR–gene contacts TSV

    Output:
        Tab-separated file with significant OCR–gene contacts annotated as
        'promoter' or 'enhancer', including Fisher combined p-values and
        ABC model scores.
    """
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('pathway_R1',  type=str, help="Replicate 1 MAGeCK pathway TSV")
    parser.add_argument('pathway_R2',  type=str, help="Replicate 2 MAGeCK pathway TSV")
    parser.add_argument('ssmd_R1',     type=str, help="Replicate 1 SSMD scores TSV")
    parser.add_argument('ssmd_R2',     type=str, help="Replicate 2 SSMD scores TSV")
    parser.add_argument('coords',      type=str, help="sgRNA / OCR coordinates CSV")
    parser.add_argument('abc',         type=str, help="Gzipped ABC model scores TSV")
    parser.add_argument('output_file', type=str, help="Output significant OCR–gene contacts TSV")
    parser.add_argument('--coord-col', default='Coord_RWPE',
                        help="Coordinate column in the coords CSV (e.g. Coord_RWPE, Coord_HMEC)")
    parser.add_argument('--sig-bed', default=None,
                        help="Optional output BED of significant OCRs (input for TF enrichment)")
    args = parser.parse_args()

    data = load_data(
        args.pathway_R1,
        args.pathway_R2,
        args.ssmd_R1,
        args.ssmd_R2,
        args.coords,
        args.abc,
    )

    concordance, concordance_coord = compute_concordance(data, coord_col=args.coord_col)
    sig_df, sig_ocrs = merge_filtered_ocrs_with_abc(
        data, concordance, concordance_coord, coord_col=args.coord_col
    )
    sig_df.to_csv(args.output_file, index=False, sep='\t')
    if args.sig_bed:
        write_sig_ocr_bed(sig_ocrs, args.coord_col, args.sig_bed)


if __name__ == '__main__':
    main()
