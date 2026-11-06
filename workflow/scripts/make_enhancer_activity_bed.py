import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
import argparse


"""
Process results and compute FDR-adjusted p-values.

Usage:
    python myscript.py \
        --results df_results.tsv \
        --pvals raw_pvals.tsv \
        --coord_ref coordinates.bed \
        --output results_with_fdr.tsv

Arguments:
    -r, --results     Path to the results file (e.g., with gene scores).
    -p, --pvals       Path to the file containing p-values.
    -c, --coord_ref   Path to the coordinate reference file.
    -o, --output      Path to the output file.
"""


def get_args():
    parser = argparse.ArgumentParser(
        description="Process results and compute FDR-adjusted p-values."
        )
    parser.add_argument(
        '-r', "--results",
        required=True,
        help="Path to the results file (with JACKS gene scores)."
    )

    parser.add_argument(
        '-p', "--pvals",
        required=True,
        help="Path to the file containing p-values."
    )

    parser.add_argument(
        '-c', "--coord_ref",
        required=True,
        help="Path to the coordinate reference file."
    )

    parser.add_argument(
        '-n', "--cell_line",
        required=True, 
        help="name of cell line column name"
    )

    parser.add_argument(
        '-o', "--output",
        required=True,
        help="Path to the output file."
    )

    parser.add_argument(
        '-t', "--fdr_threshold",
        type=float,
        default=0.05,
        help="FDR threshold to use (default: 0.05)."
    )

    parser.add_argument(
        "--no_fdr_filter",
        action="store_true",
        help="Keep all OCRs (do not filter on the FDR threshold)."
    )

    return parser.parse_args()


def make_enhancer_activity_bed(results, 
                               pvals, 
                               coord_ref, 
                               cell_line, 
                               output, 
                               fdr_threshold=0.05,
                               fdr_filter=True):
    df_result = pd.read_csv(results, sep='\t')
    df_pvals = pd.read_csv(pvals, sep='\t')
    df_coord = pd.read_csv(coord_ref)

    df_result = df_result.rename(columns={col: f"{col}_score" for col in df_result.columns if col != "Gene"})
    df_pvals = df_pvals.rename(columns={col: f"{col}_pval" for col in df_pvals.columns if col != "Gene"})
    df_result = pd.merge(df_result, df_pvals, on='Gene', how='inner')

    pvals = df_result[f'{cell_line}_D14_pval']
    _, fdr, _, _ = multipletests(pvals, method='fdr_bh')

    df_result[f'{cell_line}_D14_fdr'] = fdr
    if fdr_filter:
        df_result = df_result[df_result[f"{cell_line}_D14_fdr"] < fdr_threshold]

    df_coord_unique = df_coord[['sequence', 'open_reg', f'Coord_{cell_line}']].drop_duplicates(subset='open_reg').set_index('open_reg')

    def parse_coord(coord_str):
        try:
            chrom, pos = coord_str.split(':')
            start, end = pos.split('-')
            return chrom, int(start), int(end)
        except:
            return np.nan, np.nan, np.nan
    
    rows = []

    for _, row in df_result.iterrows():
        coord = df_coord_unique.loc[row.Gene][f'Coord_{cell_line}']
        if pd.isna(coord):
            continue

        if not pd.isna(coord):
            chrom, start, end = parse_coord(coord)
        
        rows.append([chrom, start, end, row.Gene, np.abs(row[f'{cell_line}_D14_score'])])

    df_bed = pd.DataFrame(rows, columns=['chr', 'start', 'end', 'OCR', f'{cell_line}_activity'])
    df_bed.columns = ['#chr', 'start', 'end', 'OCR', f'{cell_line}']
    df_bed.to_csv(output, sep='\t', index=False)

def main():
    args = get_args()
    make_enhancer_activity_bed(
        results=args.results,
        pvals=args.pvals,
        coord_ref=args.coord_ref,
        cell_line=args.cell_line,
        output=args.output,
        fdr_threshold=args.fdr_threshold,
        fdr_filter=not args.no_fdr_filter
    )

if __name__ == "__main__":
    main()