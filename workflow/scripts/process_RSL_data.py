import pandas as pd
import numpy as np
import argparse

def process_RSL_data(input_file, output_file, sample_col1, sample_col2, pseudocount=1):
    """
    Processes an input dataframe of RSL data, filters, collapses by guide,
    and computes log2 fold changes with generalizable sample column names.

    Parameters:
    - input_file (str): Path to the input CSV file.
    - output_file (str): Path to save the processed output CSV file.
    - sample_col1 (str): The name of the first sample column.
    - sample_col2 (str): The name of the second sample column.
    - pseudocount (int): Pseudocount to add to avoid division by zero (default=1).
    """
    # Load data
    RSL_counts = pd.read_csv(input_file)

    # Step 1: Filter rows where the sum of counts is >= 2
    RSL_counts['row_sum'] = RSL_counts[[sample_col1, sample_col2]].sum(axis=1)
    RSL_filtered = RSL_counts[RSL_counts['row_sum'] >= 2].drop(columns=['row_sum']).reset_index(drop=True)

    # Step 2: Collapse data (convert to binary presence/absence for each sample)
    RSL_collapsed = RSL_filtered.copy()
    RSL_collapsed[sample_col1] = RSL_collapsed[sample_col1].apply(lambda x: 1 if x >= 1 else 0)
    RSL_collapsed[sample_col2] = RSL_collapsed[sample_col2].apply(lambda x: 1 if x >= 1 else 0)

    # Step 3: Group by 'guide.set' and count RSLs per sample
    collapsed_RSL_counts = RSL_collapsed.groupby('guide.set')[[sample_col1, sample_col2]].sum().reset_index()

    # Step 4: Adjust counts to exclude a specific RSL and add pseudocount
    collapsed_RSL_counts.iloc[:, 1:] -= 1
    collapsed_RSL_counts.iloc[:, 1:] += pseudocount

    # Step 5: Calculate log2 fold change
    collapsed_RSL_counts['log2_fold_change'] = np.log2(
        collapsed_RSL_counts[sample_col2] / collapsed_RSL_counts[sample_col1]
    )

    # Step 6: Extract open chromatin region from 'guide.set' and create a new column
    collapsed_RSL_counts['open_chromatin_region'] = collapsed_RSL_counts['guide.set'].apply(lambda x: x.split('_')[0])

    # Step 7: Reorder columns to place 'open_chromatin_region' after 'guide.set'
    column_order = ['guide.set', 'open_chromatin_region', sample_col1, sample_col2, 'log2_fold_change']
    collapsed_RSL_counts = collapsed_RSL_counts[column_order]

    # Step 8: Rank by log2_fold_change in descending order
    collapsed_RSL_counts = collapsed_RSL_counts.sort_values(by='log2_fold_change', ascending=False).reset_index(drop=True)

    # Save the result to the output file as a tab-separated file
    collapsed_RSL_counts.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process RSL data and calculate log2 fold changes.")
    parser.add_argument("input_file", help="Path to the input CSV file.")
    parser.add_argument("output_file", help="Path to the output CSV file.")
    parser.add_argument("sample_col1", help="Name of the first sample column.")
    parser.add_argument("sample_col2", help="Name of the second sample column.")
    parser.add_argument("--pseudocount", type=int, default=1, help="Pseudocount to add (default: 1).")

    args = parser.parse_args()

    process_RSL_data(
        input_file=args.input_file,
        output_file=args.output_file,
        sample_col1=args.sample_col1,
        sample_col2=args.sample_col2,
        pseudocount=args.pseudocount,
    )
