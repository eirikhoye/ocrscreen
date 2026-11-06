import argparse
import gseapy as gp
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Run GSEApy preranked analysis on a ranked gene list"
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Ranked gene list file (.rnk, two columns: gene, score, no header)"
    )
    parser.add_argument(
        "-o", "--outdir", required=True,
        help="Output directory for GSEA results"
    )
    parser.add_argument(
        "-g", "--geneset", default="KEGG_2021_Human",
        help="Gene set library (default: KEGG_2021_Human)"
    )
    parser.add_argument(
        "-p", "--permutations", type=int, default=1000,
        help="Number of permutations (default: 1000)"
    )
    args = parser.parse_args()

    # Read ranked list
    rnk = pd.read_csv(args.input, sep="\t", header=None)

    res_file = f"{args.outdir}/gsea_results.tsv"

    # Run GSEApy preranked
    pre_res = gp.prerank(
        rnk=rnk,
        gene_sets=args.geneset,
        outdir=args.outdir,
        permutation_num=args.permutations,
        seed=42
    )

    # Print summary to stdout
    pre_res.res2d.to_csv(res_file, sep="\t")


if __name__ == "__main__":
    main()
