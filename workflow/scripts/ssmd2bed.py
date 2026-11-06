import pandas as pd
import argparse


def get_args():
    parser = argparse.ArgumentParser(
        description="Convert SSMD scores to BED format (sgRNA or OCR level)."
    )

    parser.add_argument("-r", "--results", required=True, help="SSMD score file")
    parser.add_argument("-c", "--coords", required=True, help="Coordinate reference file")
    parser.add_argument("-n", "--sample", required=True, help="Sample (e.g. RWPE_R1)")
    parser.add_argument("-o", "--output", required=True, help="Output BED file")
    parser.add_argument(
        "-m", "--mode",
        choices=["sgrna", "ocr"],
        required=True,
        help="Output mode: sgrna or ocr"
    )

    return parser.parse_args()


def parse_coords(df, coord_col):
    df = df.dropna(subset=[coord_col]).copy()

    coords = df[coord_col].str.extract(r"(chr[^:]+):(\d+)-(\d+)")
    df["chrom"] = coords[0]
    df["chromStart"] = coords[1].astype(int)
    df["chromEnd"] = coords[2].astype(int)
    
    return df


def make_sgrna_bed(df_scores, df_coords, sample):

    # Merge on sgRNA
    df = df_scores.merge(df_coords, left_on="guide.set", right_on="sgRNA", how="inner")

    df = parse_coords(df, f"Coord_{sample.split('_')[0]}")

    score_col = f"{sample}_D14.Score"

    bed = df[[
        "chrom",
        "chromStart",
        "chromEnd",
        "guide.set",
        score_col
    ]].copy()

    bed.columns = ["chrom", "chromStart", "chromEnd", "name", "score"]

    return bed


def make_ocr_bed(df_scores, df_coords, sample):

    # Merge first
    df = df_scores.merge(df_coords, left_on="guide.set", right_on="sgRNA", how="inner")
    
    score_col = f"{sample}_D14.Score"

    # Extract OCR (R1234 from R1234_1)
    df["open_reg"] = df["guide.set"].str.rsplit("_", n=1).str[0]

    # Aggregate (choose strongest effect → min score)
    df_ocr = df.groupby("open_reg")[score_col].min().reset_index()

    # Get unique coordinates per OCR
    coord_map = df_coords.drop_duplicates("open_reg")[["open_reg", f"Coord_{sample.split('_')[0]}"]]

    df_ocr = df_ocr.merge(coord_map, on="open_reg", how="left")

    df_ocr = parse_coords(df_ocr, f"Coord_{sample.split('_')[0]}")

    bed = df_ocr[[
        "chrom",
        "chromStart",
        "chromEnd",
        "open_reg",
        score_col
    ]].copy()

    bed.columns = ["chrom", "chromStart", "chromEnd", "name", "score"]

    return bed


def main():
    args = get_args()

    df_scores = pd.read_csv(args.results, sep="\t")
    df_coords = pd.read_csv(args.coords)

    if args.mode == "sgrna":
        bed = make_sgrna_bed(df_scores, df_coords, args.sample)

    elif args.mode == "ocr":
        bed = make_ocr_bed(df_scores, df_coords, args.sample)

    else:
        raise ValueError("Invalid mode")

    bed.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()