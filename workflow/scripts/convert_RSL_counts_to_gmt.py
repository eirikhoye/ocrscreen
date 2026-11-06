import pandas as pd
import argparse

def annotate_NC(coord_file, peak_id_col='Peak ID', coord_col='Coord_RWPE', control_ID="NCR"):
    df = pd.read_csv(coord_file, sep='\t')#.copy()

    is_ncr = df[peak_id_col].str.startswith(f"{control_ID}", na=False)
    ncr_coords = df.loc[is_ncr, coord_col]
    unique_coords = ncr_coords.unique()
    coord_to_ncr_num = {coord: i + 1 for i, coord in enumerate(unique_coords)}

    df.loc[is_ncr, peak_id_col] = (
        df.loc[is_ncr, coord_col]
        .map(coord_to_ncr_num)
        .apply(lambda n: f"{control_ID}{n}")
    )
    return df



def create_gmt_file(input_file, coord_file, coord_col, control_ID, output_file, peak_id_col='Peak ID'):
    df = pd.read_csv(input_file, sep='\t')
    coords = annotate_NC(coord_file, peak_id_col=peak_id_col, coord_col=coord_col, control_ID=control_ID)

    lookup = coords.set_index("Guide ID")["Peak ID"]
    df["open_chromatin_region"] = df["guide.set"].map(lookup)

    guide_sets = {}

    for _, row in df.iterrows():
        region_id = row['open_chromatin_region'].split('_')[0]
        guide_rna = row['guide.set']
        
        if region_id not in guide_sets:
            guide_sets[region_id] = []
        
        guide_sets[region_id].append(guide_rna)
    
    with open(output_file, 'w') as f:
        for region_id, guides in guide_sets.items():
            # Create the GMT line for each region
            gmt_line = '\t'.join([region_id, 'open_chrom'] + guides)
            f.write(gmt_line + '\n')

    print(f"GMT file has been written to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Convert gene-summary CSV file to GMT format")
    parser.add_argument('input_file', help="Path to the input gene-summary CSV file")
    parser.add_argument('coord_file', help="Path to the guide RNA coordinates file")
    parser.add_argument('coord_col', help="Column name that contains coordinates")
    parser.add_argument('control_ID', help="Name of control")
    parser.add_argument('output_file', help="Path to the output GMT file")

    args = parser.parse_args()

    create_gmt_file(args.input_file, 
                    args.coord_file, 
                    f"Coord_{args.coord_col}", 
                    args.control_ID, 
                    args.output_file, 
                    peak_id_col='Peak ID')

if __name__ == '__main__':
    main()
