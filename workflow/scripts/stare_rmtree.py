import os
import shutil
import argparse

parser = argparse.ArgumentParser(description="Check if res folder exists")
parser.add_argument('--metadata', help="Path to metadata file")
parser.add_argument('--outdir', help='Path to stare results dir')
output = parser.parse_args()

# Check if metadata file is missing but directoyr exists:
if not os.path.exists(output.metadata) and os.path.exists(output.outdir):
    print(f"Metadata file missing, removing old output directory: {output.outdir}")
    shutil.rmtree(output.outdir)


