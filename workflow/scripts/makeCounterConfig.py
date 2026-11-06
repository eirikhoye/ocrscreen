import os
import sys
import argparse

def parse_samples(samples_file, fastqdir):
    """Parse the samples file and return a list of sample information."""
    samples_and_files = []

    with open(samples_file, 'r') as infile:
        for line in infile:
            line = line.strip()
            if not line:
                continue
            fields = line.split('\t')
            if fields[0] != "file" and fields[0].endswith("_R1_001.fastq.gz"):
                fastqpref = fields[0].split("_R1_001.fastq.gz")[0]
                R1_read = os.path.join(fastqdir, fields[0])
                R2_read = os.path.join(fastqdir, f"{fastqpref}_R2_001.fastq.gz")
                sample = fields[1]
                samples_and_files.append((sample, R1_read, R2_read))
            else:
                print(f"Warning: Entry {fields[0]} in {samples_file} is not the expected format")

    return samples_and_files

def generate_config_file(template_file, outfile, library, out_cnttable, fastq_files, samples):
    """Generate the config file based on the template and input values."""
    with open(template_file, 'r') as infile, open(outfile, 'w') as out:
        for line in infile:
            if line.startswith('design ='):
                out.write(f"design = {library}\n")
            elif line.startswith('out ='):
                out.write(f"out = {out_cnttable}\n")
            elif line.startswith('fastq ='):
                out.write(f"fastq = {','.join(fastq_files['R1'])}\n")
            elif line.startswith('umifastq ='):
                out.write(f"umifastq = {','.join(fastq_files['R2'])}\n")
            elif line.startswith('samples ='):
                out.write(f"samples = {','.join(samples)}\n")
            else:
                out.write(line)

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Generate config file for CrisprCounter.jar")
    parser.add_argument('--template', required=True, help="Path to template.properties")
    parser.add_argument('--fastqdir', required=True, help="Path to the directory containing FASTQ files")
    parser.add_argument('--samples', required=True, help="Path to samples.txt")
    parser.add_argument('--library', required=True, help="Path to library CSV file")
    parser.add_argument('--prefix', required=True, help="Prefix for file names")
    parser.add_argument('--outdir', required=True, help="Output directory for the properties file")
    
    args = parser.parse_args()

    # File paths
    outfile = os.path.join(args.outdir, f"{args.prefix.split('/')[-1].split('.')[0]}.properties")
    out_cnttable = args.prefix

    # Parse samples and files
    samples_and_files = parse_samples(args.samples, args.fastqdir)
    sample_names = [sample[0] for sample in samples_and_files]
    R1_files = [sample[1] for sample in samples_and_files]
    R2_files = [sample[2] for sample in samples_and_files]

    # Dictionary to hold fastq files
    fastq_files = {'R1': R1_files, 'R2': R2_files}

    # Generate config file
    generate_config_file(args.template, outfile, args.library, out_cnttable, fastq_files, sample_names)

if __name__ == "__main__":
    main()
