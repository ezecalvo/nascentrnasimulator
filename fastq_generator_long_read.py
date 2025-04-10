import pandas as pd
import random
import string
import argparse
import gzip
import os

def generate_random_string(length=8):
    """Generate a random string of given length."""
    return ''.join(random.choices(string.ascii_letters + string.digits, k=length))

def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""
    complement = str.maketrans('ACGTacgt', 'TGCAtgca')
    return sequence.translate(complement)[::-1]

def process_sequence(sequence, sequencing_type):
    """Process sequence based on sequencing type."""
    if sequencing_type == "RNA":
        return sequence[::-1]  # Reverse the sequence
    elif sequencing_type == "cDNA":
        return sequence if random.choice([True, False]) else reverse_complement(sequence)
    return sequence

def convert_to_fastq(df, output_file, sequencing_type):
    """Convert dataframe to FASTQ format and save to a gzipped file."""
    with gzip.open(output_file, 'wt') as f:  # 'wt' mode to write text to the gzip file
        for _, row in df.iterrows():
            read_name = read_name = f"{row['molecule_id']}"
            sequence = process_sequence(row['full_molecule_sequence'], sequencing_type)
            quality_scores = "I" * len(sequence)  # Dummy quality scores
            
            f.write(f"{read_name}\n")
            f.write(f"{sequence}\n")
            f.write(f"+\n")
            f.write(f"{quality_scores}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="convert dataframe to FASTQ file.")
    parser.add_argument("--input_df", help="full path to the input dataframe")
    parser.add_argument("--seq_type", choices=["RNA", "cDNA"], help="Type of sequencing: RNA or cDNA")
    parser.add_argument('--o', type=str, default='./', help='output path')
    args = parser.parse_args()
    

    # Generate output name and file for output file with read and substitution info
    filename = os.path.splitext(os.path.basename(args.input_df))[0]
    # Extract only the gene ID (first part before the first underscore)
    base_filename = filename.split("_")[0]
    # Construct the output filename
    output_filename = str(args.o) + base_filename +'_'+args.seq_type +".fastq.gz"

    df = pd.read_csv(args.input_df, delimiter="\t")
    
    convert_to_fastq(df, output_filename, args.seq_type)
