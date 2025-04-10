import pandas as pd
import gzip
import random
import argparse
import os

def extract_sequence_regions(df, read_length):
    #df = df.copy()
    df["read_upstream"] = df.apply(lambda row: row["full_molecule_sequence"][row["read_start"]:row["read_start"] + read_length], axis=1)
    df["read_downstream"] = df.apply(lambda row: row["full_molecule_sequence"][max(0, row["read_end"] - read_length):row["read_end"]], axis=1)
    return df

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(complement)[::-1]

def process_sequence(sequence, reverse_comp=False):
    """Process sequence by reversing and complementing if needed."""
    return reverse_complement(sequence) if reverse_comp else sequence

def convert_to_fastq(df, output_prefix, sequencing_type, strandedness):
    """Convert dataframe to FASTQ format based on sequencing type and strandedness."""
    if sequencing_type == "SE":
        output_file = f"{output_prefix}.fastq.gz"
        with gzip.open(output_file, 'wt') as f:
            for _, row in df.iterrows():
                read_name = f"{row['molecule_id']}"
                
                if strandedness == "rf":
                    sequence = process_sequence(row['read_downstream'], reverse_comp=True)
                elif strandedness == "fr":
                    sequence = process_sequence(row['read_upstream'], reverse_comp=False)
                elif strandedness == "unstranded":
                    if random.choice([True, False]):
                        sequence = process_sequence(row['read_downstream'], reverse_comp=True)
                    else:
                        sequence = process_sequence(row['read_upstream'], reverse_comp=False)
                
                quality_scores = "I" * len(sequence)  # Dummy quality scores
                
                f.write(f"@{read_name}\n")
                f.write(f"{sequence}\n")
                f.write(f"+\n")
                f.write(f"{quality_scores}\n")
    
    elif sequencing_type == "PE":
        output_file_r1 = f"{output_prefix}_R1.fastq.gz"
        output_file_r2 = f"{output_prefix}_R2.fastq.gz"
        
        with gzip.open(output_file_r1, 'wt') as f1, gzip.open(output_file_r2, 'wt') as f2:
            for _, row in df.iterrows():
                read_name = f"{row['molecule_id']}"
                
                if strandedness == "rf":
                    sequence_r1 = process_sequence(row['read_downstream'], reverse_comp=True)
                    sequence_r2 = process_sequence(row['read_upstream'], reverse_comp=False)
                elif strandedness == "fr":
                    sequence_r1 = process_sequence(row['read_upstream'], reverse_comp=False)
                    sequence_r2 = process_sequence(row['read_downstream'], reverse_comp=True)
                elif strandedness == "unstranded":
                    if random.choice([True, False]):
                        sequence_r1 = process_sequence(row['read_downstream'], reverse_comp=True)
                        sequence_r2 = process_sequence(row['read_upstream'], reverse_comp=False)
                    else:
                        sequence_r1 = process_sequence(row['read_upstream'], reverse_comp=False)
                        sequence_r2 = process_sequence(row['read_downstream'], reverse_comp=True)
                
                quality_scores_r1 = "I" * len(sequence_r1)
                quality_scores_r2 = "I" * len(sequence_r2)
                
                f1.write(f"@{read_name}\n")
                f1.write(f"{sequence_r1}\n")
                f1.write(f"+\n")
                f1.write(f"{quality_scores_r1}\n")
                
                f2.write(f"@{read_name}\n")
                f2.write(f"{sequence_r2}\n")
                f2.write(f"+\n")
                f2.write(f"{quality_scores_r2}\n")




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="convert dataframe to FASTQ file.")
    parser.add_argument("--input_df", help="full path to the input dataframe")
    parser.add_argument("--read_length", type=int,default=100,help="length of each read")
    parser.add_argument("--seq_type", choices=["PE", "SE"], type=str, default='SE', help="Type of sequencing: SE or PE")
    parser.add_argument("--s", choices=["rf", "fr","unstranded"], type=str, default='rf', help="strandedness of the simulated library")
    parser.add_argument('--o', type=str, default='./', help='output path')
    args = parser.parse_args()


# Generate output name and file for output file with read and substitution info
    filename = os.path.splitext(os.path.basename(args.input_df))[0]
    # Extract only the gene ID (first part before the first underscore)
    base_filename = filename.split("_")[0]
    # Construct the output filename
    output_prefix = str(args.o) + base_filename
    print(output_prefix)
    df = pd.read_csv(args.input_df, delimiter="\t")

    result_df = extract_sequence_regions(df, read_length=100)

    convert_to_fastq(result_df, output_prefix, sequencing_type=args.seq_type, strandedness=args.s)





