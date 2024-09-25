import os
import shutil
import argparse
from Bio import SeqIO

def check_file(file_path):
    """
    Check if a file passes the filtering criteria:
    - More than 10 individuals have more than 50% "Ns".
    - The file is not blank.
    Returns True if the file passes, False otherwise.
    """
    num_bad_alignments = 0
    total_records = 0

    try:
        for record in SeqIO.parse(file_path, "fasta"):
            total_records += 1
            sequence = str(record.seq)
            if sequence.count("N") / len(sequence) > 0.5:
                num_bad_alignments += 1
    except UnicodeDecodeError:
        print(f"Skipping file due to encoding error: {file_path}")
        return False

    if total_records == 0 or num_bad_alignments > 10:
        return False
    return True

def filter_files(source_dir, target_dir):
    """
    Scan each FASTA file in the source directory, apply filtering criteria,
    and copy the files that pass the filter to the target directory.
    """
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    for filename in os.listdir(source_dir):
        file_path = os.path.join(source_dir, filename)
        if os.path.isfile(file_path) and check_file(file_path):
            shutil.copy(file_path, os.path.join(target_dir, filename))

def main():
    parser = argparse.ArgumentParser(description="Filter FASTA files based on alignment criteria.")
    parser.add_argument("source_dir", type=str, help="The source directory containing FASTA files.")
    parser.add_argument("target_dir", type=str, help="The target directory for filtered FASTA files.")

    args = parser.parse_args()

    filter_files(args.source_dir, args.target_dir)

if __name__ == "__main__":
    main()

