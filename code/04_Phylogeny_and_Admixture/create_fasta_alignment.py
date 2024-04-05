import os
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import vcf
from collections import defaultdict
import argparse  # Import argparse

# Mapping from nucleotide pairs to IUPAC codes
iupac_codes = {
    ('A', 'C'): 'M', ('A', 'G'): 'R', ('A', 'T'): 'W',
    ('C', 'G'): 'S', ('C', 'T'): 'Y', ('G', 'T'): 'K',
    ('A', 'A'): 'A', ('C', 'C'): 'C', ('G', 'G'): 'G', ('T', 'T'): 'T',
    # Include reverse combinations
    ('C', 'A'): 'M', ('G', 'A'): 'R', ('T', 'A'): 'W',
    ('G', 'C'): 'S', ('T', 'C'): 'Y', ('T', 'G'): 'K'
}

def get_iupac_code(alleles):
    """Return the IUPAC code for a pair of alleles."""
    # Sort alleles to match the keys in the iupac_codes dictionary
    sorted_alleles = tuple(sorted(alleles))
    return iupac_codes.get(sorted_alleles, 'N')  # Default to 'N' if not found

def read_ranges(file_path):
    """Read ranges from a file and return a list of tuples (gene_name, chromosome, start, end)."""
    ranges = []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('@')
            gene_name = parts[0]
            chrom, positions = parts[1].split(':')
            start, end = map(int, positions.split('-'))
            ranges.append((gene_name, chrom, start, end))
    return ranges

def fasta_alignment_from_vcf_for_ranges(vcf_file, ranges_file, output_dir):
    """Generate FASTA files for specified ranges from a VCF file."""
    ranges = read_ranges(ranges_file)
    vcf_reader = vcf.Reader(open(vcf_file, 'rb'))

    for gene_name, chrom, start, end in ranges:
        result = defaultdict(list)
        sites = []

        for record in vcf_reader:
            if record.CHROM == chrom and start <= record.POS <= end:
                ref = record.REF
                result['ref'].append(ref)
                sites.append(record.POS)
                for sample in record.samples:
                    name = sample.sample
                    if sample.gt_bases:
                        alleles = sample.gt_bases.replace('|', '/').split('/')
                        if len(set(alleles)) == 1:
                            result[name].append(alleles[0])
                        else:
                            iupac_code = get_iupac_code(alleles)
                            result[name].append(iupac_code)
                    else:
                        result[name].append('N')

        recs = [SeqRecord(Seq(''.join(result[sample])), id=sample) for sample in result]

        fasta_output = os.path.join(output_dir, f"{gene_name}.fasta")
        with open(fasta_output, 'w') as output_handle:
            SeqIO.write(recs, output_handle, "fasta")

# Setup argparse
parser = argparse.ArgumentParser(description='Generate FASTA files from VCF for specified ranges.')
parser.add_argument('vcf_file', type=str, help='Path to the VCF file')
parser.add_argument('ranges_file', type=str, help='Path to the file containing ranges')
parser.add_argument('output_dir', type=str, help='Directory to save the output FASTA files')

# Parse arguments
args = parser.parse_args()

# Ensure output directory exists
if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

# Example usage: python script.py path/to/your/vcf_file.vcf path/to/your/ranges.txt path/to/output_dir
fasta_alignment_from_vcf_for_ranges(args.vcf_file, args.ranges_file, args.output_dir)
