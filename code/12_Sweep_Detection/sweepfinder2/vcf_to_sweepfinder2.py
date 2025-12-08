#!/usr/bin/env python3
"""
Convert polarized VCF to SweepFinder2 input format

This script converts a polarized VCF file (with INFO/AA ancestral allele annotations)
to the SweepFinder2 SFS input format.

SweepFinder2 format:
position    x    n    folded
12345       5    20   0

Where:
- position: genomic position (bp)
- x: count of derived allele
- n: total allele number (2 * number of individuals)
- folded: 0 for unfolded/polarized SFS, 1 for folded

Ancestral allele logic:
- If AA=REF: derived allele = ALT, x = AC
- If AA=ALT: derived allele = REF, x = AN - AC
- If AA=. or AA=N: skip (cannot polarize)

Usage:
    python3 vcf_to_sweepfinder2.py input.vcf.gz output.sf2.txt [--chromosome chr6]

Author: Jason Gallant Lab
Date: 2024
"""

import sys
import subprocess
import argparse


def vcf_to_sf2(vcf_path, output_path, chromosome=None):
    """
    Convert VCF with INFO/AA to SweepFinder2 format

    Args:
        vcf_path: Path to polarized VCF file (can be gzipped)
        output_path: Output file path
        chromosome: Optional chromosome filter (e.g., 'chr6')
    """

    # Build bcftools query command
    # Extract: CHROM, POS, REF, ALT, INFO/AA, INFO/AC, INFO/AN
    query_cmd = [
        'bcftools', 'query',
        '-f', '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA\t%INFO/AC\t%INFO/AN\n',
        vcf_path
    ]

    # Add chromosome filter if specified
    if chromosome:
        query_cmd.extend(['-r', chromosome])

    # Open output file
    with open(output_path, 'w') as out:
        # Write SweepFinder2 header
        out.write("position\tx\tn\tfolded\n")

        # Run bcftools query
        proc = subprocess.Popen(query_cmd, stdout=subprocess.PIPE, text=True)

        processed = 0
        skipped_no_aa = 0
        skipped_missing = 0
        skipped_aa_mismatch = 0

        for line in proc.stdout:
            fields = line.strip().split('\t')

            # Handle cases where fields might be missing
            if len(fields) < 7:
                skipped_missing += 1
                continue

            chrom, pos, ref, alt, aa, ac, an = fields

            # Skip if not polarizable (AA is missing or N)
            if aa == '.' or aa == 'N' or aa == '':
                skipped_no_aa += 1
                continue

            # Skip if missing allele count data
            if an == '.' or an == '0' or an == '' or ac == '.' or ac == '':
                skipped_missing += 1
                continue

            # Convert to integers
            try:
                ac = int(ac)
                an = int(an)
            except ValueError:
                skipped_missing += 1
                continue

            # Determine derived allele count based on ancestral allele
            if aa == ref:
                # Ancestral = REF, derived = ALT
                x = ac
            elif aa == alt:
                # Ancestral = ALT, derived = REF
                x = an - ac
            else:
                # AA doesn't match REF or ALT (could be multiallelic or error)
                # Skip these sites
                skipped_aa_mismatch += 1
                continue

            # Write SweepFinder2 format: position, x, n, folded
            # folded=0 indicates unfolded/polarized SFS
            out.write(f"{pos}\t{x}\t{an}\t0\n")
            processed += 1

        # Wait for bcftools to finish
        proc.wait()

    # Print summary statistics
    print(f"Conversion complete: {output_path}")
    print(f"  Processed sites: {processed}")
    print(f"  Skipped (no AA): {skipped_no_aa}")
    print(f"  Skipped (missing data): {skipped_missing}")
    print(f"  Skipped (AA mismatch): {skipped_aa_mismatch}")
    print(f"  Total attempted: {processed + skipped_no_aa + skipped_missing + skipped_aa_mismatch}")


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description='Convert polarized VCF to SweepFinder2 format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('vcf', help='Input polarized VCF (can be gzipped)')
    parser.add_argument('output', help='Output SweepFinder2 format file')
    parser.add_argument('-c', '--chromosome',
                       help='Restrict to chromosome (e.g., chr6)')

    args = parser.parse_args()

    # Run conversion
    vcf_to_sf2(args.vcf, args.output, args.chromosome)


if __name__ == "__main__":
    main()
