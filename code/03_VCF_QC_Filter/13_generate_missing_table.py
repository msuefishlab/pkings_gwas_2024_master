#!/usr/bin/env python3

import sys
import os
import glob

# Replace 'genotypes.txt' with your actual genotype filename
input_file = 'output_data/03_QC/MIC4273_HAP2_NEWREF_MSU618_all_fish_all_sites.FILTERED_BIALLELIC_SNPS_ONLY.renamed.filtered.stats.tab'
output_file = 'output_data/03_QC/MIC4273_HAP2_NEWREF_MSU618_all_fish_all_sites.missingness_per_site_per_group.txt'

# Directory where your keep files are located
keep_files_directory = 'input_data/06_Association/'  # Current directory

# Pattern to match keep files
keep_files_pattern = '*.keep.txt'

# Function to read keep files and build group mapping
def build_group_mapping(keep_files):
    group_samples = {}
    for keep_file in keep_files:
        group_name = os.path.splitext(os.path.basename(keep_file))[0]
        # Optionally, remove '.keep' from group name if present
        group_name = group_name.replace('.keep', '')
        with open(keep_file, 'r') as f:
            sample_ids = [line.strip() for line in f if line.strip()]
            group_samples[group_name] = set(sample_ids)
    return group_samples

def main():
    # Use glob to find all keep files matching the pattern
    keep_files = glob.glob(os.path.join(keep_files_directory, keep_files_pattern))
    if not keep_files:
        print(f"No keep files found matching pattern '{keep_files_pattern}' in directory '{keep_files_directory}'.")
        sys.exit(1)

    group_samples = build_group_mapping(keep_files)

    # This will map sample IDs (with and without suffix) to their column indices
    sample_columns = {}
    sample_columns_no_suffix = {}

    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        header_line = f_in.readline().strip()
        headers = header_line.split('\t')

        # Build mapping from sample IDs to column indices
        for i, sample_header in enumerate(headers):
            sample_id = sample_header
            # Remove any suffix after the first '.'
            sample_id_no_suffix = sample_id.split('.')[0]
            sample_columns[sample_id] = i  # Column index of each sample
            sample_columns_no_suffix[sample_id_no_suffix] = i  # Map stripped sample ID to column index

        # Verify that all samples in keep files are present in the header
        for group_name, samples in group_samples.items():
            missing_samples = samples - set(sample_columns_no_suffix.keys())
            if missing_samples:
                print(f"Warning: The following samples in '{group_name}' are not found in the genotype file header: {', '.join(missing_samples)}")

        # Prepare group column indices
        group_column_indices = {}
        for group_name, samples in group_samples.items():
            indices = []
            for sample_id in samples:
                if sample_id in sample_columns_no_suffix:
                    idx = sample_columns_no_suffix[sample_id]
                    indices.append(idx)
                else:
                    # Sample ID not found, warning already printed
                    continue
            if not indices:
                print(f"Warning: No valid samples found for group '{group_name}'.")
            group_column_indices[group_name] = indices

        # Write header for output file
        f_out.write('CHROM\tPOS\t' + '\t'.join(sorted(group_column_indices.keys())) + '\n')

        # Now, process each line
        for line_number, line in enumerate(f_in, start=2):
            line = line.strip()
            if not line:
                continue  # Skip empty lines
            fields = line.split('\t')

            # Ensure the line has the correct number of fields
            if len(fields) != len(headers):
                print(f"Warning: Line {line_number} has {len(fields)} fields, expected {len(headers)}. Skipping line.")
                continue

            # For each group, calculate missingness
            missingness_per_group = {}
            for group_name, indices in group_column_indices.items():
                total_samples = len(indices)
                if total_samples == 0:
                    missingness_per_group[group_name] = 'NA'
                    continue
                missing_samples = sum(1 for idx in indices if fields[idx] == './.')
                missingness = missing_samples / total_samples
                missingness_per_group[group_name] = missingness

            # Write the results to the output file
            chrom = fields[0]
            pos = fields[1]
            missingness_values = [
                str(missingness_per_group[group]) if missingness_per_group[group] != 'NA' else 'NA'
                for group in sorted(group_column_indices.keys())
            ]
            f_out.write(f"{chrom}\t{pos}\t" + '\t'.join(missingness_values) + '\n')

if __name__ == "__main__":
    main()
