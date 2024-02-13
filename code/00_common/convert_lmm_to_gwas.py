import pandas as pd
import sys

def convert_lmm_to_gwas(lmm_file, gwas_file):
    # Read the LMM file
    lmm_data = pd.read_csv(lmm_file, sep='\t')

    # Rename columns
    lmm_data.rename(columns={'rs': 'SNP', 'ps': 'BP', 'p_lrt': 'P'}, inplace=True)

    # Filter out rows where 'P' is NaN or blank
    lmm_data = lmm_data[pd.notna(lmm_data['P'])]

    # Select only the required columns
    required_columns = ['chr', 'BP', 'SNP', 'P']
    gwas_data = lmm_data[required_columns]

    # Rename 'chr' to 'CHR'
    gwas_data.rename(columns={'chr': 'CHR'}, inplace=True)

    # Save to GWAS file
    gwas_data.to_csv(gwas_file, sep='\t', index=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_lmm_file> <output_gwas_file>")
        sys.exit(1)

    input_lmm_file = sys.argv[1]
    output_gwas_file = sys.argv[2]

    convert_lmm_to_gwas(input_lmm_file, output_gwas_file)
    print("Conversion complete!")
