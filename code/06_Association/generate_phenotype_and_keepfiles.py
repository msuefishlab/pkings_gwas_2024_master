import pandas as pd
import argparse

def list_unique_values(file_path):
    # Read the file
    data = pd.read_csv(file_path, sep='\t')

    # Get unique values
    unique_populations = data['POP'].dropna().unique().astype(str)
    unique_phenotypes = data['Phenotype'].dropna().unique().astype(str)
    unique_species = data['Species'].dropna().unique().astype(str)


    # Print comma-separated lists
    print("Unique Populations:", ",".join(unique_populations))
    print("Unique Phenotypes:", ",".join(unique_phenotypes))
    print("Unique Species:", ",".join(unique_species))

def filter_data(file_path, output_directory, include_pop=None, exclude_pop=None, include_phenotype=None, exclude_phenotype=None, include_species=None, exclude_species=None):
    # Read the file
    data = pd.read_csv(file_path, sep='\t')

    # Apply filters
    if include_pop:
        data = data[data['POP'].isin(include_pop)]
    if exclude_pop:
        data = data[~data['POP'].isin(exclude_pop)]
    if include_phenotype:
        data = data[data['Phenotype'].isin(include_phenotype)]
    if exclude_phenotype:
        data = data[~data['Phenotype'].isin(exclude_phenotype)]
    if include_species:
        data = data[data['Species'].isin(include_species)]
    if exclude_species:
        data = data[~data['Species'].isin(exclude_species)]

    # Get unique phenotypes
    unique_phenotypes = data['Phenotype'].dropna().unique()

    # Dictionary to store user input
    phenotype_values = {}

    # Ask user for input for each unique phenotype
    for phenotype in unique_phenotypes:
        while True:
            try:
                user_input = int(input(f"Select a value (0, 1, 2, or -9) for phenotype '{phenotype}': "))
                if user_input in [0, 1, 2, -9]:
                    phenotype_values[phenotype] = int(user_input)
                    break
                else:
                    print("Invalid input. Please enter 0, 1, 2, or -9.")
            except ValueError:
                print("Invalid input. Please enter a number.")

    # Map the user input values to the data
    data['Mapped_Value'] = data['Phenotype'].map(phenotype_values).fillna(-9).astype(int)

   # Prompt for output filename
    outfilename = input("Enter the name for the output files (without extension): ")

    # Prepare data for .keep.txt file
    keep_data = data[['entity:participant_id']].copy()
    keep_file_path = f"{output_directory}/{outfilename}.keep.txt"
    keep_data.to_csv(keep_file_path, sep='\t', index=False, header=False)

    # Prepare data for .pheno.txt file
    pheno_data = data[['entity:participant_id','entity:participant_id', 'Mapped_Value']].copy()
    pheno_data.columns = ['FID', 'IID', 'pheno']
    pheno_file_path = f"{output_directory}/{outfilename}.pheno.txt"
    pheno_data.to_csv(pheno_file_path, sep='\t', index=False, header=True)

    print(f"Files saved: {keep_file_path} and {pheno_file_path}")

def parse_arguments():
    parser = argparse.ArgumentParser(description='Filter data or list unique values based on POP, Phenotype, and Species.')
    parser.add_argument('file_path', type=str, help='Path to the data file')
    parser.add_argument('--list', action='store_true', help='List unique values of POP, Phenotype, and Species')
    parser.add_argument('--include_pop', nargs='*', help='POP values to include')
    parser.add_argument('--exclude_pop', nargs='*', help='POP values to exclude')
    parser.add_argument('--include_phenotype', nargs='*', help='Phenotype values to include')
    parser.add_argument('--exclude_phenotype', nargs='*', help='Phenotype values to exclude')
    parser.add_argument('--include_species', nargs='*', help='Species to include')
    parser.add_argument('--exclude_species', nargs='*', help='Species to exclude')
    parser.add_argument('--output_directory', type=str, help='Directory to save output files')
    return parser.parse_args()

def main():
    args = parse_arguments()

    if args.list:
        list_unique_values(args.file_path)
    else:
        filtered_data = filter_data(args.file_path, args.output_directory,
                                    include_pop=args.include_pop, 
                                    exclude_pop=args.exclude_pop,
                                    include_phenotype=args.include_phenotype, 
                                    exclude_phenotype=args.exclude_phenotype,
                                    include_species=args.include_species,
                                    exclude_species=args.exclude_species)

if __name__ == "__main__":
    main()
