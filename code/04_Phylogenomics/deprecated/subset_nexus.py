import argparse
from Bio.Nexus import Nexus

def subset_nexus(nexus_file, individuals_file, output_file):
    # Read the list of individuals from the text file
    with open(individuals_file, 'r') as f:
        individuals_to_include = set(line.strip() for line in f)

    # Parse the Nexus file
    nexus_data = Nexus.Nexus()
    nexus_data.read(nexus_file)

    # Filter the taxa to only include the ones in the individuals list
    taxa_to_keep = [taxon for taxon in nexus_data.taxlabels if taxon in individuals_to_include]

    # Create a new Nexus object for the subset
    subset_nexus = Nexus.Nexus()

    # Subset the matrix by filtering taxa
    for taxon in taxa_to_keep:
        subset_nexus.matrix[taxon] = nexus_data.matrix[taxon]

    # Copy the taxlabels to the subset Nexus
    subset_nexus.taxlabels = taxa_to_keep

    # Write the subset Nexus file
    with open(output_file, 'w') as outfile:
        subset_nexus.write_nexus_data(outfile)

    print(f"Subset Nexus file saved to {output_file}")

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Subset a Nexus file to include only specific individuals.")
    parser.add_argument("nexus_file", help="Path to the Nexus file")
    parser.add_argument("individuals_file", help="Path to the text file with the list of individuals to include")
    parser.add_argument("output_file", help="Path to save the subset Nexus file")

    # Parse arguments
    args = parser.parse_args()

    # Call the subset_nexus function with the arguments
    subset_nexus(args.nexus_file, args.individuals_file, args.output_file)

if __name__ == "__main__":
    main()
