import random

def read_nexus_file(file_path):
    """Reads the nexus file and extracts the charset and partition definitions."""
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # Extract header (everything before the `begin sets;` block)
    header = []
    charset_block = []
    partition_block = []
    footer = []
    
    inside_charset = False
    inside_partition = False
    inside_footer = False
    
    for line in lines:
        line = line.strip()
        
        if line.lower() == 'begin sets;':
            inside_charset = True
            charset_block.append(line + "\n")
            continue
        elif line.lower() == 'charpartition mymodels =':
            inside_charset = False
            inside_partition = True
            partition_block.append(line + "\n")
            continue
        elif inside_partition and line.lower() == 'end;':
            inside_partition = False
            inside_footer = True
            footer.append(line + "\n")
            continue

        if inside_charset:
            charset_block.append(line + "\n")
        elif inside_partition:
            partition_block.append(line + "\n")
        elif inside_footer:
            footer.append(line + "\n")
        else:
            header.append(line + "\n")
    
    return header, charset_block, partition_block, footer

def write_nexus_file(output_file, header, charset_block, partition_block, footer):
    """Writes the new Nexus file with the subset of charsets and partitions."""
    with open(output_file, 'w') as file:
        file.writelines(header)
        file.write("begin sets;\n")
        file.writelines(charset_block)
        file.write("charpartition mymodels =\n")
        file.writelines(partition_block)
        file.writelines(footer)

def parse_charset_block(charset_block):
    """Parses the charset block into a dictionary for easy access."""
    charset_dict = {}
    for line in charset_block:
        if 'charset' in line:
            parts = line.split()
            fasta_file = parts[1].strip()  # Get the fasta file name (e.g., LOC111832677.fasta)
            charset_dict[fasta_file] = line
    return charset_dict

def parse_partition_block(partition_block):
    """Parses the partition block into a list of tuples (fasta_file, partition_model)."""
    partitions = []
    for line in partition_block:
        if ":" in line:  # Find the lines that have models and file associations
            fasta_file = line.split(":")[1].split("{")[0].strip()  # Get the fasta file name
            partitions.append((fasta_file, line))
    return partitions

def randomly_subset_partitions(partitions, charset_dict, n):
    """Randomly selects n partitions and returns the corresponding charset and partitions."""
    if n >= len(partitions):
        return [charset_dict[part[0]] for part in partitions], [part[1] for part in partitions]
    
    subset_partitions = random.sample(partitions, n)
    subset_charset = [charset_dict[part[0]] for part in subset_partitions]
    
    # Ensure the last partition ends with a semicolon instead of a comma
    subset_partition = [part[1].rstrip(',') + ';' if i == len(subset_partitions) - 1 else part[1]
                        for i, part in enumerate(subset_partitions)]
    
    return subset_charset, subset_partition

# Set your input and output file paths
input_file = 'loci.best_model.repaired.nex'  # Replace with the path to your Nexus file
output_file = 'loci.best_model.repaired.subset.nex'  # Replace with your desired output file name

# Number of partitions to subset
n = 200  # Replace this with the number of partitions you want

# Read the Nexus file
header, charset_block, partition_block, footer = read_nexus_file(input_file)

# Parse the charset and partition blocks
charset_dict = parse_charset_block(charset_block)
partitions = parse_partition_block(partition_block)

# Randomly subset n partitions and their corresponding charset entries
subset_charset, subset_partition = randomly_subset_partitions(partitions, charset_dict, n)

# Write the new Nexus file with the subset
write_nexus_file(output_file, header, subset_charset, subset_partition, footer)

print(f"Successfully wrote {n} random partitions to {output_file}")
