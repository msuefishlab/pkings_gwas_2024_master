import os
import sys

def merge_count_files(input_dir, output_file):
    # Get the list of count files
    count_files = [f for f in os.listdir(input_dir) if f.endswith('.counts')]
    count_file_bases = [os.path.splitext(f)[0] for f in count_files]

    # Create a dictionary to store file handles
    file_handles = {f: open(os.path.join(input_dir, f), 'r') for f in count_files}

    # Skip the first two lines (header) in each file
    for f in file_handles.values():
        next(f)  # Skip first header line
        next(f)  # Skip second header line

    # Open the output file
    with open(output_file, 'w') as f_out:
        # Write the header
        f_out.write(f"COUNTSFILE NPOP {len(count_files)} NSITES 77107170\n")
        f_out.write("CHROM\tPOS\t" + "\t".join(count_file_bases) + "\n")

        # Process each file line by line
        for _ in range(77107170):  # Assuming the total number of sites is known
            merged_line = {}
            
            # Read one line from each file and extract data
            for count_file, handle in file_handles.items():
                line = handle.readline().strip()
                chrom, pos, counts = line.split()[:3]
                merged_line[(chrom, pos)] = merged_line.get((chrom, pos), []) + [counts]

            # Write the merged data to the output file
            for (chrom, pos), counts in merged_line.items():
                f_out.write(f"{chrom}\t{pos}\t" + "\t".join(counts) + "\n")

    # Close all file handles
    for f in file_handles.values():
        f.close()

    print(f"Merged file saved as: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python merge_counts.py <input_dir> <output_file>")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_file = sys.argv[2]

    merge_count_files(input_dir, output_file)
