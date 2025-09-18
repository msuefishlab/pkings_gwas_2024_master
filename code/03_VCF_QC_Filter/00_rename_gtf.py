import pandas as pd
import argparse
import csv

# Setup command line argument parsing
parser = argparse.ArgumentParser(description='Update chromosome names in a GTF file based on a mapping file.')
parser.add_argument('gtf_file', type=str, help='Path to the GTF file')
parser.add_argument('map_file', type=str, help='Path to the mapping file')
parser.add_argument('output_file', type=str, help='Path to save the updated GTF file')
args = parser.parse_args()

# Load the mapping file
map_df = pd.read_csv(args.map_file, sep=',', header=None, names=['old_name', 'new_name'])

# Create a dictionary from the mapping file for quick look-up
name_map = dict(zip(map_df['old_name'], map_df['new_name']))

# Load the GTF file
gtf_df = pd.read_csv(args.gtf_file, sep='\t', header=None, comment='#',
                     names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])

# Replace old chromosome names using the map
gtf_df['seqname'] = gtf_df['seqname'].map(name_map).fillna(gtf_df['seqname'])

# Save the updated GTF file, ensuring no extra quotes are added
with open(args.output_file, 'w', newline='') as out_file:
    gtf_df.to_csv(out_file, sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE, escapechar='\\')
