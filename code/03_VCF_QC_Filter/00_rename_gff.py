import pandas as pd
import argparse

# Setup command line argument parsing
parser = argparse.ArgumentParser(description='Update chromosome names in a GFF file based on a mapping file.')
parser.add_argument('gff_file', type=str, help='Path to the GFF file')
parser.add_argument('map_file', type=str, help='Path to the mapping file')
parser.add_argument('output_file', type=str, help='Path to save the updated GFF file')
args = parser.parse_args()

# Load the mapping file
map_df = pd.read_csv(args.map_file, sep=',', header=None, names=['old_name', 'new_name'])

# Create a dictionary from the mapping file for quick look-up
name_map = dict(zip(map_df['old_name'], map_df['new_name']))

# Load the GFF file
gff_df = pd.read_csv(args.gff_file, sep='\t', header=None, comment='#', names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])

# Replace old chromosome names using the map
gff_df['seqname'] = gff_df['seqname'].map(name_map).fillna(gff_df['seqname'])

# Save the updated GFF file
gff_df.to_csv(args.output_file, sep='\t', index=False, header=False)
