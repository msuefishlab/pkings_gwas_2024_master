import pandas as pd
import os

# Function to replace S3 path but keep the filename
def replace_s3_path(s3_path):
    if pd.isna(s3_path):
        return s3_path
    bucket_name = os.environ.get('bucket')
    print(bucket_name)
    return f"gs://{bucket_name}/reads/" + s3_path.split('/')[-1]

# Load the CSV file
df = pd.read_csv('./data_model/Deduplicated Inventory-GWAS_Set-2.csv')

# Load the TSV file
fish_data = pd.read_csv('./data_model/fish_data_2024.txt', sep='\t')

# Filter out rows where 'file_readno' or 'S3 Path' is missing
df = df.dropna(subset=['file_readno', 'S3 Path'])

# Create a pivot table
pivot_df = df.pivot_table(index=['Individual_ID', 'barcode', 'lane', 'fcid'],
                          columns='file_readno', 
                          values='S3 Path', 
                          aggfunc='first').reset_index()

# Replace S3 paths in the R1 and R2 columns and rename them
pivot_df['S3 Path R1'] = pivot_df['R1'].apply(replace_s3_path)
pivot_df['S3 Path R2'] = pivot_df['R2'].apply(replace_s3_path)

# Drop the original R1 and R2 columns
pivot_df.drop(columns=['R1', 'R2'], inplace=True)

# Merge with fish_data to get the PN column
merged_df = pivot_df.merge(fish_data, left_on='Individual_ID', right_on='SpecNo')

# Check for duplicate 'participant' columns and remove if necessary
if 'participant' in merged_df.columns:
    # Select the first 'participant' column
    merged_df['participant'] = merged_df.loc[:, 'participant'].iloc[:, 0]

# Add new columns
merged_df['entity:sample_id'] = merged_df['fcid'] + "_" + merged_df['lane'].astype(str) + "_" + merged_df['PN']
merged_df['ID'] = merged_df['fcid'] + "_" + merged_df['lane'].astype(str)
merged_df['PL'] = "Illumina"
merged_df['PU'] = merged_df['fcid'] + "_" + merged_df['lane'].astype(str) + "_" + merged_df['barcode']

# Rename columns
merged_df.rename(columns={'fcid': 'LB', 'PN': 'participant', 'S3 Path R1':'R1',"S3 Path R2": 'R2'}, inplace=True)

# Assign 'participant' to 'SN'
merged_df['SN'] = merged_df['participant']

# Reorder columns
final_columns = ['entity:sample_id', 'ID', 'Individual_ID', 'LB', 'lane', 'PL', 'PU', 'R1', 'R2', 'barcode', 'participant', 'SN']
final_df = merged_df[final_columns]

# Save the result to a new TSV file
final_df.to_csv('./data_model/samples_2024.tsv', sep='\t', index=False)
