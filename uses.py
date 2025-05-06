import requests
import pandas as pd


list_file = "data/cpdat_v4.0/cpdat_v4.0_list_presence_data.csv"
functional_file = "data/cpdat_v4.0/cpdat_v4.0_functional_use_data.csv"
output_file = "data/cpdat_merged_files.csv"
# Load both files
list_presence_df = pd.read_csv(list_file, low_memory=False)
functional_df = pd.read_csv(functional_file, low_memory=False)

# Select relevant columns
functional_cols = functional_df[[
    'data_source', 'data_group', 'curated_chemical_name', 'curated_casrn', 'dtxsid',
    'raw_functional_use', 'function_category', 'component_name'
]]

list_cols = list_presence_df[[
    'data_source', 'curated_chemical_name', 'curated_casrn', 'dtxsid',
    'data_group', 'organization', 'component_name', 'keyword_set'
]]
# Rename to avoid confusion between data_groups
list_df = pd.read_csv(list_file, usecols=list_cols, low_memory=False)
list_df = list_df.rename(columns={
    'data_group': 'list_name',
    'organization': 'list_organization',
    'data_source': 'list_data_source'
})

# Prepare output file
first_chunk = True
chunk_size = 5_000
reader = pd.read_csv(functional_file, usecols=functional_cols, chunksize=chunk_size, low_memory=False)

# Process and merge in chunks
for i, func_chunk in enumerate(reader):
    print(f"Processing chunk {i + 1}...")

    merged_chunk = pd.merge(func_chunk, list_df, on='dtxsid', how='left')

    # Save the merged chunk
    merged_chunk.to_csv(output_file, mode='a', index=False, header=first_chunk)
    first_chunk = False

print(f"✅ Merged file saved at: {output_file}")
