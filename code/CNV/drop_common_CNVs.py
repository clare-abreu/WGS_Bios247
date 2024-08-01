#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 21:46:59 2024

@author: clare
"""
import pandas as pd
from collections import Counter
import os
import sys

# List of genome identifiers
genomes = ['WGS-March23-143-02-04', 'WGS-March23-143-03-33', 'WGS-March23-143-03-72', 'WGS-March23-143-04-24', 'WGS-March23-143-04-81']

# Dictionary to store the dataframes
bedgrphs = {}

# Load the data for each genome
for genome in genomes:
    file_path = f'{genome}.gt.bedgraph'
    if os.path.exists(file_path):
        bedgrphs[genome] = pd.read_csv(file_path, header=None, sep='\t', names=['chromosome', 'start', 'end', 'coverage'])
    else:
        print(f"File {file_path} not found.")
        
# Check if data was loaded
if not bedgrphs:
    print("No data loaded. Exiting.")
    sys.exit()

# Count the occurrences of each row across all dataframes
counter = Counter()
for genome, df in bedgrphs.items():
    counter.update(df[['chromosome', 'start', 'end']].apply(tuple, axis=1))

# Identify rows that occur in more than 3 dataframes
common_rows = {row for row, count in counter.items() if count > 3}
print(f"Number of common rows found: {len(common_rows)}")

# Create new dataframes excluding common rows and save them
for genome, df in bedgrphs.items():
    df['combined'] = df[['chromosome', 'start', 'end']].apply(tuple, axis=1)
    filtered_df = df[~df['combined'].isin(common_rows)]
    filtered_df = filtered_df.drop(columns='combined')

    # Load chromosome mapping
    mapping = pd.read_csv('renameChrs_SGD_to_Ensembl.txt', sep='  ', header=None, names=['original', 'new'])

    # Convert chromosome names
    mapping_dict = dict(zip(mapping['original'], mapping['new']))
    filtered_df['chromosome'] = filtered_df['chromosome'].map(mapping_dict).fillna(filtered_df['chromosome'])

    # Drop rows with "Mito" chromosome names
    filtered_df = filtered_df[~filtered_df['chromosome'].str.contains("Mito")]

    # Save the filtered data with column labels
    output_file = f'{genome}.drop_ancestral.csv'
    filtered_df.to_csv(output_file, index=False, header=True)
    print(f"Saved filtered data to {output_file}")
