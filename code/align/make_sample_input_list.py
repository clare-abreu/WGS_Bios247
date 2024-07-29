#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 20:32:59 2023

@author: clare
"""

# Make list of all samples
# This will be the input for yeast_alignment.sbatch

# Create a list of samples with the format WGS-March23-143-A-B
with open(r'yeast_alignment_samples_all.inp', 'w') as fp:
    for a in range(1, 6):  # A ranges from 01 to 05
        for b in range(1, 97):  # B ranges from 01 to 96
            a_str = f"{a:02}"  # Format A with leading zeros
            b_str = f"{b:02}"  # Format B with leading zeros
            sample_name = f"WGS-March23-143-{a_str}-{b_str}"
            # Write the sample name
            fp.write(sample_name + "\n")

    print('Done')