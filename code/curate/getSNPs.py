#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 14 23:17:09 2023

@author: clare
"""
#%% Import:
import numpy as np
import pandas as p
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import colors
# from adjustText import adjust_text
from Bio import SeqIO
from Bio.Seq import Seq

import gzip
from cyvcf2 import VCF

from itertools import combinations

sns.set_style('white')
sns.set_style('ticks')
sns.set_color_codes()

# %% Functions:
    
def get_vcf_names(vcf_path):
    with gzip.open(vcf_path, "rt") as ifile:
        for line in ifile:
            if line.startswith("#CHROM"):
                vcf_names = [x.strip('\n') for x in line.split('\t')]
                break
    ifile.close()
    return vcf_names

def hard_filter(variant,hard_filters):
    
    variant_annotations = {annotation:variant.INFO.get(annotation) for annotation in hard_filters.keys()}
    variant_annotations = {key:(val if val!=None else np.nan) for key,val in variant_annotations.items()}
    
    filter_out = np.any([variant_annotations['QD']<hard_filters['QD'],
            (variant_annotations['FS']>hard_filters['FS']),
            (variant_annotations['SOR']>hard_filters['SOR']),
            (variant_annotations['MQ']<hard_filters['MQ']),
            (variant_annotations['MQRankSum']<hard_filters['MQRankSum']),
            (variant_annotations['ReadPosRankSum']<hard_filters['ReadPosRankSum'])])
    
    return filter_out

#%% Read vcf:
file_name = 'WGS_AllSamples.vcf.gz'  ## Insert your high-coverage vcf file name here

names = get_vcf_names(file_name)
vcf = p.read_csv(file_name,compression='gzip', comment='#', delim_whitespace=True, header=None, names=names)


#%%
### HOMOZYGOUS ONLY
### UPDATE THE VCF FILE WITH OUR FILTERS ADDED and manual curation!!!


variable_sites_df = p.read_csv(f'{file_name.strip(".vcf.gz")}_manuallyCurated_edited.csv')

from cyvcf2 import VCF, Writer

fname = "WGS_highcoverage_withFilterInfo_Homozygous.vcf"
vcf = VCF(file_name)

vcf.add_info_to_header({'ID': 'to_include', 'Description': 'If we want to include in downstream analysis (removing hard filters, ancestral genotypes)',
    'Type':'Float', 'Number': '1'})
w = Writer(fname, vcf)

variable_sites_str = [str(v) for v in variable_sites_df['variant_string'].values]
place=0
for v in vcf:
    stripped_v = str(v).strip('\n')
    str_v = str(v)
    if str_v in variable_sites_df['variant_string'].values:
        
        if variable_sites_df[variable_sites_df['variant_string']==str_v]['manually_verified'].values[0]=='TRUE':
            if variable_sites_df[variable_sites_df['variant_string']==str_v]['Heterozygous?'].values[0]!=True:
                v.INFO["to_include"] = 1
                print('including')
                place+=1
        else:
            v.INFO["to_include"] = 0
            #print(variable_sites_df[variable_sites_df['variant_string']==str_v]['manually_verified'].values[0])
            #print(variable_sites_df[variable_sites_df['variant_string']==str_v]['Heterozygous?'].values[0])
    else:
        v.INFO["to_include"] = 0
    w.write_record(v)

w.close(); vcf.close()

print(place)


#%%
### HETEROZYGOUS ONLY
### UPDATE THE VCF FILE WITH OUR FILTERS ADDED and manual curation!!!


variable_sites_df = p.read_csv(f'{file_name.strip(".vcf.gz")}_manuallyCurated_edited.csv')

from cyvcf2 import VCF, Writer

fname = "WGS_highcoverage_withFilterInfo_Heterozygous.vcf"
vcf = VCF(file_name)

vcf.add_info_to_header({'ID': 'to_include', 'Description': 'If we want to include in downstream analysis (removing hard filters, ancestral genotypes)',
    'Type':'Float', 'Number': '1'})
w = Writer(fname, vcf)

variable_sites_str = [str(v) for v in variable_sites_df['variant_string'].values]
place=0
for v in vcf:
    stripped_v = str(v).strip('\n')
    str_v = str(v)
    if str_v in variable_sites_df['variant_string'].values:
        
        if variable_sites_df[variable_sites_df['variant_string']==str_v]['manually_verified'].values[0]=='TRUE':
            if variable_sites_df[variable_sites_df['variant_string']==str_v]['Heterozygous?'].values[0]==True:
                v.INFO["to_include"] = 1
                print('including')
                place+=1
        else:
            v.INFO["to_include"] = 0
            #print(variable_sites_df[variable_sites_df['variant_string']==str_v]['manually_verified'].values[0])
            #print(variable_sites_df[variable_sites_df['variant_string']==str_v]['Heterozygous?'].values[0])
    else:
        v.INFO["to_include"] = 0
    w.write_record(v)

w.close(); vcf.close()

print(place)

#%%
# HERE, you must annotate in the command line. 
