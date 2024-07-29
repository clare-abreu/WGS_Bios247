#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 14 14:38:10 2023

@author: clare

This code was adapted from Grant Kinsler's code.

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
    ## Reads the header of the gzipped VCF file to extract sample names.
    with gzip.open(vcf_path, "rt") as ifile:
        for line in ifile:
            if line.startswith("#CHROM"):
                vcf_names = [x.strip('\n') for x in line.split('\t')]
                break
    ifile.close()
    return vcf_names

def hard_filter(variant,hard_filters):
    ## Applies filters to distinguish between true genetic variants and artifacts.
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
file_name = 'WGS_AllSamples_20Xcoverage.vcf.gz'  ## Insert your high-coverage vcf file name here

names = get_vcf_names(file_name)
vcf = p.read_csv(file_name,compression='gzip', comment='#', delim_whitespace=True, header=None, names=names)

#%%
# the Reference genome has weird names for the chromosomes
# so here's a way to map between these names and human-readable ones

chromosome_map = {'chrI':'ref|NC_001133|',
'chrII':'ref|NC_001134|',
'chrIII':'ref|NC_001135|',
'chrIV':'ref|NC_001136|',
'chrV':'ref|NC_001137|',
'chrVI':'ref|NC_001138|',
'chrVII':'ref|NC_001139|',
'chrVIII':'ref|NC_001140|',
'chrIX':'ref|NC_001141|',
'chrX':'ref|NC_001142|',
'chrXI':'ref|NC_001143|',
'chrXII':'ref|NC_001144|',
'chrXIII':'ref|NC_001145|',
'chrXIV':'ref|NC_001146|',
'chrXV':'ref|NC_001147|',
'chrXVI':'ref|NC_001148|',
'chrMito':'ref|NC_001224|'}

chromosome_map_rev = {val:key for (key,val) in chromosome_map.items()}

#%%
hard_thresh_annotations = ['QD', 'FS', 'SOR', 'MQ', 'MQRankSum','ReadPosRankSum']
hard_thresh_lists = {annotation:[] for annotation in hard_thresh_annotations}
vcf = VCF(file_name)
for variant in vcf:
    for annotation in hard_thresh_annotations:
        
        hard_thresh_lists[annotation].append(variant.INFO.get(annotation))
        
#%% The VCF file we have so far is unfiltered. In order to get a tractible number of SNPs more likely to be real, we need to impose some hard filters on the data. 

### HARD FILTERS
# I made some informed changes based on the plots below. 
# I started with the GATK default filters and changed accordingly

## GATK DEFAULT FILTERS
QD_thresh = 2 # get rid of less than; GATK default is 2
FS_thresh = 60 # get rid of greater than; GATK default is 60
SOR_thresh = 3 # get rid of greater than; GATK default is 3
MQ_thresh = 40 # get rid of less than; GATK default is 40
MQRankSum_thresh = -12.5 # get rid of less than; GATK default is -12.5
ReadPosRankSum_thresh = -8.0 # get rid of less than; GATK default is -8.0

## MORE CONSERVATIVE FILTERS
#QD_thresh = 5 # get rid of less than; GATK default is 2
#FS_thresh = 60 # get rid of greater than; GATK default is 60
#SOR_thresh = 3 # get rid of greater than; GATK default is 3
#MQ_thresh = 50 # get rid of less than; GATK default is 40
#MQRankSum_thresh = -3.0 # get rid of less than; GATK default is -12.5
#ReadPosRankSum_thresh = -5 # get rid of less than; GATK default is -8.0

# put them in a dictionary that will be useful later
hard_filters = {'QD':QD_thresh,
                'FS':FS_thresh,
                'SOR':SOR_thresh,
                'MQ':MQ_thresh,
                'MQRankSum':MQRankSum_thresh,
                'ReadPosRankSum':ReadPosRankSum_thresh 
               }     

# %% Summarize filtering
hard_thresh_df = p.DataFrame.from_dict(hard_thresh_lists)

hard_filtered_out = hard_thresh_df[(hard_thresh_df['QD']<QD_thresh) |
                                   (hard_thresh_df['FS']>FS_thresh) |
                                   (hard_thresh_df['SOR']>SOR_thresh) |
                                   (hard_thresh_df['MQ']<MQ_thresh) |
                                   (hard_thresh_df['MQRankSum']<MQRankSum_thresh) |
                                   (hard_thresh_df['ReadPosRankSum']<ReadPosRankSum_thresh)]

hard_filtered = [True if entry in hard_filtered_out.index else False for entry in hard_thresh_df.index]
hard_thresh_df['hard_filtered'] = hard_filtered

print('TOTAL BEFORE FILTER:',len(hard_filtered))
print('TOTAL FILTERED OUT:',sum(hard_filtered))
print('TOTAL REMAINING AFTER FILTER:',len(hard_filtered)-sum(hard_filtered))

#%%
### CALL ANCESTRAL GENOTYPES

GQ_thresh = 70
n_high_quality = 5


filtered_out = []
ancestral_variants = []
variable_sites = []
low_quality_sites = []

vcf = VCF(file_name)
for variant in vcf: # or VCF('some.bcf')
    
    # exclude variants that hard filtered out
    if hard_filter(variant,hard_filters):
        continue 
    
    # exclude mitochondrial variants
    if chromosome_map_rev[variant.CHROM] == 'chrMito':
        continue
    
    # if all variants are the same, it's for sure an ancestral genotype!
    if len(set(variant.gt_types)) == 1:
        ancestral_variants.append(variant)
        continue
    else: 
        # it's a bit more complicated... we should find "majority" genotypes 
        high_quality_indices = np.where(variant.gt_quals > GQ_thresh)[0]
        
        if len(high_quality_indices) > n_high_quality:
            if len(set(variant.gt_types[high_quality_indices])) == 1:
                ancestral_variants.append(variant)
                continue
            elif np.sum(variant.gt_types != 0) >= 0.95*len(variant.gt_types): ## all the sites are non-reference
                ancestral_variants.append(variant)
                continue
            else:
                variable_sites.append(variant)
                continue
        else:
            low_quality_sites.append(variant)
            continue

# %% 
### output file for manual curation

variable_sites_dict = {
                        'IGV_search':[],
                       'manually_verified':[],
                       'location':[],
                       'REF':[],
                       'ALT':[],
                       'alt_calls':[],
                       'variant_string':[],
                        }

sample_names = np.array(names[9:])

for v,variant in enumerate(variable_sites):
    if len(str(variant))<20000:  # Sometimes there are strange long entries
        if variant.end > 50: # If the site is less than 50bp from the beginning of the chromosome, you cannot subtract 50
            variable_sites_dict['IGV_search'].append(f'{chromosome_map_rev[variant.CHROM]}:{variant.end-50}-{variant.end+50}')
            variable_sites_dict['location'].append(f'{chromosome_map_rev[variant.CHROM]}:{variant.end}')
            variable_sites_dict['REF'].append(variant.REF)
            variable_sites_dict['ALT'].append(variant.ALT)

            alt_locs = np.where(variant.gt_types != 0)[0]
            alt_sample_calls = [f'{name}:{bases} ({qual}, {ref}/{alt})' for name,bases,qual,ref,alt in zip(sample_names[alt_locs],
                                                    variant.gt_bases[alt_locs],
                                                    variant.gt_quals[alt_locs],
                                                    variant.gt_ref_depths[alt_locs],
                                                    variant.gt_alt_depths[alt_locs])]

            variable_sites_dict['alt_calls'].append('\t'.join(alt_sample_calls))


            variable_sites_dict['variant_string'].append(str(variant))
            variable_sites_dict['manually_verified'].append('')
        else:
            variable_sites_dict['IGV_search'].append(f'{chromosome_map_rev[variant.CHROM]}:{variant.end-(variant.end-1)}-{variant.end+50}')
            variable_sites_dict['location'].append(f'{chromosome_map_rev[variant.CHROM]}:{variant.end}')
            variable_sites_dict['REF'].append(variant.REF)
            variable_sites_dict['ALT'].append(variant.ALT)

            alt_locs = np.where(variant.gt_types != 0)[0]
            alt_sample_calls = [f'{name}:{bases} ({qual}, {ref}/{alt})' for name,bases,qual,ref,alt in zip(sample_names[alt_locs],
                                                    variant.gt_bases[alt_locs],
                                                    variant.gt_quals[alt_locs],
                                                    variant.gt_ref_depths[alt_locs],
                                                    variant.gt_alt_depths[alt_locs])]

            variable_sites_dict['alt_calls'].append('\t'.join(alt_sample_calls))


            variable_sites_dict['variant_string'].append(str(variant))
            variable_sites_dict['manually_verified'].append('')

variable_sites_df = p.DataFrame(variable_sites_dict)
    
variable_sites_df.to_csv(f'{file_name.strip(".vcf.gz")}_manuallyCurated.csv')

# %% Generate scripts for IGV to save images of SNPs

# Add essentially randomly chosen controls that will appear at the bottom of each image for comparison
# Select >80X coverage samples from various starting genotypes:
samples_curated = p.read_csv(f'{file_name.strip(".vcf.gz")}_manuallyCurated.csv',index_col=0)
#print(samples_curated)
control_files = ['WGS_ReSeq-04-95', # Source: H2O2/NaCl
                 'WGS_ReSeq-04-89', # Source: Glu/Lac
                 'WGS_ReSeq-04-40', # Source: Gal
                 'WGS_ReSeq-04-08', # Source: Gal
                 'WGS_ReSeq-04-04'] # Source: NaCl

name_prefix = '/Users/clare/Documents/WGS_bam_files_5.23/'
name_suffix = '.sorted.rg.md.bam'

snapshot_dir = '/Users/clare/Documents/GitHub/FluctuatingEvolution/code/processing/WGS_data_processing_reseq_May2023/CurateAnnotate/IGV_snapshots/'

# I used intervals of 100, which is the number of images IGV will generate at once.
# The number corresponds to the number of rows in the {}_manuallyCurated.csv file above.
# intervals should go up to the last full set of 100:
intervals = [0,100,200,300,400]
for k in range(len(intervals)):
    with open(f'IGV_manualCuration_batchInput_{intervals[k]}.txt','w') as igv_batch:
        igv_batch.write('new\n')
        igv_batch.write('genome sacCer3\n')
    #     for index in variable_sites_df.index:
        for index in samples_curated.index[intervals[k]:intervals[k]+100]:
    #         igv_search = variable_sites_df['IGV_search'].values[index]
    #         these_names = [call.split(':')[0] for call in variable_sites_df['alt_calls'].values[index].split('\t')]
            igv_search = samples_curated['IGV_search'].values[index]
            these_names = [call.split(':')[0] for call in samples_curated['alt_calls'].values[index].split('\t')]
            #print(these_names)

            for control in control_files:
                if control not in these_names:
                    these_names.append(control)

            name_joined = ','.join([f'{name_prefix}{name}{name_suffix}' for name in these_names])

            igv_batch.write('new\n')
            igv_batch.write(f'load {name_joined}\n')
            igv_batch.write(f'snapshotDirectory {snapshot_dir}\n')
            igv_batch.write(f'goto {igv_search}\n')
            igv_batch.write(f'collapse\n')
            igv_batch.write(f'snapshot ind{index}_{igv_search}.png\n')

# Last few- the last incomplete set; in this case, it is after 500:
with open(f'IGV_manualCuration_batchInput_500.txt','w') as igv_batch:
        igv_batch.write('new\n')
        igv_batch.write('genome sacCer3\n')
    #     for index in variable_sites_df.index:
        for index in samples_curated.index[500:]:
    #         igv_search = variable_sites_df['IGV_search'].values[index]
    #         these_names = [call.split(':')[0] for call in variable_sites_df['alt_calls'].values[index].split('\t')]
            igv_search = samples_curated['IGV_search'].values[index]
            these_names = [call.split(':')[0] for call in samples_curated['alt_calls'].values[index].split('\t')]
            #print(these_names)

            for control in control_files:
                if control not in these_names:
                    these_names.append(control)

            name_joined = ','.join([f'{name_prefix}{name}{name_suffix}' for name in these_names])

            igv_batch.write('new\n')
            igv_batch.write(f'load {name_joined}\n')
            igv_batch.write(f'snapshotDirectory {snapshot_dir}\n')
            igv_batch.write(f'goto {igv_search}\n')
            igv_batch.write(f'collapse\n')
            igv_batch.write(f'snapshot ind{index}_{igv_search}.png\n')
