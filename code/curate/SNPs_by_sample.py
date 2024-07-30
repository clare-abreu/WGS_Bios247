#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 11:44:26 2023

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


def get_sample_info(names,calledSNPs):
    sample_names = np.array(names[9:])
    
    sample_snps = {}
    sample_allinfo = {}
    
    all_alt_calls = []
    
    max_snps_per_sample = 0
    
    for sample in sample_names:
        
        this_sample_snps = calledSNPs[calledSNPs[sample]!='0/0']
        
        this_sample_snps['LOC'] = [f'chr{chrom}:{loc}' for (chrom,loc) in zip(this_sample_snps['CHROM'],this_sample_snps['POS'])]
        
        if len(this_sample_snps)> max_snps_per_sample:
            max_snps_per_sample = len(this_sample_snps)
        
        for (chrom,pos,ref,alt,gene,effect,hgvs_c,hgvs_p,call) in zip(this_sample_snps['CHROM'].values,
                                                                      this_sample_snps['POS'].values,
                                                                      this_sample_snps['REF'].values,
                                                                      this_sample_snps['ALT'].values,
                                                                      this_sample_snps['GENE'].values,
                                                                      this_sample_snps['EFFECT'].values,
                                                                      this_sample_snps['HGVS_C'].values,
                                                                      this_sample_snps['HGVS_P'].values,
                                                                      this_sample_snps[sample].values):
            
            all_alt_calls.append([sample,chrom,pos,ref,alt,gene,effect,hgvs_c,hgvs_p,call])
        
        
        sample_snps[sample] = [f'{gene}:{effect}:{mutation_zygosity(call)}' for (gene,effect,call) in zip(this_sample_snps['GENE'].values,
                                                                                       this_sample_snps['EFFECT'].values,
                                                                                       this_sample_snps[sample].values)]
        
        
        sample_allinfo[sample] = '~'.join([f'{chrom}:{pos}:{ref}:{alt}:{gene}:{effect}:{hgvs_c}:{hgvs_p}:{call}' for (chrom,pos,ref,alt,gene,effect,hgvs_c,hgvs_p,call) in zip(
                                                                                       this_sample_snps['CHROM'].values,
                                                                                        this_sample_snps['POS'].values,
                                                                                        this_sample_snps['REF'].values,
                                                                                        this_sample_snps['ALT'].values,
                                                                                       this_sample_snps['GENE'].values,
                                                                                       this_sample_snps['EFFECT'].values,
                                                                                       this_sample_snps['HGVS_C'].values,
                                                                                        this_sample_snps['HGVS_P'].values,
                                                                                           this_sample_snps[sample].values)])
        return(all_alt_calls,sample_snps,sample_allinfo)



def mutation_zygosity(call):
    
    if (call == '0/1') or (call == '1/0') or (call == '0|1') or (call == '1|0') or (call == '0/2') or (call == '0|2'):
        return 'HET'
    elif (call == '1/1') or (call == '1|1') or (call == '2/2') or (call == '2|2'):
        return 'HOM'

#%% Read vcf:
file_name = 'WGS_AllSamples_20Xcoverage.vcf.gz'  ## Insert your high-coverage vcf file name here

names = get_vcf_names(file_name)

#%% HOMOZYGOUS:

sample_columns = names[9:]
header_list = ['CHROM','POS','REF','ALT','GENE','EFFECT','HGVS_C','HGVS_P'] + sample_columns
calledSNPs = p.read_csv('WGS_Homozygous_calledSNPs.tab',delimiter='\t',names=header_list)
calledSNPs = calledSNPs.drop(0,0)

calledSNPs.to_csv('WGS_Homozygous_calledSNPs.csv',index=False)
df_calledSNPs = calledSNPs

#%%

sample_names = np.array(names[9:])

sample_snps = {}
sample_allinfo = {}

all_alt_calls = []

max_snps_per_sample = 0

for sample in sample_names:
    
    this_sample_snps = calledSNPs[calledSNPs[sample]!='0/0']
    
    this_sample_snps['LOC'] = [f'chr{chrom}:{loc}' for (chrom,loc) in zip(this_sample_snps['CHROM'],this_sample_snps['POS'])]
    
    if len(this_sample_snps)> max_snps_per_sample:
        max_snps_per_sample = len(this_sample_snps)
    
    for (chrom,pos,ref,alt,gene,effect,hgvs_c,hgvs_p,call) in zip(
                                                                                   this_sample_snps['CHROM'].values,
                                                                                    this_sample_snps['POS'].values,
                                                                                    this_sample_snps['REF'].values,
                                                                                    this_sample_snps['ALT'].values,
                                                                                   this_sample_snps['GENE'].values,
                                                                                   this_sample_snps['EFFECT'].values,
                                                                                   this_sample_snps['HGVS_C'].values,
                                                                                    this_sample_snps['HGVS_P'].values,
                                                                                       this_sample_snps[sample].values):
        
        all_alt_calls.append([sample,chrom,pos,ref,alt,gene,effect,hgvs_c,hgvs_p,call])
    
    
    sample_snps[sample] = [f'{gene}:{effect}:{mutation_zygosity(call)}' for (gene,effect,call) in zip(this_sample_snps['GENE'].values,
                                                                                   this_sample_snps['EFFECT'].values,
                                                                                   this_sample_snps[sample].values)]
    
    
    sample_allinfo[sample] = '~'.join([f'{chrom}:{pos}:{ref}:{alt}:{gene}:{effect}:{hgvs_c}:{hgvs_p}:{call}' for (chrom,pos,ref,alt,gene,effect,hgvs_c,hgvs_p,call) in zip(
                                                                                   this_sample_snps['CHROM'].values,
                                                                                    this_sample_snps['POS'].values,
                                                                                    this_sample_snps['REF'].values,
                                                                                    this_sample_snps['ALT'].values,
                                                                                   this_sample_snps['GENE'].values,
                                                                                   this_sample_snps['EFFECT'].values,
                                                                                   this_sample_snps['HGVS_C'].values,
                                                                                    this_sample_snps['HGVS_P'].values,
                                                                                       this_sample_snps[sample].values)])
#%%
verbose_calls = p.DataFrame(all_alt_calls,columns=['sample','chrom','pos','ref','alt','gene','effect','hgvs_c','hgvs_p','call'])
verbose_calls.sort_values(['chrom','pos','ref','alt','sample']).to_csv('WGS_Homozygous_SNPs_verbose_by_call.csv',index=False)

sample_oriented_snps = p.DataFrame(columns=['sample','mutations','all_information'])
sample_oriented_snps['sample'] = sample_snps.keys()
sample_oriented_snps['mutations'] = sample_snps.values()
sample_oriented_snps['all_information'] = sample_allinfo.values()

sample_oriented_snps.to_csv('WGS_Homozygous_SNPs_by_sample.csv',index=False)


#%% HETEROZYGOUS:

sample_columns = names[9:]
header_list = ['CHROM','POS','REF','ALT','GENE','EFFECT','HGVS_C','HGVS_P'] + sample_columns
hetero_calledSNPs = p.read_csv('WGS_Heterozygous_calledSNPs.tab',delimiter='\t',names=header_list)
hetero_calledSNPs = hetero_calledSNPs.drop(0,0)

hetero_calledSNPs.to_csv('WGS_Heterozygous_calledSNPs.csv',index=False)

[all_alt_calls,sample_snps,sample_allinfo] = get_sample_info(names,hetero_calledSNPs)

verbose_calls = p.DataFrame(all_alt_calls,columns=['sample','chrom','pos','ref','alt','gene','effect','hgvs_c','hgvs_p','call'])
verbose_calls.sort_values(['chrom','pos','ref','alt','sample']).to_csv('WGS_Heterozygous_SNPs_verbose_by_call.csv',index=False)

sample_oriented_snps = p.DataFrame(columns=['sample','mutations','all_information'])
sample_oriented_snps['sample'] = sample_snps.keys()
sample_oriented_snps['mutations'] = sample_snps.values()
sample_oriented_snps['all_information'] = sample_allinfo.values()

sample_oriented_snps.to_csv('WGS_Heterozygous_SNPs_by_sample.csv',index=False)





