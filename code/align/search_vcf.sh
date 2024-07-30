#!/bin/bash

filelist=yeast_alignment_samples_all.inp

# Remove existing output files
rm -f found_vcf.tmp
rm -f missing_vcf.tmp

for file in $(cat "$filelist"); do
  if [ -e "Output/${file}.g.vcf.gz" ]; then
    echo "${file:0:21}	Output/${file}.g.vcf.gz" >> found_vcf.tmp
  else
    echo "${file:0:21}" >> missing_vcf.tmp
  fi
done