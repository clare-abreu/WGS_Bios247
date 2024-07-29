#!/bin/bash

filelist=yeast_alignment_samples_all.inp

# Remove existing output files
rm -f found_bam.tmp
rm -f missing_bam.tmp

for file in $(cat "$filelist"); do
  if [ -e "Output/${file}.sorted.rg.md.bam" ]; then
    echo "${file:0:21}	Output/${file}.sorted.rg.md.bam" >> found_bam.tmp
  else
    echo "${file:0:21}" >> missing_bam.tmp
  fi
done