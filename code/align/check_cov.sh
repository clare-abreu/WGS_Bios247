#!/bin/bash

# Find .rg.md.bam files in the Output/ directory
files=$(find Output/ -name "*rg.md.bam" | sort)

# Remove existing output files
rm -f hi_cov_aves.tmp
rm -f hi_cov_file_names.tmp
rm -f lo_cov_aves.tmp
  
for file in $files; do
  echo $file 
  sum=$(samtools depth $file | awk '{sum+=$3} END { print sum }')
  NR=$(samtools depth $file | awk '{sum+=$3} END { print NR }')
  avg=$(bc -l <<< "$sum / $NR")
  
  if [ $NR -gt 0 ] && (( $(bc -l <<< "$avg > 20.0") )); then
    printf "%s   %.2f\n" "$file" "$avg" >> hi_cov_aves.tmp
    printf "%s\n" "${file:7:21}" >> hi_cov_file_names.tmp  # Updated slice to account for "Output/" prefix
  else
    printf "%s   %.2f\n" "$file" "$avg" >> lo_cov_aves.tmp
  fi
done
