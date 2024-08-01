#!/bin/bash
# Search all files in current directory with the rg.md.bam suffix:
files=$(find . -name "*rg.md.bam" | sort)

for file in $files; do
  echo $file 
  # Calculate the sum of coverage:
  sum=$(samtools depth $file | awk '{sum+=$3} END { print sum }')
  # Calculate the number of positions:
  NR=$(samtools depth $file | awk '{sum+=$3} END { print NR }')
  # Calculate the average coverage across the whole genome:
  avg=$(bc -l <<< "$sum / $NR")
  # Set a threshold for coverage N-fold above the average (N*$avg)
  mult_avg=$(bc -l <<< "8*$avg")
	
  # Make a bedgraph file only if average coverage is above a threshold N ($avg > N; usually N~20)
  if [ $NR -gt 0 ] && (( $(bc -l <<< "$avg > 17.0") )); then
    bedtools genomecov -ibam $file -bg > data.bedgraph

    # Extract file name for output
    base_name=$(basename "$file" .rg.md.bam)
    output_file="${base_name}.gt.bedgraph"

    # Check if coverage (column 4) is greater than threshold, output to file if yes:
    awk -v mult_avg="$mult_avg" '$4 > mult_avg' data.bedgraph > "$output_file"

  fi

  # Delete the temporary file
  rm -f data.bedgraph
done