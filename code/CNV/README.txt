To check for copy number variation:

1. Run the shell script record_hicov_regions.sh in the same directory where your rg.md.bam files are located. Set the N-fold threshold above average coverage for looking for copy-number variation. This will produce bedgraph files, where the first column is the chromosome, the second column is the start position, the third column the end position, and the fourth column the coverage.

2. In order to rule out copy number variation that is common to all mutants and hence from the ancestor, we run drop_common_CNVs.py to eliminate this data. This code produces .csv files (examples can be found in ../data/CNV_candidates