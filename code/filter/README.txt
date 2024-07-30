To filter data:

(1) Filter high-coverage samples to find SNPs
Run filter_data.py, which reads the zipped highcoverage.vcf file, filters SNPs to exclude ancestral genotypes as well as other filters, and generates scripts for IGV. You must run these scripts locally with IGV, which means you must download all of the sorted.rg.md.bam files to a local folder or a remote folder that is accessible from the desktop (they can be up to 500MB each), along with their index files (ending in sorted.rg.md.bam.bai). filter_data.py will also save a csv file with the list of filtered SNPs called {}_manuallyCurated.csv, which you will save a copy of and edit during the curation step.

(2) Curate the SNPs
In IGV, go to the dropdown menu Tools -> Run Batch Script and run each of the txt files generated from above. Note that the code for generating these txt files will need to be modified based on how many candidate SNPs you have. Proceed to the curate/ folder and README.