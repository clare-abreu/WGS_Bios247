# WGS pipeline

This WGS pipeline was adapted from Grant Kinsler's code, who based it on pipelines from others in the Petrov lab.

You should start with the align/ folder, where we align to the reference genome. Then move to the filter/ folder, where we filter genotypes using GATK filters and also exclude mitochondrial variants and ancestral genotypes. Finally, we generate scripts to run in IGV, and move to the curate/ folder, where we curate SNPs into a list based on our confirmations in IGV.

Note: The code in the CNV folder looks for copy number variation. You do not need to call variants to do this- just align to the reference genome. 
