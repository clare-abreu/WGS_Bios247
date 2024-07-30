To curate and annotate SNPs:

(1) Curate SNPs in IGV
	1.1 In IGV, go to the dropdown menu Tools -> Run Batch Script and run each of the txt files generated previously in the Filter folder. IGV should save the images in the IGV_snapshots folder (which you created). 

	1.2 Save a copy of the {}_manuallyCurated.csv file as {}_manuallyCurated_edited.csv. This is the spreadsheet that you will edit. Open ~100 images at a time and verify that the SNP looks real. A common reason for an erroneous SNP is that it appears in a repetitive region of the genome. Mark the real SNPs as TRUE in the manually_verified column of the spreadsheet. Add a column called "Heterozygous?" for marking SNPs that appear heterozygous (about half one nucleotide, half another) as TRUE. Add another column called "Notes" for noting anything else.


(2) Update the VCF file with filters and curation
	2.1 Run the first few sections of getSNPs.py, up until the part where you save a new VCF file called "WGS_AllSamples_highcoverage_withFilterInfo_Homozygous.vcf." You must go to Sherlock to annotate the SNPs with this file. Transfer WGS_AllSamples_highcoverage_withFilterInfo_Homozygous.vcf to your Output folder, as well as the file renameChrs_SGD_to_Ensembl.txt, and run the following commands:

ml biology; ml bcftools

bcftools view -i'to_include=1' Output/WGS_highcoverage_withFilterInfo_Homozygous.vcf > Output/WGS_Homozygous_filtered.vcf

bcftools annotate --rename-chrs renameChrs_SGD_to_Ensembl.txt Output/WGS_Homozygous_filtered.vcf > Output/WGS_Homozygous_filtered.renamedChrs.vcf

bcftools norm -m "-" Output/WGS_Homozygous_filtered.renamedChrs.vcf > Output/WGS_Homozygous_filtered.renamedChrs.biAllelic.vcf


Then, we use snpEff to annotate the SNPs (note that we do this a second time for SNPs we labeled heterozygous during IGV curation):

ml java/11.0.11

java -Xmx8g -jar /home/groups/dpetrov/SOFTWARE/snpEff/snpEff.jar -v -c /home/groups/dpetrov/SOFTWARE/snpEff/snpEff.config R64-1-1.99 Output/WGS_Homozygous_filtered.renamedChrs.biAllelic.vcf > Output/WGS_Homozygous_filtered.renamedChrs.biAllelic.ann.vcf

java -jar /home/groups/dpetrov/SOFTWARE/snpEff/SnpSift.jar extractFields Output/WGS_Homozygous_filtered.renamedChrs.biAllelic.ann.vcf CHROM POS REF ALT "ANN[0].GENE" "ANN[0].EFFECT" "ANN[0].HGVS_C" "ANN[0].HGVS_P" "GEN[*].GT" > WGS_Homozygous_calledSNPs.tab


(3) Make a database of SNPs

Run SNPs_by_sample.py