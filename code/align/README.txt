To align WGS data:

(1) Move sequencing data (fastq R1 and R2 files) to Sherlock and do any helpful renaming.

(2) Put your reference genome in your directory. 

(3) Process files using gatk:
	
	3.1 Run yeast_alignment.sbatch. 
Split up if necessary into plates (96 parallel runs per sbatch), e.g. yeast_alignment_P1.sbatch with input file yeast_alignment_P1.inp. These input files should contain the names of the fastq files with no suffixes using the script make_sample_input_list.py, and each list has 96 (x2 for forward/reverse) files, so I set the sbatch to -array=1-96.

	3.2 Check for completion.
This may be an issue with Sherlock timing out (but still telling you the job completed). Not all files may have .g.vcf.gz and .sorted.rg.md.bam output files. Use the scripts search_vcf.sh and search_bam.sh to search for all files in the Output folder (modified to look for your particular files). Redo alignment by running yeast_alignment_redo.sbatch, which uses missing_bam.tmp as input. If you're missing some of both types of files, combine your missing samples into one list to use as input for yeast_alignment_redo.sbatch. Continue to run and re-run (adjusting #SBATCH --array=1-N each time) until this process is complete and you have no more missing output files.

	3.3 Run yeast_mergeGVCFs+callgenotypes.sbatch. This merges all the files together and uses GATK to call genotypes on the full set.
There input yeast_samples.sample_map is required here- it lists all of the files in Output/. I formatted found_bam.tmp, so you can just copy found_bam.tmp to yeast_samples.sample_map. Also required is intervals.list. The output is WGS_AllSamples.vcf.gz. Note: Do NOT make a directory  called genomicsdb in Output/ before running this. If you get an error suggesting this, it may be because you don't have intervals.list in the current directory.

	3.4 Calculate average coverage for each sample, keeping only samples above your favorite sample depth (I kept everything above 20X coverage on average). The script check_cov.sh saves the results. Note: before running it you must load ncurses with ml system; ml ncurses.

	3.5 Once you have the list of high-coverage files, include only those in the output vcf file using bcftools. 
First load biology and bcftools (ml biology; ml bcftools). Then this command: 
bcftools view -S hi_cov_file_names.tmp Output/WGS_AllSamples.vcf.gz > Output/WGS_AllSamples_20Xcoverage.vcf
And zip:
gzip Output/WGS_AllSamples_20Xcoverage.vcf

(4) Filter data
Proceed to Filter folder and README.