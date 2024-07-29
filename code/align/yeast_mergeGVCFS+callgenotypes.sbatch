#!/usr/bin/env bash
#SBATCH -J gtak_mergeall
#SBATCH -p hns,dpetrov,normal,owners
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -t 2-00:00
#SBATCH --array=1-1
#SBATCH --mem-per-cpu=2G
#SBATCH --requeue
#SBATCH -o SlurmFiles/slurm-%A_%a_%x.out

ml system
ml ncurses/6.0
ml biology
ml bwa/0.7.17
ml java

ref_genome=$"/home/groups/dpetrov/clare/S288C_reference/S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fasta"
dir=$"RawData/"
out_dir=$"Output/"

### Merge GVCF files

# (version for many files) - need to make a sample map
/home/groups/dpetrov/SOFTWARE/gatk-4.2.0.0/gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport -L intervals.list --sample-name-map yeast_samples.sample_map --genomicsdb-workspace-path "$out_dir"genomicsdb


### Call genotypes for all samples at the same time.
/home/groups/dpetrov/SOFTWARE/gatk-4.2.0.0/gatk --java-options "-Xmx4g" GenotypeGVCFs -R "$ref_genome" -V gendb://"$out_dir"/genomicsdb -O "$out_dir"WGS_AllSamples.vcf.gz
