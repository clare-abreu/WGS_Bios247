#!/usr/bin/env bash
#SBATCH -J bwa_reseq
#SBATCH -p hns,dpetrov,normal,owners
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -t 2-00:00
#SBATCH --array=1-96
#SBATCH --mem-per-cpu=2G
#SBATCH --requeue
#SBATCH -o SlurmFiles/slurm-%A_%a_%x.out

bn=$(sed -n "$SLURM_ARRAY_TASK_ID"p yeast_alignment_samples.inp)

ml system
ml ncurses/6.0
ml biology
ml bwa/0.7.17
ml java

ref_genome=$"/home/groups/dpetrov/clare/S288C_reference/S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fasta"
dir=$"../RawData/WGS_ReSeq/"
out_dir=$"Output/"

### create index of reference genome
bwa index "$ref_genome"
/home/groups/dpetrov/SOFTWARE/gatk-4.2.0.0/gatk CreateSequenceDictionary -R "$ref_genome"
samtools faidx "$ref_genome"

### align reads to reference, using 16 threads
bwa mem -t 16 "$ref_genome" "$dir$bn"_R1.fastq.gz "$dir$bn"_R2.fastq.gz | samtools view -bS - | samtools sort > "$out_dir$bn".sorted.bam

### fix read groups
java -jar /home/groups/dpetrov/SOFTWARE/picard-tools/picard.jar AddOrReplaceReadGroups I="$out_dir$bn".sorted.bam O="$out_dir$bn".sorted.rg.bam RGID="$bn" RGLB=lib"$bn" RGPL=illumina RGPU=unit"$bn" RGSM="$bn"
samtools index "$out_dir$bn".sorted.rg.bam

### mark duplicates
java -jar /home/groups/dpetrov/SOFTWARE/picard-tools/picard.jar MarkDuplicates I="$out_dir$bn".sorted.rg.bam O="$out_dir$bn".sorted.rg.md.bam M="$out_dir$bn".md.metrics.txt
samtools index "$out_dir$bn".sorted.rg.md.bam

### call variants
/home/groups/dpetrov/SOFTWARE/gatk-4.2.0.0/gatk --java-options "-Xmx4g" HaplotypeCaller -R "$ref_genome" -I "$out_dir$bn".sorted.rg.md.bam -ERC GVCF -O "$out_dir$bn".g.vcf.gz
