# Scripts

This folder contains scripts used in this study.

These scripts have been modified slightly from when they were run to remove the sbatch arguments and to add more comments to describe input and other relevant information.

Within each of the headings below, the scripts were run in the order listed. Some scripts were run on each species independently, others were run for all samples.

## Preparing Reference Genomes

**rw_index_genome.sh** - Indexes NARW reference genome
**srw_index_genome.sh** - Indexes SRW reference genome

## North Atlantic Right Whale Preprocessing

1. **rw_trim_array.sh** - Trim raw fastq files and run fastqc on trimmed reads
2. **rw_map_min1Mb.sh** - Map trimmed reads to reference genome of scaffold 1Mbp or longer
3. **rw_sortsam.sh** - Sort the sam file and save as a .bam file
4. **rw_merge_mark.sh** - Merges the .bam files across lanes and marks duplicates
5. **qualimap.sh** - Runs qualimap to assess mapping performance across all samples
6. **rw_haplotypecaller_gvcf.sh** - Run GATK's HaplotypeCaller in GVCF mode
7. *GENOMICSDB* - Creates a genomics database for each scaffold. For this step, I had to create a script for each scaffold and run it from the scaffold directory.
8. **rw_genotypegvcfs.sh** - Calls haplotypes across individuals into a multi-individual vcf

## Southern Right Whale Preprocessing
1. *Trimming performed on all samples in the same script rw_trim_array.sh*
2. **srw_map_min1Mb.sh** - Map trimmed reads to reference genome of scaffold 1Mbp or longer
3. **srw_sortsam.sh** - Sort the sam file and save as a .bam file
4. *Merging .bamfiles and Marking Dups done for all right whales with rw_merge_mark.sh above*
5. *Qualimap done for all samples using qualimap.sh above*
6. **srw_haplotypecaller_gvcf.sh** - Run GATK's HaplotypeCaller in GVCF mode
7. *GENOMICSDB* - Creates a genomics database for each scaffold. For this step, I had to create a script for each scaffold and run it from the scaffold directory.
8. **srw_genotypegvcfs.sh** - Calls haplotypes across individuals into a multi-individual vcf


### Bowhead Preprocessing
1. **bowhead_sra_to_fastq.sh** - Download and unpack Bowhead SRA into fastq.gz files
2. *Trimming performed on all samples in the same script rw_trim_array.sh*
3a. **bowhead_map_min1Mb.sh** - Map trimmed reads to reference genome of NARW
3b. **bowhead_srwmap_min1Mb.sh** - Map trimmed reads to reference genome of SRW
4a. **bowhead_sortsam.sh** - Sort the sam file and save as a .bam file of bowhead on NARW
4b. **bowhead_sortsam_srw** - Sort the sam file and save as a .bam file of bowhead on SRW


## Analyses

