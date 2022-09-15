# Scripts

This folder contains scripts used in this study.

Each of these scripts may have also been run on another species, but using the nearly identical code. These scripts have been modified slightly from when they were run to remove the sbatch arguments and to add more comments to describe input and other relevant information.

Within each of the headings below, the scripts were run in the order listed. Some scripts were run on each species independently, others were run for all samples.

## Preprocessing Reads

**01_rw_trim_array.sh** - Trim raw fastq files and run fastqc on trimmed reads  
**02_rw_align_genome.sh** - Indexes NARW reference genome  
**03_rw_map_min1Mb.sh** - Map trimmed reads to reference genome of scaffold 1Mbp or longer  
**04_rw_sortsam.sh** - Sort the sam file and save as a .bam file  
**05_rw_merge_mark.sh** - Merges the .bam files across lanes and marks duplicates  
**06_qualimap.sh** - Runs qualimap to assess mapping performance across all samples  
**07_rw_haplotypecaller_gvcf.sh** - Run GATK's HaplotypeCaller in GVCF mode  
**08_narw_genomicsDBimport_1.sh** - Creates a genomics database for each scaffold. For this step, I had to create a script for each scaffold and run it from the scaffold directory.  
**09_genotypegvcf_on_narw.sh** - Calls haplotypes across individuals into a multi-individual vcf  
  
## Filtering VCF and Preliminary Data Exploration  

**10_vcf_filtering.sh** - Perform initial quality filtering on vcf files  
**11_KING_PLINK_Relatedness.sh** - Uses KING to calculate kinship coefficients between individuals  
**12_filter_relatedness.sh** - Remove related individuals and sex scaffolds  
**13_LDdecay.sh** - Calculate the decay in linkage with distance  
**14_PLINK_popstructure.sh** - Create the eigen vectors for population structure and inputs for ADMIXTURE  
**15_narw_admixture.sh** - Run ADMIXTURE for population sturcture analyses  
**16_polarizeVCF.py** - A custom python script by Josquin Daron to assign ancestral alleles (Python 3.8.1) *CANNOT BE USED WITHOUT PERMISSION*  
**16_run_polarize.sh** - Running the polarize python script on the cluster  
**17_roh_python.py** - Python script to calculate FROH and their abundance (Python 3.8.1)  
**17_run_roh_python.sh** - Bash script to run python script for FROH calculations  
**17_len_roh_python.py** - Python script to calculate the length of each ROH (Python 3.8.1)   
**17_run_len_roh.sh** - Bash script to run python script for sength of ROH  
  
## Demogrphic History Analyses

**18_phasing_vcf.sh** - Phase a VCF using SHAPEIT2  
**19_create_bed_vcf_inputs4msmc.sh** - Creates the mappability bed files and final vcf file used by generatemultihetsep.py needs for MSMC2  
**20_create_msmc_input_narw_set2** - Runs generatemultihetsep.py on a few samples to create MSMC2 input files  
**21_narw_stairwayplot_blueprint** - Sample blueprint file for stairway plot  
**21_create_stairwayplot_blueprint.sh** - Code to create bash script for stairway plot to run  
**21_run_stairwayplot_blueprint.sh** - Run Stairway plot as a scheduled job  
**22_IBDSeq.sh** - Run IBDSeq to generate IBD blocks  
**23_ibdne_narw.sh** - Run IBDNe  
**24_run_msmc2_srw.sh** - Run MSMC2 on a few samples  
**24_run_msmc2_all_between.sh** - Sample script for running MSMC on between species comparisons ahead of running MSMC-IM  
**25_combine_msmc-im.sh** - Run MSMC-IM   
