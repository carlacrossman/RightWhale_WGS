## Notes on WGS mapping ##


June 6, 2022

I am startig a proper log to act as a digital lab notebook.

To date, I have done the following on the WGS dataset:

**Trimming**
- I have trimmed 1 sample (EGL254-1) using TRIMMOMATIC. I would like to set it up using an array, but I don't know how to write it as a two variable array. I used the following code with a for loop, but I suggest the subsequent for my loop for multiple samples where lane_num is a file that contains lane numbers (L002, L003, L004), with some obvious changes to use different sample names.

```
#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --account=def-frasiert
#SBATCH --job-name=trim_sample
#SBATCH --mem-per-cpu=8000M
#SBATCH --cpus-per-task=8

SAMPLE="EGL254-1"

for LANE in L002 L003 L004;do \
java -Xmx2G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE \
	-threads 8 -phred33 \
	raw_fastq/$SAMPLE_2-2252940_S1_${LANE}_R1_001.fastq.gz \
	raw_fastq/$SAMPLE_2-2252940_S1_${LANE}_R2_001.fastq.gz \
	trimmed_reads/$SAMPLE-$LANE-FP.fq.gz \
	trimmed_reads/$SAMPLE-$LANE-FUP.fq.gz \
	trimmed_reads/$SAMPLE-$LANE-RP.fq.gz \
	trimmed_reads/$SAMPLE-$LANE-RUP.fq.gz \
	ILLUMINACLIP:reference/adapters.fa:2:30:15 LEADING:20 SLIDINGWINDOW:5:20 \
	AVGQUAL:30 MINLEN:36 2> trimmed_reads/$SAMPLE-$LANE.trim.out; done
```
- - - - - - - - - - - - - - - - - - - -
```
#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --account=def-frasiert
#SBATCH --job-name=trim_sample
#SBATCH --array=1-3
#SBATCH --mem-per-cpu=8000M
#SBATCH --cpus-per-task=8

SAMPLE="EGL254-1"
LANE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" lane_num)

java -Xmx2G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE \
	-threads 8 -phred33 \
	raw_fastq/$SAMPLE_2-2252940_S1_${LANE}_R1_001.fastq.gz \
	raw_fastq/$SAMPLE_2-2252940_S1_${LANE}_R2_001.fastq.gz \
	trimmed_reads/$SAMPLE-$LANE-FP.fq.gz \
	trimmed_reads/$SAMPLE-$LANE-FUP.fq.gz \
	trimmed_reads/$SAMPLE-$LANE-RP.fq.gz \
	trimmed_reads/$SAMPLE-$LANE-RUP.fq.gz \
	ILLUMINACLIP:reference/adapters.fa:2:30:15 LEADING:20 SLIDINGWINDOW:5:20 \
	AVGQUAL:30 MINLEN:36 2> trimmed_reads/$SAMPLE-$LANE.trim.out
```

**Quality Plots**
- I have performed FastQC on the 6 trimmed paired reads files for EGL254-1. I would like to build this into the script for the trimming, as I ran these in an interactive shell. This code ran for 30-45 mins per sample.

```
fastqc -o trimQC/ EGL254-1-L002.FP.q.gzf
```

I saved the plots to my computer in the RightWhale/WGS/test_sample/ folder as .html files. Today I compared them and while all of the general quality scores still looked great, the biggest difference was the overrepresented sequences were removed - so we are trimming the correct adapter sequences!

**Mapping**
- I am working on the mapping. I have run the mapping code a few times and I have ended with an out of memory error. It was run as an array for the three lanes and I had given the script 20 cores and 40 GB each. I set a run time to 48 hours, but it looks like BWA MEM only took about 24 hours. There are .bam files but according to samtools quickcheck, they are missing EOF and are therefore incomplete. The error looks to be from SortSam in GATK. Here is the tail of the slurm log output:
```
slurmstepd: error: Detected 1 oom-kill event(s) in StepId=35491929.batch. Some of your processes may have been killed by the cgroup out-of-memory handler.
```

I had tested the code using a pipe on a subset of the .fq file, and it ran fine. I was able to open it in IGV and all looked great. 

I have increased the number of cores to 30. I also gave each job 45GB of memory. I also eliminated the pipe so I saved the .sam file and used it as the input as well. I started this at ~9:30am. This finished running in about 3.5 hours using the following code:

```
#!/bin/bash
#SBATCH --job-name=$SAMPLE_alignment
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-3
#SBATCH --cpus-per-task=30
#SBATCH --mem=45G
#SBATCH --time=48:00:00

SAMPLE="EGL254-1"
LANE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" lane_num)

READGROUP="@RG\tID:${LANE}\tSM:${SAMPLE}\tLB:RW\tPU:${SAMPLE}\tCN:MCGILL\tPL:ILLUMINA"

bwa mem -M -t 30 \
  -R $READGROUP \
  reference/Eubalaena_glacialis_HiC.fasta \
  trimmed_reads/$SAMPLE-$LANE-FP.fq.gz \
  trimmed_reads/$SAMPLE-$LANE-RP.fq.gz > saved_results/$SAMPLE-$LANE.sam

java -Xmx2G -jar ${GATK_JAR} SortSam \
  -I saved_results/$SAMPLE-$LANE.sam \
  -O saved_results/$SAMPLE-$LANE-aln-sorted.bam \
  -SO coordinate \
  --CREATE_INDEX true --MAX_RECORDS_IN_RAM 500000
```

- I should find out why the pipe didn't work, but I could just run it this way and add rm saved_results/$SAMPLE-$LANE.sam file.

- I started to run the rw_align_min1Mb.sh which maps the reads to a reduced genome file. This was batch job *35628302* 
```
#!/bin/bash
#SBATCH --job-name=$SAMPLE_alignment
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-3
#SBATCH --cpus-per-task=30
#SBATCH --mem=20G
#SBATCH --time=06:00:00

LANE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" lane_num)
SAMPLE="EGL254-1"

bwa mem -M -t 30 \
  -R '@RG\tID:${LANE}\tSM:${SAMPLE}\tLB:RW\tPU:${SAMPLE}\tCN:MCGILL\tPL:ILLUMINA' \
  reference/Eubalaena_glacialis_HiC_min1Mb.fasta \
  trimmed_reads/$SAMPLE-$LANE-FP.fq.gz \
  trimmed_reads/$SAMPLE-$LANE-RP.fq.gz > saved_results/$SAMPLE-$LANE-aln_min1Mb.sam

java -Xmx2G -jar ${GATK_JAR} SortSam \
  -I saved_results/$SAMPLE-$LANE-aln_min1Mb.sam \
  -O saved_results/$SAMPLE-$LANE-aln-sorted_min1Mb.bam \
  -SO coordinate \
  --CREATE_INDEX true --MAX_RECORDS_IN_RAM 500000
```

**Inspecting mapping results**

- I saved EGL254-1-L002-aln-sorted.bam and its index file to my test_sample/ directory. In IGV the bam file looks correct (and like there will be plenty of variable sites!).

- As this looked good, I went ahead and deleted the intermediate .sam files.

### June 7, 2022 ###

**Read Mapping**
- I went to inspect the results from the read mapping to the reduced genome file. It doesn't appear that a .sam file was created. I wonder if I deleted the others accidentally as they were being written. Not sure what could have happened here. I am running that code again this morning. rw_align_min1Mb.sh 

I noticed it didnt have the correct assignment for READGROUP, so I fixed that before running. The job started as job: *35708635*.

- The job failed in less than 3 hours from an out of memory handler in SortSam. I think if I increased the memory to 45GB again, it would have been fine, 20GB was not enough. As I wrote this to same the itermediate .sam file, I will run GATKs SortSam on its own to create the sorted bam file. I started this as job: *35718380*

**QUALIMAP**

- I want to test how qualimap works on one .bam file. I am laoding the qualimap module from genpipes and then running it on an interactive node on EGL254-1-L002.

From the saved results folder:

```
qualimap bamqc -bam EGL254-1-L002-aln-sorted.bam -outfile EGL254-1-L002_bamqc.pdf
```

This took way too long to run in an interactive shell. Will assess program needs again and schedule the job to run later.

Running a shell script with qualimap multi-bamqc on the 6 bam files I have using the following code:

```
#!/bin/bash
#SBATCH --job-name=qualimap_testing
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=30
#SBATCH --mem=10G
#SBATCH --time=06:00:00

qualimap multi-bamqc -d multiqc_input \
        -outdir multiqc_test/ -nr 10000 \
        --java-mem-size=8G -r
```

Started job *35728110*

Downloaded .html report and images for EGL254-1-L002-aln-sorted. Looks like each one is taking about an hour to run. The smaller scaffold ones should go faster. This should be done around 8pm tonight. Hopefully it will have compiled the summary document too.


### June 8, 2022 ###

Today I am going to:
	- compare qualimap results
	- download the bowhead whale raw fastq files
	- merge .bam files
	- mark duplicates

**Comparing qualimap results**

Yesterday I ran qualimap on 1 sample for 3 lanes against the full and reduced reference genome. Mean coverage for each was around 17 - 18 suggesting we will get good depth once we merge lanes even after the duplicates are removed. The mean mapping quality for the whole reference genome was 16.7 whereas it jumped to 55 with the min1Mb reference genome file. These testing reports are saved to my local computer.

**Downloading bowhead genome files**

The genome and SRA etc information for the bowhead can be found on NCBI herhttps://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=run_browser&run=SRR1685383 or under the accession: SRR1685383

First I need to download the .sra file into a bowhead folder in the raw_fastq directory and load the sratool kit to unpack the fastq files from the sra file.

```
wget 	https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR1685383/SRR1685383

module load mugqic/sratoolkit/2.8.2-1
```

I will now use SRATool kit to convert the .sra file to fastq files that I can process alongside all of our other samples. We will need to decide whether or not we should align this to each right whale reference genome or not.

I will start with the following code. I am trying to run in an interactive shell for 2 hours with 8 cores.

```
fastq-dump --split-files --gzip -v SRR1685383.sra
```

This was taking too long, so I wrote the following script and ran it to unpack. Started as job *35882729*

```
#!/bin/bash
#SBATCH --job-name=bowhead_unpacking
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=10
#SBATCH --mem=10G
#SBATCH --time=24:00:00

fastq-dump --split-files --gzip -v SRR1685383.sra
```

This will give 2 fastq files, 1 for each read direction and I can sort them ad hoc by lane number which is stored as the 4th element in the fastq header.

**Merging .bam files**

I want to merge the .bam files for the test sample with both reference genomes. I will then mark duplicates, and run qualimap once more to keep statistics to back up why we dropped the small scaffolds.

I will use GATK MergeSamFiles for this. It claims to be the same as samtools merge and it is probably best to stay within the same toolkit for minor steps like this where it likely doesn't matter. Started the following code as a batch file for 6 hours as job *35885316*.

```
java -Xmx2G -jar ${GATK_JAR} MergeSamFiles \
  --INPUT EGL254-1-L002-aln-sorted.bam \
  --INPUT EGL254-1-L003-aln-sorted.bam \
  --INPUT EGL254-1-L004-aln-sorted.bam \
  --OUTPUT EGL254-1-merged.bam
  
java -Xmx2G -jar ${GATK_JAR} MergeSamFiles \
  --INPUT EGL254-1-L002-aln-sorted_min1Mb.bam \
  --INPUT EGL254-1-L003-aln-sorted_min1Mb.bam \
  --INPUT EGL254-1-L004-aln-sorted_min1Mb.bam \
  --OUTPUT EGL254-1-merged_min1Mb.bam
```

This finished in approx 2 hours.

I deleted .sam files in the saved_results/ directory.


**Marking duplicates**

I am marking duplicates in the meged .bam files.

```
#!/bin/bash
#SBATCH --job-name=bowhead_unpacking
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=10
#SBATCH --mem=10G
#SBATCH --time=6:00:00

SAMPLE="EGL254-1"

java -Xmx2G -jar ${GATK_JAR} MarkDuplicates \
  --REMOVE_DUPLICATES false --CREATE_INDEX true \
  -I saved_results/$SAMPLE-merged.bam \
  -O alignment/$SAMPLE-merged-nodup.bam \ --METRICS_FILE=alignment/$SAMPLE-merged.dup.metrics
  
java -Xmx2G -jar ${GATK_JAR} MarkDuplicates \
  --REMOVE_DUPLICATES false --CREATE_INDEX true \
  -I saved_results/$SAMPLE-merged_min1Mb.bam \
  -O alignment/$SAMPLE-merged-nodup_min1Mb.bam \ --METRICS_FILE=alignment/$SAMPLE-merged_min1Mb.dup.metrics
```

Started as job *35908418*. Note the name of the job was accidentally left as bowhead_unpack.

### June 9, 2022 ###

Today I would like to:
	- finish marking duplicates
	- prepare sample name files
	- move Eau fastq files and EGL254-1
	- make sure fq for EGL254-1 in raw_fastq folder
	- fine tuning scripts for array jobs
	
**Marking Duplicates**

It appears that after 30 minutes, the first mark duplicates file encountered an out of memory error (it had been given 10G). I am rerunning the same code (with a corrected name - mark_dups_moremem) and providing 45G of memory. I requested 6hrs for the job, but it will not likely take that long.

The job started as *36009649*.

This job ended after 90 mins again because of an outof memory handler. I reran the script, just with the reduced min1Mb dataset and with 60G and 20 cores. Started job *36025925*.

This looked like it was gettign hung up. I started it again as job *36033042* after increasing the amount of memory allowed by java to 20G. The job timed out after 3 hours, but was not complete. In the script, I had allowed 60GB of memory in the sbatch settings and it looks like that was the bottleneck. The code used 21GB of memory, was given 20 cores and cpu efficiency of 18%. I restarted with 20G and 6 hours (*36054251*). I also started the same code for marking duplicates for the whole genome (*36054267*). Oom error.

Tried again, with 60GB again *36058052* and *36058082*.

**Sample Directories**

In the projects/def-frasiert/RW_WGS/ directory, I saved a file with a list of sample names as *samples*. I also saved a file called *lane_num* which has a list of the lane numbers used for the right whales. There will need to be a different file for lane number for the bowhead run.

I then created directories for each sample in the raw_fastq where I will move the fastq.gz files form the raw_data folder. There should be 6 files (3 lanes by 3 directions in each directory).

I moved the files to the appropriate directory using the following code for each sample:
```
mv raw_reads/Q009281_Right-Whales/*/*/EGL308-1a_2-*R*_001.fastq.gz raw_fastq/EGL308-1a
```
I then checked that we had 6 files per sample as expected.

**Preparing scripts for multiple samples**

Here I will work to fine tune the scripts I have working for one sample and start laying out the issues I need to resolve for scaling up to multiple samples.

I will need 1 file for each of the following:

	- Trimming and Fastqc
	- Mapping to min1Mb genome and Sorting/Indexing Sam
	- Merging .bam files
	- Marking Duplicates
	- Qualimap

Below I will copy and paste the script I have written so far for each step and highlight the parts I know I will need to change. Each of these have been saved locally and I will copy to Cedar when I am have with them. As the lane numbers are different for the bowhead, I have to decide how I should account for this going forward and how that affects the following scripts.

*TRIMMING AND FASTQC*

The time and memory requirements for this script will need to be changed. i also need to change the input for TRIMMOMATIC to auto detect number of threads. I also need to make this a multi-variable array for sample name and take into account that I want to eliminate a middle portion of the fastq file name. This file is saved locally as *rw_trim_array.sh*.
```
#!/bin/bash
#SBATCH --account=def-frasiert
#SBATCH --job-name=trim_sample
#SBATCH --array=1-3
#SBATCH --mem-per-cpu=8000M
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00

SAMPLE="EGL254-1"
LANE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" lane_num)

java -Xmx2G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE \
    -threads 8 -phred33 \
    ~/projects/def-frasiert/RW_WGS/raw_fastq/$SAMPLE/${SAMPLE}_2-2252940_S1_${LANE}_R1_001.fastq.gz \
	~/projects/def-frasiert/RW_WGS/raw_fastq/$SAMPLE/${SAMPLE}_2-2252940_S1_${LANE}_R2_001.fastq.gz \
	~/scratch/trimmed_reads/$SAMPLE-$LANE-FP.fq.gz \
    ~/scratch/temp/$SAMPLE-$LANE-FUP.fq.gz \
    ~/scratch/trimmed_reads/$SAMPLE-$LANE-RP.fq.gz \
    ~/scratch/temp/$SAMPLE-$LANE-RUP.fq.gz \
    ILLUMINACLIP:~/projects/def-frasiert/RW_WGS/reference/adapters.fa:2:30:15 LEADING:20 SLIDINGWINDOW:5:20 \
    AVGQUAL:30 MINLEN:36 2> ~projects/def-frasiert/QC/std_out/$SAMPLE-$LANE.trim.out

fastqc -o ~projects/def-frasiert/QC/trimQC/ ~/scratch/trimmed_reads/$SAMPLE-$LANE-FP.fq.gz
fastqc -o ~projects/def-frasiert/QC/trimQC/ ~/scratch/trimmed_reads/$SAMPLE-$LANE-RP.fq.gz

rm ~/scratch/temp/$SAMPLE-$LANE-*.fq.gz
```

*MAPPING AND SORTING*

In this script, I know that for a single sample, 20GB was not enough, so I ran the SortSam as a separate script, but with enough memory, it should be fine. My SortSam test script alone, ran fine with 45G of memory. I also could not make this work with a pipe, likely due to the memory requirements, so I will write the intermediate .sam file and delete it at the end of the script.

This file will need to be corrected for memory and time needs, automating thread count usage in bwa and it will require tweaking to allow for a multivariable array. This script is saved locally as rw_map_min1Mb.sh

Note: In order to run this script for bowhead, the RG will need to be changed!
```
#SBATCH --job-name=map_reads
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-3
#SBATCH --cpus-per-task=30
#SBATCH --mem=20G
#SBATCH --time=06:00:00

LANE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" lane_num)
SAMPLE="EGL254-1"
READGROUP="@RG\tID:${LANE}\tSM:${SAMPLE}\tLB:RW\tPU:${SAMPLE}\tCN:MCGILL\tPL:ILLUMINA"

bwa mem -M -t 30 \
  -R $READGROUP \
  ~projects/def-frasiert/RW_WGS/reference/Eubalaena_glacialis_HiC_min1Mb.fasta \
  ~/scratch/trimmed_reads/$SAMPLE-$LANE-FP.fq.gz \
  ~/scratch/trimmed_reads/$SAMPLE-$LANE-RP.fq.gz > \
  ~/scratch/temp_mapping/$SAMPLE-$LANE-aln_min1Mb.sam \
  2> ~/projects/def-frasiert/RW_WGS/QC/std_out/$SAMPLE-$LANE-bwa.err

java -Xmx2G -jar ${GATK_JAR} SortSam \
  -I ~/scratch/temp_mapping/$SAMPLE-$LANE-aln_min1Mb.sam \
  -O ~/scratch/mapped_reads/$SAMPLE-$LANE-aln-sorted_min1Mb.bam \
  -SO coordinate \
  --CREATE_INDEX true --MAX_RECORDS_IN_RAM 500000 \
  2> ~/projects/def-frasiert/RW_WGS/QC/std_out/$SAMPLE-$LANE-SortSam.out

rm ~/scratch/temp_mapping/$SAMPLE-$LANE-*.sam
```

*MERGING BAM*

Again, need to correct time and memory needs. Maybe this could be piped into marked duplcates, but I didn't test that code yet. It was also not run with any array parameters so those may need to be checked. This is saved locally as *rw_merge_bam.sh*.

```
#!/bin/bash
#SBATCH --job-name=merge_bam
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-22%5
#SBATCH --cpus-per-task=20
#SBATCH --mem=10G
#SBATCH --time=06:00:00

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples)

java -Xmx2G -jar ${GATK_JAR} MergeSamFiles \
  --INPUT ~/scratch/mapped_reads/$SAMPLE-L002-aln-sorted_min1Mb.bam \
  --INPUT ~/scratch/mapped_reads/$SAMPLE-L003-aln-sorted_min1Mb.bam \
  --INPUT ~/scratch/mapped_reads/$SAMPLE-L004-aln-sorted_min1Mb.bam \
  --OUTPUT ~/scratch/merged_reads/$SAMPLE-merged_min1Mb.bam \
  2> ~/projects/def-frasiert/RW_WGS/QC/std_out/$SAMPLE-merge.out

```

*MARKING DUPLICATES*
 
Still to come


*QUALIMAP*

In order for this to run properly, it required a tab delimited input file I will call qualimap_input that describes the data.

e.g.
EGL254-1-L002   EGL254-1-L002-aln-sorted.bam    whole
EGL254-1-L003   EGL254-1-L003-aln-sorted.bam    whole
EGL254-1-L004   EGL254-1-L004-aln-sorted.bam    whole
EGL254-1-L002_min       EGL254-1-L002-aln-sorted_min1Mb.bam     min1Mb
EGL254-1-L003_min       EGL254-1-L003-aln-sorted_min1Mb.bam     min1Mb
EGL254-1-L004_min       EGL254-1-L004-aln-sorted_min1Mb.bam     min1Mb

```
#!/bin/bash
#SBATCH --job-name=qualimap_testing
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=30
#SBATCH --mem=10G
#SBATCH --time=06:00:00

qualimap multi-bamqc -d qualimap_input \
        -outdir qualimap_test/ -nr 10000 \
        --java-mem-size=8G -r
```

**Bowhead Reference**

The 92G 200 PE file didn't even complete after 24 hours. There were several other sra files available from the same project. I am downloading the Bowhead MP 10Kb genome now instead. It is only 17G. I will unpack in the same way, but hopefully it will be far more efficient.

```
get 	https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR1685387/SRR1685387

module load mugqic/sratoolkit/2.8.2-1
```

I will now use SRATool kit to convert the .sra file to fastq files that I can process alongside all of our other samples. We will need to decide whether or not we should align this to each right whale reference genome or not.

```
#!/bin/bash
#SBATCH --job-name=bowhead_unpacking
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=30
#SBATCH --mem=10G
#SBATCH --time=24:00:00

fastq-dump --split-files --gzip -v SRR1685387.sra
```

Submitted as job *36025535*.


**Splitting bowhead fastq into lanes**

Based on the NCBI database, I suspect that the bowhead files were run in 3 lanes (L003, L004 and L006). 

I found a one-line awk command that grabs the third variable after a : in the fastq headers which would be the lane field. It then takes the header from those lines and the following three lines and writes a new fastq file. I tested it on the small fastq locally.

```
awk 'BEGIN {FS = ":"} {lane=$4 ; print > "bowhead."lane".fastq" ; for (i = 1; i <= 3; i++) {getline ; print > "bowhead."lane".fastq"}}' < EGL254-1_tosplit.fastq
```

The files I have will be compressed, so I need to unzip and pipe into awk and then ideally write zip files. To tackle the first problem, I have the following I will run on the actual 2 read files:

```
#!/bin/bash
#SBATCH --job-name=bowhead_fastqsplit
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=20
#SBATCH --mem=10G
#SBATCH --time=6:00:00

zcat SRR1685387_1.fastq.gz | awk 'BEGIN {FS = ":"} {lane=$4 ; print > "SRR1685387.L00"lane"_R1_001.fastq" ; for (i = 1; i <= 3; i++) {getline ; print > "SRR1685387.L00"lane"_R1_001.fastq"}}' -

zcat SRR1685387_2.fastq.gz | awk 'BEGIN {FS = ":"} {lane=$4 ; print > "SRR1685387.L00"lane"_R2_001.fastq" ; for (i = 1; i <= 3; i++) {getline ; print > "SRR1685387.L00"lane"_R2_001.fastq"}}' -
``` 

As this provides, different output files, I want to be sure I am not zipping them together. I will run this and then gzip the output fastq files. Started as job *36055235*. As tbis was running, it didn't look right. I unzipped the start of the fastq file and it was not a true fastq file. Maybe it wasn't unpacked correctly, but I suspect something else is wrong with the files themselves.

I moved back to unpack the large .sra in the raw_fastq folder and started a job *36057178* with --gzip removed and I will be able to run the awk command right away without the pipe.


### June 10, 2022 ###

	- split bowhead fq files into lanes
	- check qualimap again
	- test haplotypecaller
	
**Splitting bowhead into lanes**

The bowhead .sra finished unpacking overnight. There were the same number of 'spots' written as it listed on the NCBI page, so it appears to have unpacked correctly. I started to split the unzipped fastq files by lane this morning (job: *36094310*).

This should produce 6 files. A fastq for each lane of each read direction. Once these are split, I will gz the files and they will be in the same structure as the right whale samples. Gzip started in a job (*36110001*).

**Marking duplicates ran without issue**

Now I am rerunning qualimap on these two files and beginning to test haplotype caller (but just on the reduced genome file). Qualimap on marked dups started as job *36095208*. 

**Testing Haplotype caller**

I set up a file in scratch/ called scaffolds list. It just lists the names of the 23 scaffolds in the reduced genome file. I believe these are "HiC_scaffold_1", "HiC_scaffold_2" etc. I think this will start an array job that will run haplptype caller on each scaffold independently and I can merge the gvf files after for each individual if this works.

Started Haploytpe caller with the following code as job *36096638*.

```
#!/bin/bash
#SBATCH --job-name=test-haplotype-gvcf
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-23:4
#SBATCH --cpus-per-task=4
#SBATCH --mem=25G
#SBATCH --time=24:00:00

SCAFFOLD=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scaffold_list)
SAMPLE="EGL254-1"

java -Xms20G -Xmx20G -jar ${GATK_JAR} HaplotypeCaller \
  -R reference/Eubalaena_glacialis_HiC_min1Mb.fasta \
  -I alignment/$SAMPLE-merged-nodup_min1Mb.bam \ 
  -O alignment/$SAMPLE.g.vcf.gz \
  -ERC GVCF \
  -L $SCAFFOLD
```

I noticed as soon as it started running that in the array line, as it reads above, it is telling the program to read every 4th element. I meant for it to limit itself to running 4 at a time. Oh well. This will be a good test for now. It is running on scaffolds: 1, 5, 9, 13, 17 and 21. It should have been written as [1-23]%4. To complete the missing scaffolds, I should be able to run the same code as I am currently running, except with --array=[2-4,6-8,10-12,14-16,18-20, 22,23]%4.

There was an error with all of the runs at the same time. It said that a non-standard base was encountered in the reference. There is likely a problem from when the fasta file was subsetted in R into the large scaffolds only.

I started the job as *36099843* on the whole genome file for just a few scaffolds to test variant calling.  Until I realized there was another problem - I was saving to the same output file. I corrected the problem and restarted as job *36110303* using the code below.

```
#!/bin/bash
#SBATCH --job-name=test-haplotype-gvcf
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=[1-23]:6
#SBATCH --cpus-per-task=4
#SBATCH --mem=25G
#SBATCH --time=6:00:00

SCAFFOLD=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scaffold_list)
SAMPLE="EGL254-1"

java -Xms20G -Xmx20G -jar ${GATK_JAR} HaplotypeCaller \
  -R reference/Eubalaena_glacialis_HiC.fasta \
  -I alignment/$SAMPLE-merged-nodup.bam \
  -O alignment/$SAMPLE-$SCAFFOLD-whole.g.vcf.gz \
  -ERC GVCF \
  -L $SCAFFOLD
```

I will run GATKs QCref to check the min1Mb reference file. If there is a problem, I will recreate it, check it, and use that new file in subsequent analyses. There was a problem with the lane endings. I corrected locally in Notepad++ and scp to Cedar in scratch/.

**Module List**

To expedite module loading and ensure all necessary modules are loaded, I saved a module collection that includes the following:

```
module load mugqic/java/openjdk-jdk1.8.0_72 \
mugqic/bvatools/1.6 \
mugqic/fastqc/0.11.5 \
mugqic/trimmomatic/0.36 \
mugqic/samtools/1.9 \
mugqic/bwa/0.7.17 \
mugqic/GenomeAnalysisTK/4.1.0.0 \
mugqic/R_Bioconductor/3.5.0_3.7 \
mugqic/MultiQC/1.10.1 \
mugqic/qualimap/2.2.1 
```

It can be loaded using :
```
module restore rw_mods
```

### June 11, 2022 ###

I am just checking in on the results from yesterday.

**Bowhead**

The bowhead samples finished zipping. 

** Haplotype Caller **

Haplotype caller ran sucessfully I think one a few scaffolds. The run took 3 hours each and used 8.5G on 4 nodes. Still low CPU efficiency.

I finished uploading the min1Mb genome file yesterday and indexing it in the scratch/ and indexed it using bwa. I now started haplotype caller as job *36183948*.

I downloaded one of the whole vcf files to my local computer to check on it. Looks good to me - or at least it looks like a complete vcf file.

** Qualimap **

Qualimap didn't run properly. It didn't finish running on either .bam file. Could it be because the bam file were not sorted again after marking duplicates? No.

** Trimming all ** 

I started a job trimming all samples as job *36207505* after a few kinks were worked out.
I used the following code:
```
#!/bin/bash
#SBATCH --account=def-frasiert
#SBATCH --job-name=trim_sample
#SBATCH --array=1-69%9
#SBATCH --mem-per-cpu=8000M
#SBATCH --cpus-per-task=12
#SBATCH --time=04:00:00

SAMPLELANE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/projects/def-frasiert/RW_WGS/fastq_list)

java -Xmx6G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE \
    -threads 12 -phred33 \
    ~/projects/def-frasiert/RW_WGS/raw_fastq/${SAMPLELANE}_R1_001.fastq.gz \
	~/projects/def-frasiert/RW_WGS/raw_fastq/${SAMPLELANE}_R2_001.fastq.gz \
    ~/scratch/trimmed_reads/$SAMPLELANE-FP.fq.gz \
    ~/scratch/temp/$SAMPLELANE-FUP.fq.gz \
    ~/scratch/trimmed_reads/$SAMPLELANE-RP.fq.gz \
    ~/scratch/temp/$SAMPLELANE-RUP.fq.gz \
    ILLUMINACLIP:~/projects/def-frasiert/RW_WGS/reference/adapters.fa:2:30:15 LEADING:20 SLIDINGWINDOW:5:20 \
    AVGQUAL:30 MINLEN:36 2> ~/projects/def-frasiert/RW_WGS/QC/std_out/$SAMPLELANE.trim.out

fastqc -o ~/projects/def-frasiert/RW_WGS/QC/trimQC/ ~/scratch/trimmed_reads/$SAMPLELANE-FP.fq.gz
fastqc -o ~/projects/def-frasiert/RW_WGS/QC/trimQC/ ~/scratch/trimmed_reads/$SAMPLELANE-RP.fq.gz

rm ~/scratch/temp/$SAMPLELANE-*.fq.gz
```

### June 12, 2022 ###

After 13 hours it looks good. A few runs timed out before the second fastqc could run, so I will double check the slurm output of each to ensure they worked well. And I will manually complete the few that remain unfinished.
I noticed in the output file that trimmomatic couldn't find an absolute directory and needed a relative path. So I updated this and I am rerunning. This should all finish tomorrow evening. In the meantime I am prepping the SRW genome for large scaffolds only.

Hopefully I can start the mapping tomorrow evening.

### June 13, 2022 ###

Today I will:
	- Check Trimmomatic
	- upload .fasta reference sequences with the correct unix line endings
	- index reference genomes for NARW and SRW
	
** Trimmomatic **

Everything appeared to be runnning well and the std_out files indicate that the adapter files were now located without issue. These are running as job *36262239*. The array job will be finished around 8pm. I can then submit the missing/unfinished jobs.

Samples that trimmomatic didn't finish:
	EGL183-1 L002
	EGL183-1 L003
	EGL183-1 L004
	Eau034A L002
	Eau034A L003
	Eau10b L002
	Eau10b L004
	SID179132 L002
	SRR1685383 L002
	SRR1685383 L003
	SRR1685383 L004
	
Resubmitted with 6 hours as job *36262239*.

Finished SID179132 L004 fastqc in an interactive shell.

** Reference sequences **

I opened the reference sequences locally after the small scaffolds were removed and saved them in Notepad++ with UNIX line endings. I then transfered these to the ~/projects/def-frasiert/RW_WGS/reference/. 

I ran the rw_align_genome.sh and srw_align_genome.sh scripts to index the reference files.

** Mapping and Sorting Scripts **

Since the NARW and SRW will be aligned to different genomes, I am mapping them with separate scripts. I subsetting the fastq_list file by species to create 3 files that list sample information by species for input for bwa. Once the trimming is completed, I will submit these.

*NORTH ATLANTIC RIGHT WHALE*

```
#!/bin/bash
#SBATCH --job-name=rw_map_reads
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-36%8
#SBATCH --cpus-per-task=28
#SBATCH --mem=15G
#SBATCH --time=05:00:00

SAMPLELANE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/projects/def-frasiert/RW_WGS/fastq_list)
SAMPLE=$(grep -oP '(?<=\b).*?(?=_2)' <<< "$SAMPLELANE")
LANE=${SAMPLELANE: -4}
READGROUP="@RG\tID:${LANE}\tSM:${SAMPLE}\tLB:RW\tPU:${SAMPLE}\tCN:MCGILL\tPL:ILLUMINA"

bwa mem -M -t 28 \
  -R $READGROUP \
  ~/projects/def-frasiert/RW_WGS/reference/Eubalaena_glacialis_HiC_min1Mb.fasta \
  ~/scratch/trimmed_reads/$SAMPLELANE-FP.fq.gz \
  ~/scratch/trimmed_reads/$SAMPLELANE-RP.fq.gz > \
  ~/scratch/temp_mapping/$SAMPLELANE-aln_min1Mb.sam \
  2> ~/projects/def-frasiert/RW_WGS/QC/std_out/$SAMPLELANE-bwa.err

java -Xmx12G -jar ${GATK_JAR} SortSam \
  -I ~/scratch/temp_mapping/$SAMPLELANE-aln_min1Mb.sam \
  -O ~/scratch/mapped_reads/$SAMPLELANE-aln-sorted_min1Mb.bam \
  -SO coordinate \
  --CREATE_INDEX true --MAX_RECORDS_IN_RAM 500000 \
  2> ~/projects/def-frasiert/RW_WGS/QC/std_out/$SAMPLELANE-SortSam.out

rm ~/scratch/temp_mapping/$SAMPLELANE-*.sam
```

*SOUTHERN RIGHT WHALE*
```
#!/bin/bash
#SBATCH --job-name=srw_map_reads
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-30%8
#SBATCH --cpus-per-task=28
#SBATCH --mem=15G
#SBATCH --time=05:00:00

SAMPLELANE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/projects/def-frasiert/RW_WGS/srw_fastq_list)
SAMPLE=$(grep -oP '(?<=\b).*?(?=_2)' <<< "$SAMPLELANE")
LANE=${SAMPLELANE: -4}
READGROUP="@RG\tID:${LANE}\tSM:${SAMPLE}\tLB:RW\tPU:${SAMPLE}\tCN:MCGILL\tPL:ILLUMINA"

bwa mem -M -t 28 \
  -R $READGROUP \
  ~/projects/def-frasiert/RW_WGS/reference/RWref_HiC_min1Mb.fasta \
  ~/scratch/trimmed_reads/$SAMPLELANE-FP.fq.gz \
  ~/scratch/trimmed_reads/$SAMPLELANE-RP.fq.gz > \
  ~/scratch/temp_mapping/$SAMPLELANE-aln_min1Mb.sam \
  2> ~/projects/def-frasiert/RW_WGS/QC/std_out/$SAMPLELANE-bwa.err

java -Xmx12G -jar ${GATK_JAR} SortSam \
  -I ~/scratch/temp_mapping/$SAMPLELANE-aln_min1Mb.sam \
  -O ~/scratch/mapped_reads/$SAMPLELANE-aln-sorted_min1Mb.bam \
  -SO coordinate \
  --CREATE_INDEX true --MAX_RECORDS_IN_RAM 500000 \
  2> ~/projects/def-frasiert/RW_WGS/QC/std_out/$SAMPLELANE-SortSam.out

rm ~/scratch/temp_mapping/$SAMPLELANE-*.sam
```

*BOWHEAD WHALE*
```
#!/bin/bash
#SBATCH --job-name=bowhead_map_reads
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-3%1
#SBATCH --cpus-per-task=28
#SBATCH --mem=15G
#SBATCH --time=05:00:00

SAMPLELANE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/projects/def-frasiert/RW_WGS/bowhead_fastq_list)
SAMPLE="SRR1685383"
LANE=${SAMPLELANE: -4}B
READGROUP="@RG\tID:${LANE}\tSM:${SAMPLE}\tLB:bowhead\tPU:${SAMPLE}\tCN:BOWHEADGENOMEPROJECT\tPL:ILLUMINA"

bwa mem -M -t 28 \
  -R $READGROUP \
  ~/projects/def-frasiert/RW_WGS/reference/Eubalaena_glacialis_HiC_min1Mb.fasta \
  ~/scratch/trimmed_reads/$SAMPLELANE-FP.fq.gz \
  ~/scratch/trimmed_reads/$SAMPLELANE-RP.fq.gz > \
  ~/scratch/temp_mapping/$SAMPLELANE-aln_min1Mb.sam \
  2> ~/projects/def-frasiert/RW_WGS/QC/std_out/$SAMPLELANE-bwa.err

java -Xmx12G -jar ${GATK_JAR} SortSam \
  -I ~/scratch/temp_mapping/$SAMPLELANE-aln_min1Mb.sam \
  -O ~/scratch/mapped_reads/$SAMPLELANE-aln-sorted_min1Mb.bam \
  -SO coordinate \
  --CREATE_INDEX true --MAX_RECORDS_IN_RAM 500000 \
  2> ~/projects/def-frasiert/RW_WGS/QC/std_out/$SAMPLELANE-SortSam.out

rm ~/scratch/temp_mapping/$SAMPLELANE-*.sam
```

### June 14, 2022 ###

Today I submitted the bwa jobs using the code above as the following job numers:
	- bowhead: *36376703*
	- NARW: *36376704*
	- SRW: *36376705*
	
I started these jobs and the RD info looks great. I thought I these called on the scaffold list, which was not correct, so I killed these jobs right away and restarted them. Then I realized the scaffold list doesn't matter until calling haploytpes. The job numbers above are the correct restarted jobs. The jobs started at ~9:45am. Each job was given 5 hours to run, but I will check on the status at noon. I suspect that BWA will be done at that time and SortSam will be starting on the first samples. 

SortSam didn't have enough memory to run. The jobs are quitting after ~2 hours with an out of memory error. When the array job is finished, BWA MEM will have run successfully and there will be saved .sam files. Unfortunatly, the out of memory error, killed SortSam, but let the script still run, so the .sam files were being deleted. I killed all of the jobs and I am restarting them with the following (more memory).
```
#!/bin/bash
#SBATCH --job-name=bowhead_map_reads
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-3%1
#SBATCH --cpus-per-task=28
#SBATCH --mem=45G
#SBATCH --time=05:00:00

SAMPLELANE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/projects/def-frasiert/RW_WGS/bowhead_fastq_list)
SAMPLE="SRR1685383"
LANE=${SAMPLELANE: -4}B
READGROUP="@RG\tID:${LANE}\tSM:${SAMPLE}\tLB:bowhead\tPU:${SAMPLE}\tCN:BOWHEADGENOMEPROJECT\tPL:ILLUMINA"

bwa mem -M -t 28 \
  -R $READGROUP \
  ~/projects/def-frasiert/RW_WGS/reference/Eubalaena_glacialis_HiC_min1Mb.fasta \
  ~/scratch/trimmed_reads/$SAMPLELANE-FP.fq.gz \
  ~/scratch/trimmed_reads/$SAMPLELANE-RP.fq.gz > \
  ~/scratch/temp_mapping/$SAMPLELANE-aln_min1Mb.sam \
  2> ~/projects/def-frasiert/RW_WGS/QC/std_out/$SAMPLELANE-bwa.err

java -Xmx42G -jar ${GATK_JAR} SortSam \
  -I ~/scratch/temp_mapping/$SAMPLELANE-aln_min1Mb.sam \
  -O ~/scratch/mapped_reads/$SAMPLELANE-aln-sorted_min1Mb.bam \
  -SO coordinate \
  --CREATE_INDEX true --MAX_RECORDS_IN_RAM 500000 \
  2> ~/projects/def-frasiert/RW_WGS/QC/std_out/$SAMPLELANE-SortSam.out
```

Started at 12pm as jobs:
	- bowhead: *36385801*
	- NARW: *36385798*
	- SRW: *36385799*
	
This was taking too much memory. Everything was crashing willy-nilly. While a few jobs were completing, it was nearly impossible to easily confirm which ones worked and which didn't and why. I decided to start over and try two things.

First, I am going to try the following code to attempt to not save the .sam file at all by using samtools.
```
bwa mem -M -t 28 \
  -R $READGROUP \
  ~/projects/def-frasiert/RW_WGS/reference/Eubalaena_glacialis_HiC_min1Mb.fasta \
  ~/scratch/trimmed_reads/$SAMPLELANE-FP.fq.gz \
  ~/scratch/trimmed_reads/$SAMPLELANE-RP.fq.gz | \
samtools sort -o ~/scratch/mapped_reads/$SAMPLELANE-aln-sorted_min1Mb.bam

samtools index ~/scratch/mapped_reads/$SAMPLELANE-aln-sorted_min1Mb.bam
```

This didn't work so I am going to try the orinal code in two parts (BWA first, then SortSam, then delete .sam), but I will only run 1 species at a time and fewer parallel jobs and delete the sam files as they are finished.

I am running an array on 9 jobs (first 3 files - EGL00252, EGL013 and EGL140).
I started the bwa job as *36400930* (started at ~3:15pm). When it is done, I will run the SortSam job on the same 9 files and then delete the .sam files. I will iterate through this process until all 69 files are mapped.

Progress will be completed in this table as the jobs run.

| Purpose | Which Samples          | Job ID   |
|---------|------------------------|----------|
| bwa mem | EGL00252 EGL013 EGL140 | 36400930 |
| sortsam | EGL00252 EGL013 EGL140 | 36411550 |
| bwa mem | EGL183 EGL254 EGL272   | 36510274 |
| sortsam | EGL183 EGL254 EGL272   | 36523230 | EGL183 oom
| bwa mem | SRW SAMPLES            | 36440905 | Not all finished
| sortsam | SRW SAMPLES            | 36440905 | Not all finished
| sortsam | 8 SRW samples          | 36509793 | 4 oom failures
| bwa mem | 18 SRW samples         | 36516565 |
| sortsam | 18 SRW samples         | 36516565 |
| bwa mem | EGL276 EGL308 EGL312   | 36528642 |
| sortsam | EGL276 EGL308 EGL312   | 36537292 |

## June 15, 2022 ##

I am having memory problems with BWA and SortSam. Most things seem to be working, but it is slower than expected and some runs have timed out, bwa seemed to unload itself at a certain point and memory had been exceeded. 

I am slowing down the process and going step by step to ensure it gets run. 
First thing this morning, I am completing the SortSam for the files where the .sam was created, but not the .bam. (8 SRW files)

I also restarted the mapping only on the 9 NARW files (10-18) I had run last night and accidentally overwrote. I know this step will only take 2-3 hours. I also ran the sortsam on these files.

Progress:
	- EGL336, SID179132, SID181803 and EGL183 are the only NARW left and SortSam is running now as jobs: *36550595* and *36551310*.
	- Eau19 and Eau10b need to run srw_map_min1Mb.sh with more memory (started at  as *36558794*) and tmp_dir.
	- SRR (bowhead) needs to run in full with more memory (started at  as *36558760*) and tmp_dir.
	- EGL183 still to be sortsam using --TMP_DIR


**Merging .bam Script**

The output .bam files have the following naming convention:
```~/scratch/mapped_reads/$SAMPLELANE-aln-sorted_min1Mb.bam```
where
```$SAMPLELANE = Eau034A_2-2253009_S20_L004``` 

I am going to use the following code to merge all 3 lanes for each of the NARW and SRW samples where bam_list is a list of the samplelane info with the lane number dropped.

```
#!/bin/bash
#SBATCH --job-name=merge_bam
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-22%12
#SBATCH --cpus-per-task=20
#SBATCH --mem=10G
#SBATCH --time=03:00:00

BAMNAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/projects/def-frasiert/RW_WGS/bam_list)
SAMPLE=$(grep -oP '(?<=\b).*?(?=_2)' <<< "$BAMNAME")

java -Xmx8G -jar ${GATK_JAR} MergeSamFiles \
  --INPUT ~/scratch/mapped_reads/${BAMNAME}_L002-aln-sorted_min1Mb.bam \
  --INPUT ~/scratch/mapped_reads/${BAMNAME}_L003-aln-sorted_min1Mb.bam \
  --INPUT ~/scratch/mapped_reads/${BAMNAME}_L004-aln-sorted_min1Mb.bam \
  --OUTPUT ~/projects/def-frasiert/RW_WGS/merged_reads/$SAMPLE-merged_min1Mb.bam \
  --CREATE_INDEX true \
  2> ~/projects/def-frasiert/RW_WGS/QC/std_out/$SAMPLE-merge.out
```
 
 ** Scaffold List for GATK **
 
 I created two files in R and saved to to the RW_WGS/ called rw_sample_scaffolds and srw_sample_scaffolds.
 
 These files contain 2 tab-separated columns. The first is just sample names, the second is the scaffold list. There is a separate row for each combination of sample name and scaffold.
 
 I can then use a job array to call a specific line in this file and then call either scaffold for GATK or sample name for the output. 
 
### June 16, 2022 ###

Today I would like to do the following:
	- Finished mapping and sorting
	- Merge all bam files
	- Mark Duplicates
	- Start Qualimap

** Mapping and Sorting **

Writing to the temp directory takes more time, but there isn't a bottleneck at the memory, so it is working well enough. I have 5 files left to finish that are all running under the following job IDs:
    - *36586738_10* EGL183 L002 SortSam - Timed out
    - *36586738_11* EGL183 L003 SortSam - Timed out
    - *36586738_12* EGL183 L004 SortSam - Timed out
    - *36587056_3* Bowhead L006 SortSam
	- *36587560_19* Eau10 L002 SortSam
	- *36601885_10* EGL183 L002 SortSam - TIMED OUT! Soooo close
    - *36601885_11* EGL183 L003 SortSam
    - *36601885_12* EGL183 L004 SortSam

Restarting final SortSam. EGL183-1 L002 with 10 hours. *36636711*

** Bowhead Merging and Marking Dups **

I am testing the merging and marking dups script on the bowhead sample using the following code, it was started as job *36599567*.
```
#!/bin/bash
#SBATCH --job-name=bowhead_mergemark_bam
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=20
#SBATCH --mem=30G
#SBATCH --time=06:00:00

BAMNAME="SRR1685383"

java -Xmx25G -jar ${GATK_JAR} MergeSamFiles \
  --INPUT ~/scratch/mapped_reads/${BAMNAME}.L003-aln-sorted_min1Mb.bam \
  --INPUT ~/scratch/mapped_reads/${BAMNAME}.L004-aln-sorted_min1Mb.bam \
  --INPUT ~/scratch/mapped_reads/${BAMNAME}.L006-aln-sorted_min1Mb.bam \
  --OUTPUT ~/scratch/merged_reads/$BAMNAME-merged_min1Mb.bam \
  2> ~/projects/def-frasiert/RW_WGS/QC/std_out/$BAMNAME-merge.out
  
java -Xmx25G -jar ${GATK_JAR} MarkDuplicates \
  --REMOVE_DUPLICATES false --CREATE_INDEX true \
  --INPUT ~/scratch/merged_reads/$BAMNAME-merged_min1Mb.bam \
  --OUTPUT ~/projects/def-frasiert/RW_WGS/merged_bam/$BAMNAME-merged_marked.bam \
  --TMP_DIR ~/scratch/temp_marking/ \
  --METRICS_FILE=~/projects/def-frasiert/RW_WGS/QC/mergeQC/$BAMNAME-merged_marked.metrics
 ``` 

Marking Duplicates required a relative path instead of an absolute path for the metrics file. This job quit after ~85mins after it finished merging bam files. I restarted the code for marking duplicates alone as job *36602452*. This marking only was very close, but timed out. I restarted it with 8hours as job *36633920*.

I started the rw_merge_mark.sh on all of the RW samples with the exception of EGL183 that is still running through SortSam. This started as job *36606967*. It is an array job given 6 hours for 21 samples, with a limit of 12 at a time. It could run until 4am. Forgot to change the job name from bowhead.


## June 17, 2022 ##

So far, bam files have been successfully created for all samples and all the bam files have been merged and duplicates marked for all but 1 sample (EGL183-1 will finish later this morning). 

Next Steps:
	- QUALIMAP needs to be run on all of the samples when EGL183-1 is finished
	- Unmerged .bam files can be erased
	- GATK HaplotypeCaller to run on SRW
	- GATK HaplotypeCaller to run on NARW
	- GATK HaplotypeCaller to run on bowhead

** HaplotypeCaller on SRW **

Normally, I would like to run this on NARW first, but since EGL183 is still running, I will start with the SRW. I have created a file called srw_sample_scaffolds in the RW_WGS/ which contains a row for each sample name and scaffold combination. I will run the following script as an array job that will run 21 jobs at a time (meaning it will run all scaffolds for each sample at the same time). 

I had to reindex the reference sequence before starting. SRW haplotypecaller started as job *36685189*.

```
#!/bin/bash
#SBATCH --job-name=srw-haplotype-gvcf
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-210%21
#SBATCH --cpus-per-task=4
#SBATCH --mem=25G
#SBATCH --time=6:00:00

FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" srw_sample_scaffolds)

SAMPLE=$(echo $FILE | awk '{print $1}')
SCAFFOLD=$(echo $FILE | awk '{print $2}')

java -Xms20G -Xmx20G -jar ${GATK_JAR} HaplotypeCaller \
  -R ~/projects/def-frasiert/RW_WGS/reference/RWref_HiC_min1Mb.fasta \
  -I ~/projects/def-frasiert/RW_WGS/merged_bam/$SAMPLE-merged_marked.bam \ 
  -O ~/scratch/gvcf/SRW/$SCAFFOLD/$SAMPLE-$SCAFFOLD.g.vcf.gz \
  -ERC GVCF \
  -L $SCAFFOLD \
  2> ~/projects/def-frasiert/RW_WGS/QC/gvcfQC/$SAMPLE-$SCAFFOLD-hapcall-gvcf.out
```

** Qualimap **

First I need to create the qualimap input file I will call *qualimap_sample_file*.
I will use the following line of code to generate this list and edit the Species Name manually.

```  \ls -1 *.bam | sed 's,\(^.*\)\(-merged.*\),\1\t\1\2\tNARW,' > qualimap_sample_file ```

I will then run the following code to run Qualimap from inside the RW_WGS/merged_bam/:

```
#!/bin/bash
#SBATCH --job-name=qualimap_testing
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=30
#SBATCH --mem=25G
#SBATCH --time=24:00:00

qualimap multi-bamqc -d qualimap_sample_file \
    -outdir ~/projects/def-frasiert/QC/mapQC -nr 10000 \
    --java-mem-size=20G -r
```

This job was started as job *36693239*. I monitored the runtime and felt that some of the larger scaffolds were going to time out. I cancelled the job for the large scaffolds that were running, and all pending jobs and restarted them with 9 hours instead of 6. (The following jobs were restarted as *36701640*: 1,2,6,9,22,23,27-210).

At the end of the day, I have HaplotypeCaller running on all of the SRW samples and Qualimap running across all samples.

Tomorrow evening, I will check on the progress of everything. 
The next steps are:
	- Run GenomicsDBImport on Scaffold 6 of SRW gvcfs
	- Inspect Qualimap
	- Run GenotypeGVCFs in Scaffold 6 of SRW gvcfs
	- Run HaplotypeCaller on NARW
	- Run GenomicsDBImport on all SRW scaffolds
	
## June 18, 2022 ##

I checked everything this morning. Haplotypecaller was running well. It is on track to finish Sunday morning. 

Qualimap was taking longer than expected based on previous tests. I cancelled the job and started another with a 2 day time limit (*36779472*).

## June 19, 2022 ##

Starting GenomicsDBImport on scaffold 1 of SRW. First I need to create the sample file with the following code run in the scaffold folder.
```
 \ls -1 *.g.vcf.gz | sed 's,\(^.*\)\(-HiC.*\),\1\t\1\2,' > HiC_scaffold_1_SRW_samples
 ```
 Here is the script I used.
 ```
 #!/bin/bash
#SBATCH --job-name=srw_genomicsdb
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=6
#SBATCH --mem=25G
#SBATCH --time=48:00:00

java -Xms20G -Xmx20G -jar ${GATK_JAR} GenomicsDBImport \
    --genomicsdb-workspace-path ~/scratch/genomicsDB/SRW/ \
    -L HiC_scaffold_1 \
    --sample-name-map ~/scratch/gvcf/SRW/HiC_scaffold_1/HiC_scaffold_1_SRW_samples \
    --tmp-dir temp \
    --reader-threads 5
 ```
Note this neede to be run from the scaffold folder as absolute paths were not being accepted. Also, The workspace path needs to be a nonexistant directory. I had en empty directory and it wasn't working.

In the future, this directory I called SRW, should be more specific and refer to the scaffold name.

Not all of the files finished for haplotype caller, tomorrow I will check and sort output files.

## June 20, 2022 ##

Today I am:
	- Starting GenotypeGVCFs
	- Completing the SRW Haplotype Caller
	- Running GenomicsDBImport on remaining SRW scaffolds
	- Starting HaplotypeCaller for RW
	
**GenotypeGVCFs**

The GenomicsDBImport finished running for scaffold 1 in just under 4 hours.

Now I am starting the GenotypeGVCFs on scaffold 1 for SRW using the following code as job *36907804*:

```
#!/bin/bash
#SBATCH --job-name=srw_genotype_gvcf1
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=6
#SBATCH --mem=25G
#SBATCH --time=24:00:00

java -Xms20G -Xmx20G -jar ${GATK_JAR} GenotypeGVCFs \
   -R ~/projects/def-frasiert/RW_WGS/reference/RWref_HiC_min1Mb.fasta \
   -V gendb://SRW/ \
   -O ~/projects/def-frasiert/RW_WGS/vcf/SRW/srw_HiC_scaffold_1.vcf.gz
```

The code was from from the parent folder for the database (scratch/genomicsDB/). For the other scaffolds, I would also like to add a std error file to be saved in QC/.

This completed in about 8 hours. When I run the other SRW scaffolds, I will give them 10 hours to finish even though this was one of the largest. 

I created a small file that called the database folder names and the scaffold numbers. I then used it to run the following on the remaining genomicsDB files from the SRW_workspace/ as job *36939317*.

```
#!/bin/bash
#SBATCH --job-name=srw_genotype_gvcf
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=3-5,7,8,10-14,16-21
#SBATCH --cpus-per-task=6
#SBATCH --mem=25G
#SBATCH --time=10:00:00

FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" srw_genDB_inputs)

FOLDER=$(echo $FILE | awk '{print $1}')
SCAFFOLD=$(echo $FILE | awk '{print $2}')

java -Xms20G -Xmx20G -jar ${GATK_JAR} GenotypeGVCFs \
   -R ~/projects/def-frasiert/RW_WGS/reference/RWref_HiC_min1Mb.fasta \
   -V gendb://$FOLDER/ \
   -O ~/projects/def-frasiert/RW_WGS/vcf/SRW/srw_HiC_scaffold_$SCAFFOLD.vcf.gz
```

Not all genomicsDBimport runs were finished. Still need to run this script on scaffold 2, 6, 9 and 15.

**SRW HaplotypeCaller**

From the log files, I identified that 4 files did not run to completion with haplotype caller. The slurm output for these files had a a time out error written in them and they were not empty. By looking at the time these were written, I could identify which Sample_Scaffold they belonged to and I could check the log file in the gvcfQC/ to confirm they were not finished. I started HaploytpeCaller for the following with 12 hours to ensure they finish running. This started as job: *36910923*.

**SRW GenomicsDBImport**

Because the sample_maps needed relative paths, I think I need to run a script from each of the gvcf/scaffold/ directories. I created scripts for each scaffold and saved them in a directory within the scripts directory.

Sample map files were also created and housed in each of the gvcf/scaffold/.

Jobs:
SRW Scaffold 3 - *36918350*
SRW Scaffold 4 - *36918854*
SRW Scaffold 5 - *36918855*
SRW Scaffold 7 - *36918916*
SRW Scaffold 8 - *36918932*
SRW Scaffold 10 - *36919123*
SRW Scaffold 11 - *16919128*
SRW Scaffold 12 - *36919141*
SRW Scaffold 13 - *36919184*
SRW Scaffold 14 - *36926999*
SRW Scaffold 15 - *36927166*
SRW Scaffold 16 - *36927169*
SRW Scaffold 17 - *36927219*
SRW Scaffold 18 - *36927381*
SRW Scaffold 19 - *36927649*
SRW Scaffold 20 - *36927682*
SRW Scaffold 21 - *36927683*
SRW Scaffold 2 - *36940498*
SRW Scaffold 6 - *36939720*
SRW Scaffold 9 - *36936998*

**RW HaplotypeCaller**

I started HaplotypeCaller on the NARW samples as job *36929083*.
```
#!/bin/bash
#SBATCH --job-name=rw-haplotype-gvcf
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-276%23
#SBATCH --cpus-per-task=4
#SBATCH --mem=25G
#SBATCH --time=10:00:00

FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" rw_sample_scaffolds)

SAMPLE=$(echo $FILE | awk '{print $1}')
SCAFFOLD=$(echo $FILE | awk '{print $2}')

java -Xms20G -Xmx20G -jar ${GATK_JAR} HaplotypeCaller \
  -R ~/projects/def-frasiert/RW_WGS/reference/Eubalaena_glacialis_HiC_min1Mb.fasta \
  -I ~/projects/def-frasiert/RW_WGS/merged_bam/$SAMPLE-merged_marked.bam \
  -O ~/scratch/gvcf/NARW/$SCAFFOLD/$SAMPLE-$SCAFFOLD.g.vcf.gz \
  -ERC GVCF \
  -L $SCAFFOLD \
  2> ~/projects/def-frasiert/RW_WGS/QC/gvcfQC/$SAMPLE-$SCAFFOLD-hapcall-gvcf.out
```

### June 21, 2022 ###

**GenotypeGVCFs**

Today I went to start the GenotypeGVCFs on the last 4 scaffolds for the SRW (2,6,9 and 15) as job *36996044*. Scaff15 didn't run for more than a few seconds and when I checked the genomicsDB output, there seemed to have been a problem. I started the GenomicsBDImport script again for scaffold 15 with a longer time limit (job *36996714*).

When this finished, I started the GenotypeGVCFs script for scaff 15 as job *37587217* on June 26.

**NARW HaplotypeCaller**

HaplotypeCaller is progressing as expected for the NARW data. By 9:45am, approx 1/2 of the files have finished running. I will check on the progress this evening. 

### June 22, 2022 ###

**GenomicsBDImport RW**

Today I will be running GenomicsDBImport of the NARW scaffolds. 

Starting GenomicsDBImport on scaffold 1 of SRW. First I need to create the sample file with the following code run in the scaffold folder.
```
 \ls -1 *.g.vcf.gz | sed 's,\(^.*\)\(-HiC.*\),\1\t\1\2,' > HiC_scaffold_1_NARW_samples
 ```
 Here is an example of the script I used.
 ```
#!/bin/bash
#SBATCH --job-name=narw_genomicsdb_1
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=6
#SBATCH --mem=25G
#SBATCH --time=6:00:00

java -Xms20G -Xmx20G -jar ${GATK_JAR} GenomicsDBImport \
    --genomicsdb-workspace-path ~/scratch/genomicsDB/NARW_workspace/scaff1/ \
    -L HiC_scaffold_1 \
    --sample-name-map ~/scratch/gvcf/NARW/HiC_scaffold_1/HiC_scaffold_1_NARW_samples \
    --tmp-dir ~/scratch/temp \
    --reader-threads 5
 ```
Note this needed to be run from the scaffold folder as absolute paths were not being accepted. Also, The workspace path needs to be a nonexistant directory. I had en empty directory and it wasn't working.

Started as jobs:
NARW Scaffold 1 - *37088572* - Didn't finish
NARW Scaffold 2 - *37088870*
NARW Scaffold 3 - *37088871*
NARW Scaffold 4 - *37088872*
NARW Scaffold 5 - *37088873*
NARW Scaffold 6 - *37088875*
NARW Scaffold 7 - *37089373*
NARW Scaffold 8 - *37088879* - Didn't finish
NARW Scaffold 9 - *37088880*
NARW Scaffold 10 - *37089085*
NARW Scaffold 11 - *37089087*
NARW Scaffold 12 - *37089100*
NARW Scaffold 13 - *37089104*
NARW Scaffold 14 - *37089199*
NARW Scaffold 15 - *37089202*
NARW Scaffold 16 - *37089205* - Didn't finish
NARW Scaffold 17 - *37089209*
NARW Scaffold 18 - *37089210*
NARW Scaffold 19 - *37089212*
NARW Scaffold 20 - *37089225* - Didn't finish
NARW Scaffold 21 - *37089245*
NARW Scaffold 22 - *37089275*
NARW Scaffold 111 - *37089299*

For the jobs that didn't finish, I erased the scaff/ directories and started again with 10 hours. 

NARW Scaffold 1 *37127835*
NARW Scaffold 8 *37127839*
NARW Scaffold 16 *37127841*
NARW Scaffold 20 *37127842*


### June 23, 2022 ###

**GenotypeGVCFs**

I created a small file that called the database folder names and the scaffold numbers. I then used it to run the following on the genomicsDB files from the NARW_workspace/ as job *37187015*.

```
#!/bin/bash
#SBATCH --job-name=narw_genotype_gvcf
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-23
#SBATCH --cpus-per-task=6
#SBATCH --mem=25G
#SBATCH --time=10:00:00

FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" narw_genDB_inputs)

FOLDER=$(echo $FILE | awk '{print $1}')
SCAFFOLD=$(echo $FILE | awk '{print $2}')

java -Xms20G -Xmx20G -jar ${GATK_JAR} GenotypeGVCFs \
   -R ~/projects/def-frasiert/RW_WGS/reference/Eubalaena_glacialis_HiC_min1Mb.fasta \
   -V gendb://$FOLDER/ \
   -O ~/projects/def-frasiert/RW_WGS/vcf/NARW/narw_HiC_scaffold_$SCAFFOLD.vcf.gz \
   2> ~/projects/def-frasiert/RW_WGS/QC/gvcfQC/NARW_$SCAFFOLD-gvcf.out
```

**Testing Filtration Parameters**

In order to apply the base recalibration, we need to select a set of variants we are confident in to train the model. I have come up with two sets of filters ('tight' and 'loose') that both will produce variants I am confident with, but we will assess the convergence of the models on these two datasets.

I am testing this on two contigs of SRW - HiC_scaffold_1 and HiC_scaffold_8.

Tight Filters:
Number of Alternate Alleles = 1
Depth of all samples > 15
Mapping Quality > 40
Genotype Quality of all samples > 40
Data missing from >3 samples

Loose Filters:
Number of Alternate Alleles <=3
Depth of all samples > 10
Mapping Quality > 30
Genotype Quality of all samples > 30
Data missing from <= 5 samples

```
bcftools filter -i 'N_ALT=1 && FORMAT/DP[0-9]>15 && INFO/MQ>40 && FORMAT/GQ[0-9]>40 && N_MISSING<3' -Oz -o ~/scratch/temp/tight_filtered_srw_HiC_scaffold_1.vcf.gz ~/projects/def-frasiert/RW_WGS/vcf/SRW/srw_HiC_scaffold_1.vcf.gz

bcftools filter -i 'N_ALT=1 && FORMAT/DP[0-9]>15 && INFO/MQ>40 && FORMAT/GQ[0-9]>40 && N_MISSING<3' -Oz -o ~/scratch/temp/tight_filtered_srw_HiC_scaffold_8.vcf.gz ~/projects/def-frasiert/RW_WGS/vcf/SRW/srw_HiC_scaffold_8.vcf.gz

bcftools filter -i 'N_ALT<=3 && FORMAT/DP[0-9]>10 && INFO/MQ>30 && FORMAT/GQ[0-9]>30 && N_MISSING<=5' -Oz -o ~/scratch/temp/loose_filtered_srw_HiC_scaffold_1.vcf.gz ~/projects/def-frasiert/RW_WGS/vcf/SRW/srw_HiC_scaffold_1.vcf.gz

bcftools filter -i 'N_ALT<=3 && FORMAT/DP[0-9]>10 && INFO/MQ>30 && FORMAT/GQ[0-9]>30 && N_MISSING<=5' -Oz -o ~/scratch/temp/loose_filtered_srw_HiC_scaffold_8.vcf.gz ~/projects/def-frasiert/RW_WGS/vcf/SRW/srw_HiC_scaffold_8.vcf.gz
```
First I need to index the .vcf.gz input files.
```
java -Xmx20G -jar ${GATK_JAR} IndexFeatureFile \
	-F ~/scratch/temp/loose_filtered_srw_HiC_scaffold_8.vcf.gz
```

I will apply the base recalibrator using these variant sets with the following code for both loose and tight sets on 4 samples only.

```
#!/bin/bash
#SBATCH --job-name=tight_srw_baserecalibrate
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-4
#SBATCH --cpus-per-task=4
#SBATCH --mem=25G
#SBATCH --time=6:00:00

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" srw_samples)

java -Xmx20G -jar ${GATK_JAR} BaseRecalibrator \
  -R ~/projects/def-frasiert/RW_WGS/reference/RWref_HiC_min1Mb.fasta \
  --known-sites ~/scratch/temp/tight_filtered_srw_HiC_scaffold_1.vcf.gz \
  -L HiC_scaffold_1 \
  -O ~/scratch/temp/recalibration/tight_$SAMPLE_HiC_scaffold_1_recalibration_report.grp \
  -I ~/projects/def-frasiert/RW_WGS/merged_bam/$SAMPLE-merged_marked.bam
  
java -Xmx20G -jar ${GATK_JAR} BaseRecalibrator \
  -R ~/projects/def-frasiert/RW_WGS/reference/RWref_HiC_min1Mb.fasta \
  --known-sites ~/scratch/temp/tight_filtered_srw_HiC_scaffold_8.vcf.gz \
  -L HiC_scaffold_8 \
  -O ~/scratch/temp/recalibration/tight_$SAMPLE_HiC_scaffold_8_recalibration_report.grp \
  -I ~/projects/def-frasiert/RW_WGS/merged_bam/$SAMPLE-merged_marked.bam 
```

I then have to apply the BQSR to the bam files. 

```
#!/bin/bash
#SBATCH --job-name=tight_srw_baserecalibrate
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-4
#SBATCH --cpus-per-task=4
#SBATCH --mem=25G
#SBATCH --time=6:00:00

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" srw_samples)

java -Xmx20G -jar ${GATK_JAR} ApplyBQSR \
  -R ~/projects/def-frasiert/RW_WGS/reference/RWref_HiC_min1Mb.fasta \
  -L HiC_scaffold_1 \
  -bqsr ~/scratch/temp/recalibration/tight_${SAMPLE}_HiC_scaffold_1_recalibration_report.grp \
  -O ~/scratch/temp/recalibration/tight_${SAMPLE}_HiC_scaffold_1_recal.bam \
  -I ~/projects/def-frasiert/RW_WGS/merged_bam/$SAMPLE-merged_marked.bam
  
java -Xmx20G -jar ${GATK_JAR} ApplyBQSR \
  -R ~/projects/def-frasiert/RW_WGS/reference/RWref_HiC_min1Mb.fasta \
  -L HiC_scaffold_1 \
  -bqsr ~/scratch/temp/recalibration/tight_${SAMPLE}_HiC_scaffold_8_recalibration_report.grp \
  -O ~/scratch/temp/recalibration/tight_${SAMPLE}_HiC_scaffold_8_recal.bam \
  -I ~/projects/def-frasiert/RW_WGS/merged_bam/$SAMPLE-merged_marked.bam
```

Recalibrate the new bam with the known-sites vcf
```
#!/bin/bash
#SBATCH --job-name=tight_srw_baserecalibrate
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-4
#SBATCH --cpus-per-task=4
#SBATCH --mem=25G
#SBATCH --time=1:30:00

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" srw_samples)

java -Xmx20G -jar ${GATK_JAR} BaseRecalibrator \
  -R ~/projects/def-frasiert/RW_WGS/reference/RWref_HiC_min1Mb.fasta \
  --known-sites ~/scratch/temp/tight_filtered_srw_HiC_scaffold_1.vcf.gz \
  -L HiC_scaffold_1 \
  -O ~/scratch/temp/recalibration/after_tight_${SAMPLE}_HiC_scaffold_1_recalibration_report.grp \
  -I ~/scratch/temp/recalibration/tight_${SAMPLE}_HiC_scaffold_1_recal.bam
  
java -Xmx20G -jar ${GATK_JAR} BaseRecalibrator \
  -R ~/projects/def-frasiert/RW_WGS/reference/RWref_HiC_min1Mb.fasta \
  --known-sites ~/scratch/temp/tight_filtered_srw_HiC_scaffold_8.vcf.gz \
  -L HiC_scaffold_8 \
  -O ~/scratch/temp/recalibration/after_tight_${SAMPLE}_HiC_scaffold_8_recalibration_report.grp \
  -I ~/scratch/temp/recalibration/tight_${SAMPLE}_HiC_scaffold_8_recal.bam
```
Create plots for before and after recalibration.
```
#!/bin/bash
#SBATCH --job-name=Covariate_plots
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-4
#SBATCH --cpus-per-task=4
#SBATCH --mem=25G
#SBATCH --time=1:30:00

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" srw_samples)

java -Xms7G -Xmx7G -jar ${GATK_JAR}  AnalyzeCovariates \
     -before ~/scratch/temp/recalibration/tight_${SAMPLE}_HiC_scaffold_8_recalibration_report.grp \
	 -after ~/scratch/temp/recalibration/after_tight_${SAMPLE}_HiC_scaffold_8_recalibration_report.grp\
     -plots Tight_8_AnalyzeCovariates-$SAMPLE.pdf
	 
java -Xms7G -Xmx7G -jar ${GATK_JAR}  AnalyzeCovariates \
     -before ~/scratch/temp/recalibration/tight_${SAMPLE}_HiC_scaffold_1_recalibration_report.grp \
	 -after ~/scratch/temp/recalibration/after_tight_${SAMPLE}_HiC_scaffold_1_recalibration_report.grp\
     -plots Tight_1_AnalyzeCovariates-$SAMPLE.pdf
	 
java -Xms7G -Xmx7G -jar ${GATK_JAR}  AnalyzeCovariates \
     -before ~/scratch/temp/recalibration/loose_${SAMPLE}_HiC_scaffold_8_recalibration_report.grp \
	 -after ~/scratch/temp/recalibration/after_loose_${SAMPLE}_HiC_scaffold_8_recalibration_report.grp\
     -plots Loose_8_AnalyzeCovariates-$SAMPLE.pdf
	 
java -Xms7G -Xmx7G -jar ${GATK_JAR}  AnalyzeCovariates \
     -before ~/scratch/temp/recalibration/loose_${SAMPLE}_HiC_scaffold_1_recalibration_report.grp \
	 -after ~/scratch/temp/recalibration/after_loose_${SAMPLE}_HiC_scaffold_1_recalibration_report.grp\
     -plots Loose_1_AnalyzeCovariates-$SAMPLE.pdf
```
I saved the plots to my local computer and there seemed to be a problem on the 'after' for scaffold 8. But the results looked better for scaffold 1. The tight and the loose filters performed very similarly.

I saved these plots locally.

Based on this I would feel more comfortable using the looser filters.

All this however, gives me question about how it works with repsect to read group. The original and recalled .bam files all appear to have the correct RG in the header, but it doesn't appear to be used properly by BaseRecalibrator as it identifies the SM instead of ID as the ReadGroup.

**Testing pca on Scaffold 1 for NARW**

I will be going through this a bit more tomorrow, but just to see if I am doing this cirrectly, I am using the multisample vcf for NARW scaffold 1 to conduct a PCA using plink and plotted in R.

```
plink --vcf narw_HiC_scaffold_1.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 \
--out narw_scaff1

plink --vcf narw_HiC_scaffold_1.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--extract narw_scaff1.prune.in \
--make-bed --pca --out narw_scaff1
```

In R:
```
setwd("C:/Users/c_cro/Documents/PhD/RightWhale/WGS/test_recal_filters/pca/")
library(tidyverse)
library(ggplot2)


pca<-read.table("narw_scaff1.eigenvec")
eigenval<-scan("narw_scaff1.eigenval")
pca<-pca[,-1]
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC",1:(ncol(pca)-1))
pve <- data.frame(PC = 1:12, pve = eigenval/sum(eigenval)*100)

a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained")

# plot pca
b <- ggplot(pca, aes(PC1, PC2)) + geom_point(size = 3)
b <- b + scale_colour_manual(values = c("red", "blue"))
b <- b + coord_equal()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
b
```

### June 24, 2022 ###

Today I checked and 6 scaffolds did not finish running genotype gvcfs. I restarted them with a 24hr time limit to ensure they finish. This started as job *37300552*.

**PCA plots**

I tested plink and R to make a few more PCA plots to see what a few other scaffolds look like. This runs super quickly on an interactive node.
 
```
plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/NARW/narw_HiC_scaffold_2.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 \
--out narw_scaff2

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/NARW/narw_HiC_scaffold_2.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--extract narw_scaff2.prune.in \
--make-bed --pca --out narw_scaff2

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/NARW/narw_HiC_scaffold_3.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 \
--out narw_scaff3

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/NARW/narw_HiC_scaffold_3.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--extract narw_scaff3.prune.in \
--make-bed --pca --out narw_scaff3
```

### June 26, 2022 ###

After speaking with Michael, we need the largest set of variants possible to train the recalibration model. Therefore, instead of just using 1 scaffold, we will use the whole genome.

At first I tried to merge the vcf files for each scaffold using the following:
```
vcf-merge srw_HiC_scaffold_1.vcf.gz srw_HiC_scaffold_2.vcf.gz srw_HiC_scaffold_3.vcf.gz srw_HiC_scaffold_4.vcf.gz srw_HiC_scaffold_5.vcf.gz srw_HiC_scaffold_6.vcf.gz srw_HiC_scaffold_7.vcf.gz srw_HiC_scaffold_8.vcf.gz srw_HiC_scaffold_9.vcf.gz srw_HiC_scaffold_10.vcf.gz srw_HiC_scaffold_11.vcf.gz srw_HiC_scaffold_12.vcf.gz srw_HiC_scaffold_13.vcf.gz srw_HiC_scaffold_14.vcf.gz srw_HiC_scaffold_15.vcf.gz srw_HiC_scaffold_16.vcf.gz srw_HiC_scaffold_17.vcf.gz srw_HiC_scaffold_18.vcf.gz srw_HiC_scaffold_19.vcf.gz srw_HiC_scaffold_20.vcf.gz srw_HiC_scaffold_21.vcf.gz | bgzip -c > merged_srw_all_scaffold.vcf.gz
```

The code was correct, but the job didn't finish and flagged that srw_scaffold_15 had an index older than the actual file. I will recreate this file to ensure it is correct. Also of note, the merged VCF did not create the files in the structure we would want. I think I would need to concatenate them. 

I will fix scaffold 15, and then try concatenating them using bcftools concat.
```
bcftools concat srw_HiC_scaffold_1.vcf.gz srw_HiC_scaffold_2.vcf.gz srw_HiC_scaffold_3.vcf.gz srw_HiC_scaffold_4.vcf.gz srw_HiC_scaffold_5.vcf.gz srw_HiC_scaffold_6.vcf.gz srw_HiC_scaffold_7.vcf.gz srw_HiC_scaffold_8.vcf.gz srw_HiC_scaffold_9.vcf.gz srw_HiC_scaffold_10.vcf.gz srw_HiC_scaffold_11.vcf.gz srw_HiC_scaffold_12.vcf.gz srw_HiC_scaffold_13.vcf.gz srw_HiC_scaffold_14.vcf.gz srw_HiC_scaffold_15.vcf.gz srw_HiC_scaffold_16.vcf.gz srw_HiC_scaffold_17.vcf.gz srw_HiC_scaffold_18.vcf.gz srw_HiC_scaffold_19.vcf.gz srw_HiC_scaffold_20.vcf.gz srw_HiC_scaffold_21.vcf.gz -Oz -o merged_srw_all_scaffold.vcf.gz

java -Xmx20G -jar ${GATK_JAR} IndexFeatureFile \
	-F ~/projects/def-frasiert/RW_WGS/vcf/SRW/merged_srw_all_scaffold.vcf.gz
```

First, I need to call variants on this merged haplotype file. This will be done in two steps 1) GATK's variant filtration which will mark variants as pass or fail and 2) SelectVariants.

```
#!/bin/bash
#SBATCH --job-name=srw_variant_filtration
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=10
#SBATCH --mem=25G
#SBATCH --time=12:00:00

java -Xmx20G -jar ${GATK_JAR} VariantFiltration \
   -R ~/projects/def-frasiert/RW_WGS/reference/RWref_HiC_min1Mb.fasta \
   -V ~/projects/def-frasiert/RW_WGS/vcf/SRW/merged_srw_all_scaffold.vcf.gz \
   -O ~/scratch/temp/filtervariants/srw_variants.vcf.gz \
   --filter-name "mapping_quality" \
   --filter-expression "MQ > 30" \
   --genotype-filter-name "depth" \
   --genotype-filter-expression "DP > 10" \
   --genotype-filter-name "genotype_quality" \
   --genotype-filter-expression "GQ > 30" 
```
### June 27, 2022 ###

**NARW Haplotype Caller**

Yesterday, when trying to concatenate the vcf for NARW, not all indivudals were present in each scaffold file. I double checked and it looks as though the haplotypes were not called for the following:

Scaffold 5 - Both SID missing from gvcf folder
Scaffold 6 - sid181 missing from gvcf
Scaffold 12 - sid17 missing from gvcf
Scaffold 16 - sid17 missing gvcfs
Scaffold 17 - sid17 missing gvcf


**Variant Filtration**

It is unclear to me which direction the < or > operators should be facing - as in, it isn't clear whether you want the good quality sites to match, or fail your filter expressions. I am running both here to check where:

srw_variants.vcf.gz -> MQ > 30 , GQ > 30 , DP >10
srw_variants2.vcf.gz -> MQ < 30 , GQ < 30 , DP < 10

I think it is the latter. I will run SelectVariants to exclude filtered sites. Then I will inspect my VCF files.

```
#!/bin/bash
#SBATCH --job-name=srw_variant_filtration
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=10
#SBATCH --mem=25G
#SBATCH --time=12:00:00

java -Xmx20G -jar ${GATK_JAR} SelectVariants \
   -R ~/projects/def-frasiert/RW_WGS/reference/RWref_HiC_min1Mb.fasta \
   -V ~/scratch/temp/filtervariants/srw_variants.vcf.gz \
   -O ~/scratch/temp/filtervariants/filtered_srw_variants2.vcf.gz \
   --exclude-filtered
```
In the second one, good variants are retained (at least for mapping quality).

Alternatively, I could try filtering the data differently.
```
bcftools filter -i 'COUNT(FORMAT/DP>10)>8 && INFO/MQ>30 && COUNT(FORMAT/GQ>30)>8' -Oz -o ~/scratch/temp/filtervariants/bcf_test_vcf_filter.vcf.gz ~/projects/def-frasiert/RW_WGS/vcf/SRW/merged_srw_all_scaffold.vcf.gz

java -Xmx7G -jar ${GATK_JAR} IndexFeatureFile \
	-F ~/scratch/temp/filtervariants/bcf_test_vcf_filter.vcf.gz 
```

GATK was making it very difficult to filter out format tags (depth per sample, genotype score per sample etc.).
This code for bcftools works very well and produced 11M variants.

**Reheader Bam File**

I also need to reheader my bam files. When I created my bam files, I used the sample name as the PU identifier which is meant to be an instrument/flow cell barcode type name. This prevented BQSR from using the read group appropriately.

I will test all of this on 1 sample with the PU flag removed from the header. If this was effective, we should be able to see all three lanes listed on the BQSR reports. 
```
samtools view -H ~/projects/def-frasiert/RW_WGS/merged_bam/Eau019-merged_marked.bam > header_eau019.sam
sed i "s,^\(@RG.*\)\(\tPU:.*\)\(\tCN:.*\),\1\3," header_eau019.sam
samtools reheader header_eau019.sam ~/projects/def-frasiert/RW_WGS/merged_bam/Eau019-merged_marked.bam > Eau019-merged_marked_newhead.bam
```

**Base Recalibration**

I will proceed with the BaseRecalibrator using the bcf_test_vcf_filter.vcf.gz as known sites and using the Eau019-merged_marked_newhead.bam.

I will use the following as a script to run GATK's BaseRecalibrator followed by the Apply BQSR.
```
#!/bin/bash
#SBATCH --job-name=srw_baserecalibration
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=4
#SBATCH --mem=25G
#SBATCH --time=24:00:00

java -Xmx20G -jar ${GATK_JAR} BaseRecalibrator \
  -R ~/projects/def-frasiert/RW_WGS/reference/RWref_HiC_min1Mb.fasta \
  --known-sites ~/scratch/temp/filtervariants/bcf_test_vcf_filter.vcf.gz \
  -O ~/scratch/temp/filtervariants/recalibration/Eau019_recalibration_report.grp \
  -I ~/scratch/temp/filtervariants/Eau019-merged_marked_newhead.bam

java -Xmx20G -jar ${GATK_JAR} ApplyBQSR \
  -R ~/projects/def-frasiert/RW_WGS/reference/RWref_HiC_min1Mb.fasta \
  -bqsr ~/scratch/temp/filtervariants/recalibration/Eau019_recalibration_report.grp \
  -O ~/scratch/temp/filtervariants/recalibration/Eau019_recal.bam \
  -I ~/scratch/temp/filtervariants/Eau019-merged_marked_newhead.bam
  
java -Xmx20G -jar ${GATK_JAR} BaseRecalibrator \
  -R ~/projects/def-frasiert/RW_WGS/reference/RWref_HiC_min1Mb.fasta \
  --known-sites ~/scratch/temp/filtervariants/bcf_test_vcf_filter.vcf.gz \
  -O ~/scratch/temp/filtervariants/recalibration/Eau019_recalibration_report_after.grp \
  -I ~/scratch/temp/filtervariants/recalibration/Eau019_recal.bam 
```

**NARW fixing gvcfs after haplotypecaller mistakes**

I reran haplotype caller on 6 files that ahd not completed in full before.

I now need to run GenomicsDBImport and GenotypeGVCFs on:

SRW Scaffold 5
SRW Scaffold 6
SRW Scaffold 12
SRW Scaffold 16
SRW Scaffold 17

GenomicsBDImport running. All finished at 8pm. Started GenotypeGVCFS as job *37768837*.

### June 28, 2022 ###

*Summary of where things stand*

All NARW samples have been processed through to gvcf files. (Except HiC_scaffold_15)
All SRW samples have been processed through to gvcf files.
Bowhead has been mapped.

*Currently running*

Running basereaclibration on 1 sample (Eau019) to test how BQSR worked. Job *37731915*
	This job will likely time out, so I will stop it after ApplyBQSR is comeplete and run the final BaseRecalibrator on its own.

*To do*

Finish HaplotypeCaller on SID179132 Scaffold 15
Call Haplotypes on Bowhead
GenomicsDBImport NARW scaffold 15
Reheader all merge marked bam files
GenotypeGVCFS NARW scaffold 15
Merge and filter NARW vcfs & count SNPs
AnalyzeCovariates of SRW 1 sample test


**Finishing Scaffold 15**

I started haplotype caller for SID179132 as job *37848434*

**Bowhead HaplotypeCaller**

I am not sure whether or not I will need a g.vcf file as opposed to a .vcf file. Likely I will need only the vcf for the bowhead, but as the gvcf contains more information and it is easy to convert to a vcf with bcftools gvcf2vcf, I will go ahead and create the gvcf for the bowhead sample.

Started as job *37853552* with code:
```
#!/bin/bash
#SBATCH --job-name=bowhead-haplotype-gvcf
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-23
#SBATCH --cpus-per-task=4
#SBATCH --mem=25G
#SBATCH --time=12:00:00

FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bowhead_sample_scaffolds)

SAMPLE=$(echo $FILE | awk '{print $1}')
SCAFFOLD=$(echo $FILE | awk '{print $2}')

java -Xms20G -Xmx20G -jar ${GATK_JAR} HaplotypeCaller \
  -R ~/projects/def-frasiert/RW_WGS/reference/Eubalaena_glacialis_HiC_min1Mb.fasta \
  -I ~/projects/def-frasiert/RW_WGS/merged_bam/$SAMPLE-merged_marked.bam \
  -O ~/scratch/gvcf/bowhead/$SCAFFOLD/$SAMPLE-$SCAFFOLD.g.vcf.gz \
  -ERC GVCF \
  -L $SCAFFOLD \
  2> ~/projects/def-frasiert/RW_WGS/QC/gvcfQC/$SAMPLE-$SCAFFOLD-hapcall-gvcf.out
```

A few of the files didn't finish I restarted these scaffolds with 24 hours as job *37913533*: Scaffolds 8,9,13,14,16,20

**SRW BaseRecalibrator Test**

As I noted earlier this morning, the srw base recalibration test was taking longer at each stage thatn expected. I stopped the job after it completed the BaseRecalibrator and the ApplyBQSR. I still need to apply the BaseRecalibrator to generate the after statistics. I started this as its own job *37865719*. It could take 10-12 hours. 

In the meantime, I generated the analyze covariate plots of the before set, just to see how it looked.
```
java -Xms7G -Xmx7G -jar ${GATK_JAR}  AnalyzeCovariates \
     -before ~/scratch/temp/filtervariants/recalibration/Eau019_recalibration_report.grp \
     -plots Eau019_AnalyzeCovariates_before.pdf
```

I saved this to my local computer to confirm that yes, the lane identifier for the read groups was fixed. I can now proceed with correcting the header on the other bamfiles.

**GenomicsDBImport NARW Scaffold 15**

Now that SID179132 has finished haplotype calls, I will redo the GenomicsDBImport on Scaffold 15.

This was started as job *37867193* from the gvcf/NARW/HiC_scaffold_15/ using this code:
```
#!/bin/bash
#SBATCH --job-name=narw_genomicsdb_15
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=6
#SBATCH --mem=25G
#SBATCH --time=6:00:00

java -Xms20G -Xmx20G -jar ${GATK_JAR} GenomicsDBImport \
    --genomicsdb-workspace-path ~/scratch/genomicsDB/NARW_workspace/scaff15/ \
    -L HiC_scaffold_15 \
    --sample-name-map ~/scratch/gvcf/NARW/HiC_scaffold_15/HiC_scaffold_15_NARW_samples \
    --tmp-dir ~/scratch/temp \
    --reader-threads 5
```

**New Headers on Bam Files**

I am doing this as a two step process to ensure it all runs properly. The first job ran as *37871528* and created the new header files I will use to replace the existing ones.

```
#!/bin/bash
#SBATCH --job-name=new_headers_pt1
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-23%4
#SBATCH --cpus-per-task=6
#SBATCH --mem=25G
#SBATCH --time=0:30:00

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bam_files)

samtools view -H ~/projects/def-frasiert/RW_WGS/merged_bam/$SAMPLE.bam > header_$SAMPLE.sam

sed i "s,^\(@RG.*\)\(\tPU:.*\)\(\tCN:.*\),\1\3," header_$SAMPLE.sam
```

I checked that this ran as expected and then ran the following as pt2 to replace the existing headers.
I started this as job *37874712*, but it didn't work. Nothing was being written to the output file. 
```
#!/bin/bash
#SBATCH --job-name=new_headers_pt2
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-23%4
#SBATCH --cpus-per-task=6
#SBATCH --mem=25G
#SBATCH --time=1:00:00

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bam_files)

samtools reheader header_$SAMPLE.sam ~/projects/def-frasiert/RW_WGS/merged_bam/$SAMPLE.bam > $SAMPLE_newhead.bam
```
Instead I ran them all in an interactive shell individually.
```
samtools reheader header_EGL013-3qa-merged_marked.sam ~/projects/def-frasiert/RW_WGS/merged_bam/EGL013-3qa-merged_marked.bam > EGL013-3qa-merged_marked_newhead.bam

samtools reheader header_EGL140-1-merged_marked.sam ~/projects/def-frasiert/RW_WGS/merged_bam/EGL140-1-merged_marked.bam > EGL140-1-merged_marked_newhead.bam

samtools reheader header_EGL183-1-merged_marked.sam ~/projects/def-frasiert/RW_WGS/merged_bam/EGL183-1-merged_marked.bam > EGL183-1-merged_marked_newhead.bam

samtools reheader header_EGL254-1-merged_marked.sam ~/projects/def-frasiert/RW_WGS/merged_bam/EGL254-1-merged_marked.bam > EGL254-1-merged_marked_newhead.bam

samtools reheader header_EGL272-1-merged_marked.sam ~/projects/def-frasiert/RW_WGS/merged_bam/EGL272-1-merged_marked.bam > EGL272-1-merged_marked_newhead.bam
```
**GenotypeGVCFs on NARW scaff15**

I started the following code from the NARW_workspace/ as job *37871737*.
```
#!/bin/bash
#SBATCH --job-name=narw_genotype_gvcf
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=15
#SBATCH --cpus-per-task=6
#SBATCH --mem=25G
#SBATCH --time=24:00:00

FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" narw_genDB_inputs)

FOLDER=$(echo $FILE | awk '{print $1}')
SCAFFOLD=$(echo $FILE | awk '{print $2}')

java -Xms20G -Xmx20G -jar ${GATK_JAR} GenotypeGVCFs \
   -R ~/projects/def-frasiert/RW_WGS/reference/Eubalaena_glacialis_HiC_min1Mb.fasta \
   -V gendb://$FOLDER/ \
   -O ~/projects/def-frasiert/RW_WGS/vcf/NARW/narw_HiC_scaffold_$SCAFFOLD.vcf.gz \
   2> ~/projects/def-frasiert/RW_WGS/QC/gvcfQC/narw_HiC_scaffold_$SCAFFOLD-gvcf.out
```
**Analyze Covariates**

I ran this code to generate the covariate plots and saved a copy of the plots locally.

```
java -Xms7G -Xmx7G -jar ${GATK_JAR}  AnalyzeCovariates \
     -before ~/scratch/temp/filtervariants/recalibration/Eau019_recalibration_report.grp \
	 -after ~/scratch/temp/filtervariants/recalibration/Eau019_recalibration_report_after.grp \
     -plots Eau019_AnalyzeCovariates.pdf
```
I think the model did a very good job of recalibration, but because our data was very high quality, the model brought down the base quality very strongly as it assumed every other variant not captured in our filters was an error. 

I think we could do a couple things: 1) open up the variant filters to capture more of the variable sites in the population and maybe the base quality wouldn't be as affected (The filter that requires the site to be present in 8/10 individuals is up for debate - It looks to me like 0/0 genotype is still in the vcf file, so I think this filter is good, if not, it could be removed and I could try again), 2) We can omit this step as we can't for certain say we have enough confidence in capturing the known sites and as such, put such strong weight in assuming all other sites are errors.

**Merging NARW VCFs**

Now that all of the gvcfs have run to completion, I am going to create a merged vcf with information from all samples to be used in filtering variants.

```
bcftools concat narw_HiC_scaffold_1.vcf.gz narw_HiC_scaffold_2.vcf.gz narw_HiC_scaffold_3.vcf.gz narw_HiC_scaffold_4.vcf.gz narw_HiC_scaffold_5.vcf.gz narw_HiC_scaffold_6.vcf.gz narw_HiC_scaffold_7.vcf.gz narw_HiC_scaffold_8.vcf.gz narw_HiC_scaffold_9.vcf.gz narw_HiC_scaffold_10.vcf.gz narw_HiC_scaffold_11.vcf.gz narw_HiC_scaffold_12.vcf.gz narw_HiC_scaffold_13.vcf.gz narw_HiC_scaffold_14.vcf.gz narw_HiC_scaffold_15.vcf.gz narw_HiC_scaffold_16.vcf.gz narw_HiC_scaffold_17.vcf.gz narw_HiC_scaffold_18.vcf.gz narw_HiC_scaffold_19.vcf.gz narw_HiC_scaffold_20.vcf.gz narw_HiC_scaffold_21.vcf.gz narw_HiC_scaffold_22.vcf.gz narw_HiC_scaffold_111.vcf.gz -Oz -o merged_narw_all_scaffold.vcf.gz
```

I will also filter these variants with the same filters I used for the SRW and will compare the number of variants.

```
bcftools filter -i 'COUNT(FORMAT/DP>10)>8 && INFO/MQ>30 && COUNT(FORMAT/GQ>30)>8' -Oz -o ~/scratch/temp/filtervariants/bcf_narw_vcf_filter.vcf.gz ~/projects/def-frasiert/RW_WGS/vcf/NARW/merged_narw_all_scaffold.vcf.gz
```

Filters applied to SRW: 11379204 SNPs
Filters applied to NARW: 4320935 SNPs

### June 29, 2022 ###

Status:
	- All NARW and SRW samples have gvcfs compiled
	
Currently Running:
	- The last of the bowhead scaffolds are being created with haplotype caller. 
	
To do:
	- Reheader .bam files

```
samtools reheader header_EGL276-1-merged_marked.sam ~/projects/def-frasiert/RW_WGS/merged_bam/EGL276-1-merged_marked.bam > EGL276-1-merged_marked_newhead.bam

samtools reheader header_EGL308-1a-merged_marked.sam ~/projects/def-frasiert/RW_WGS/merged_bam/EGL308-1a-merged_marked.bam > EGL308-1a-merged_marked_newhead.bam

samtools reheader header_EGL312-1a-merged_marked.sam ~/projects/def-frasiert/RW_WGS/merged_bam/EGL312-1a-merged_marked.bam > EGL312-1a-merged_marked_newhead.bam

samtools reheader header_EGL336_1b-merged_marked.sam ~/projects/def-frasiert/RW_WGS/merged_bam/EGL336_1b-merged_marked.bam > EGL336_1b-merged_marked_newhead.bam

samtools reheader header_SID179132-merged_marked.sam ~/projects/def-frasiert/RW_WGS/merged_bam/SID179132-merged_marked.bam > SID179132-merged_marked_newhead.bam

samtools reheader header_SID181803-merged_marked.sam ~/projects/def-frasiert/RW_WGS/merged_bam/SID181803-merged_marked.bam > SID181803-merged_marked_newhead.bam

samtools reheader header_SRR1685383-merged_marked.sam ~/projects/def-frasiert/RW_WGS/merged_bam/SRR1685383-merged_marked.bam > SRR1685383-merged_marked_newhead.bam

samtools reheader header_Eau017-merged_marked.sam ~/projects/def-frasiert/RW_WGS/merged_bam/Eau017-merged_marked.bam > Eau017-merged_marked_newhead.bam

samtools reheader header_Eau018-merged_marked.sam ~/projects/def-frasiert/RW_WGS/merged_bam/Eau018-merged_marked.bam > Eau018-merged_marked_newhead.bam

samtools reheader header_Eau019-merged_marked.sam ~/projects/def-frasiert/RW_WGS/merged_bam/Eau019-merged_marked.bam > Eau019-merged_marked_newhead.bam

samtools reheader header_Eau023-merged_marked.sam ~/projects/def-frasiert/RW_WGS/merged_bam/Eau023-merged_marked.bam > Eau023-merged_marked_newhead.bam

samtools reheader header_Eau029-merged_marked.sam ~/projects/def-frasiert/RW_WGS/merged_bam/Eau029-merged_marked.bam > Eau029-merged_marked_newhead.bam

samtools reheader header_Eau034A-merged_marked.sam ~/projects/def-frasiert/RW_WGS/merged_bam/Eau034A-merged_marked.bam > Eau034A-merged_marked_newhead.bam

samtools reheader header_Eau283-merged_marked.sam ~/projects/def-frasiert/RW_WGS/merged_bam/Eau283-merged_marked.bam > Eau283-merged_marked_newhead.bam

samtools reheader header_Eau10b-merged_marked.sam ~/projects/def-frasiert/RW_WGS/merged_bam/Eau10b-merged_marked.bam > Eau10b-merged_marked_newhead.bam

samtools reheader header_Eau7-merged_marked.sam ~/projects/def-frasiert/RW_WGS/merged_bam/Eau7-merged_marked.bam > Eau7-merged_marked_newhead.bam

samtools reheader header_Eau9c-merged_marked.sam ~/projects/def-frasiert/RW_WGS/merged_bam/Eau9c-merged_marked.bam > Eau9c-merged_marked_newhead.bam
```

**Bam Distribution**

I am looking into the data more to understand why the base quality is so unevenly distributed per the BaseRecalibrator reports. I am using the code below to test the one sample we have been working with.

```
java -Xms7G -Xmx7G -jar ${GATK_JAR} CollectBaseDistributionByCycle \
	-CHART ~/scratch/temp/filtervariants/recalibration/Eau019_collect_base_dist_by_cycle.pdf \
	-I ~/projects/def-frasiert/RW_WGS/merged_bam/Eau019-merged_marked_newhead.bam \
    -O ~/scratch/temp/filtervariants/recalibration/Eau019_CollectBaseDistritbutionByCycle.txt
```

**Base Scores**

Turns out base scores from a NovaSeq only produce 4 different quality calls.

### June 30, 2022 ###

Today I did some clean-up on Cedar. I made sure all relevent log files were in the QC/ and deleted superfluous ones. I also reorganized to erased empty directories and tidy extra files.

I downloaded one of the merged vcf files to inspect on IGV and in R, but I need to create an index for IGV.

### July 1, 2022 ###

Today I will:
	- Create Index for the merged NARW
	- Filter the merged NARW and merged SRW files and create indexes of those VCF
	- Download these files
	- Create plots in R of before and after filtering
	

**Filter and Index VCF files**

```
module load bcftools/1.11

bcftools filter -i 'COUNT(FORMAT/DP>10)>8 && INFO/MQ>30 && COUNT(FORMAT/GQ>30)>8' -Oz -o ~/scratch/temp/filtervariants/merged_narw_vcf_filtered.vcf.gz ~/projects/def-frasiert/RW_WGS/vcf/NARW/merged_narw_all_scaffold.vcf.gz

bcftools filter -i 'COUNT(FORMAT/DP>10)>8 && INFO/MQ>30 && COUNT(FORMAT/GQ>30)>8' -Oz -o ~/scratch/temp/filtervariants/merged_srw_vcf_filtered.vcf.gz ~/projects/def-frasiert/RW_WGS/vcf/SRW/merged_srw_all_scaffold.vcf.gz
```

```
java -Xmx7G -jar ${GATK_JAR} IndexFeatureFile \
	-F ~/projects/def-frasiert/RW_WGS/vcf/NARW/merged_narw_all_scaffold.vcf.gz
	
java -Xmx7G -jar ${GATK_JAR} IndexFeatureFile \
	-F ~/scratch/temp/filtervariants/merged_narw_vcf_filtered.vcf.gz 
	
java -Xmx7G -jar ${GATK_JAR} IndexFeatureFile \
	-F ~/scratch/temp/filtervariants/merged_srw_vcf_filtered.vcf.gz
```
```
module load plink/1.9b_6.21-x86_64

plink --vcf ~/scratch/temp/filtervariants/merged_srw_vcf_filtered.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 \
--out srw_plink

plink --vcf ~/scratch/temp/filtervariants/merged_srw_vcf_filtered.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--extract srw_plink.prune.in \
--make-bed --pca --out srw_plink

plink --vcf ~/scratch/temp/filtervariants/merged_narw_vcf_filtered.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 \
--out narw_plink

plink --vcf ~/scratch/temp/filtervariants/merged_narw_vcf_filtered.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--extract narw_plink.prune.in \
--make-bed --pca --out narw_plink
```

### July 3, 2022 ###

**Bedtools Intersect**

I ran the following code to pull out variants that overlapped with the masked repeats .gff file form DNAZoo.

I saved this to a simple file called out.vcf, ran gzip on the file and downloaded it locally, I then reran my plots on it for inspection.

```
bedtools intersect -v -a bcf_narw_vcf_filter.vcf.gz -b ~/scratch/Eubalaena_glacialis_HiC.repeatmasker.trf.windowmasker.gff.gz -wa -header > out.vcf

gzip out.vcf
```

1334916 variants remain

There still appears to be a lot of alleles with very low frequency. I am going to try to filter based on AB at each locus.


### July 4, 2022 ###

**New Filtering**

The plots I produced suggest that MAF is a problem. I am going to retry variant filtration in a twp step process using GATK. First I will filter and select variants based on the format fields that will filter out individual calls where depth is too low and allelic balance is off (a heterozygote call where the allele balance falls outside 0.2-0.8. 
```
#!/bin/bash
#SBATCH --job-name=narw_variant_filtration
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=10
#SBATCH --mem=25G
#SBATCH --time=12:00:00

java -Xmx20G -jar ${GATK_JAR} VariantFiltration \
   -R ~/projects/def-frasiert/RW_WGS/reference/Eubalaena_glacialis_HiC_min1Mb.fasta \
   -V ~/projects/def-frasiert/RW_WGS/vcf/NARW/merged_narw_all_scaffold.vcf.gz \
   -O ~/scratch/temp/filtervariants/narw_select_variants_ind.vcf.gz \
   --genotype-filter-name "depth" \
   --genotype-filter-expression "DP < 10" \
   --genotype-filter-name "genotype_quality" \
   --genotype-filter-expression "GQ < 30" \
   --genotype-filter-name "allelic_balance" \
   --genotype-filter-expression "isHet == 1 && 0.8 < AD[1]/DP < 0.2" \
   

java -Xmx20G -jar ${GATK_JAR} SelectVariants \
   -R ~/projects/def-frasiert/RW_WGS/reference/Eubalaena_glacialis_HiC_min1Mb.fasta \
   -V ~/scratch/temp/filtervariants/narw_select_variants_ind.vcf.gz \
   -O ~/scratch/temp/filtervariants/narw_filtered_variants_ind.vcf.gz \
   --set-filtered-gt-to-nocall
   
```

The final genotype filter wasn't working. But after discussion with Michael et al. We don't need to worry about allele frequency. Poor calls will be filtered out by Genotype Quality.

I will now filter as follows:

1) Use bedtools intersect to mask repetitive regions
	- Save masked regions too (if possible)
2) Plot metrics
3) Apply filters
		- vcftools to mask DP and GQ filters
		- bcftools to mask MQ, missing data and max DP filters 
4) Plot metrics

Each of these steps do not take very long, but bedtools requires more memory. So I will run bedtools as a scheduled job and then the filtering on an interactive node:

```
module load mugqic/bedtools/2.30.0 bcftools/1.11 mugqic/vcftools/0.1.14


bedtools intersect -v -a ~/projects/def-frasiert/RW_WGS/vcf/NARW/merged_narw_all_scaffold.vcf.gz -b ~/scratch/Eubalaena_glacialis_HiC.repeatmasker.trf.windowmasker.gff.gz -wa -header > ~/scratch/temp/filtervariants/narw_masked.vcf

gzip narw_masked.vcf

vcftools --vcf ~/scratch/temp/filtervariants/narw_masked.vcf.gz --out ~/scratch/temp/filtervariants/narw_masked_formatfilter.vcf.gz --minDP 10 --minGQ 30

bcftools filter -i 'COUNT(FORMAT/GT="mis")<=3 && INFO/MQ>30 && INFO/DP<790' -Oz -o ~/scratch/temp/filtervariants/narw_masked_filtered.vcf.gz   ~/scratch/temp/filtervariants/narw_masked_formatfilter.vcf.gz
```

SNP Counts:
Raw VCF - 5620253
Masked VCF - 1526834
VCFtools  - 1526834 (only marked genotypes)
BCFTools - 1407025

** Started Mapping Bowhead to SRW **

Started mapping bowhead southern right whale genome as job *38548442*.

** Started Filtering on SRW **

I started the masking with bedtools using the repeat masking file from DNAZoo.
Started as job *38548782*.

```
vcftools --gzvcf ~/scratch/temp/filtervariants/srw_masked.vcf.gz --out ~/scratch/temp/filtervariants/srw_masked_formatfilter.vcf.gz --minDP 10 --minGQ 30

bcftools filter -i 'COUNT(FORMAT/GT="mis")<3 && INFO/MQ>30 && INFO/DP<558' -Oz -o ~/scratch/temp/filtervariants/srw_masked_filtered.vcf.gz   ~/scratch/temp/filtervariants/srw_masked.vcf.gz
```

SNP Counts:
Raw VCF - 21094633
Masked VCF - 8320997
VCFtools  - 8320997 (only marked genotypes)
BCFTools - 8139503

### July 5, 2022 ###

**Bowhead Progress**

Yesterday, I finishing mapping bowhead reads to the Southern Right Whale genome. Today I am merging the three lanes and marking duplicates. Started as job *38600113*.

Started haplotype calling as job *38619570*.

**Bedtools Intersect**

I'm not entirely sure that bedtools is doing what I want it to. I tried removing the -wa filter and got the same results, so I am confident in the repeat masking with bedtools.

I also want to extract SNPs from the masked regions to look at the summary statistics. I will do this using bedtools subtract
```
bedtools subtract -a ~/projects/def-frasiert/RW_WGS/vcf/NARW/merged_narw_all_scaffold.vcf.gz \
        -b ~/scratch/Eubalaena_glacialis_HiC.repeatmasker.trf.windowmasker.gff.gz \
        -header > ~/scratch/temp/filtervariants/narw_repeats.vcf
```

RepeatMasker creates a file of repeatitive regions in the genome. 
To mask the repeats, we need use bedtools intersect to identify the regions where these do not overlap. We also need to use the unique option to extract the VCF sites only once as the repeat file can have overlapping regions. 

-v creates an output of the vcf that omits any overlap with the .gff.

I saved a copy of the good filtered .vcf.gz in the ~/projects/.../vcf/NARW/ and SRW/
The NARW repeat vcf is saved as ~/scratch/temp/filtervariants/test_wau_narw_masked.vcf.gz

**PCA on filtered VCFs**

```
module load plink/1.9b_6.21-x86_64

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/NARW/narw_masked_filtered.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 \
--out narw_plink

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/NARW/narw_masked_filtered.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--extract narw_plink.prune.in \
--make-bed --pca --out narw_plink

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/SRW/srw_masked_filtered.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 \
--out srw_plink

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/SRW/srw_masked_filtered.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--extract srw_plink.prune.in \
--make-bed --pca --out srw_plink
```

Pruned for LD

In NARW:
	- variants before: 1407025
	- variants after pruning: 48382
In SRW:
	- variants before: 8139503
	- variants after pruning: 223255


### July 6, 2022 ###

--depth Mean depth per individual
--TsTv-summary transition and transversions
--site-pi nucleotide diversity per site
--het heterozygosity per individual
--hardy pvalue for HWE 
--TajimaD <integer> Tajima's D in bins
--missing-indv
--missing-site

```
vcftools --gzvcf ~/scratch/temp/filtervariants/narw_masked_filtered.vcf.gz --out narw_vcftools_stats --depth

vcftools --gzvcf ~/scratch/temp/filtervariants/narw_masked_filtered.vcf.gz --out narw_vcftools_stats --TsTv-summary 

vcftools --gzvcf ~/scratch/temp/filtervariants/narw_masked_filtered.vcf.gz --out narw_vcftools_stats --site-pi

vcftools --gzvcf ~/scratch/temp/filtervariants/narw_masked_filtered.vcf.gz --out narw_vcftools_stats --het

vcftools --gzvcf ~/scratch/temp/filtervariants/narw_masked_filtered.vcf.gz --out narw_vcftools_stats --hardy

vcftools --gzvcf ~/scratch/temp/filtervariants/narw_masked_filtered.vcf.gz --out narw_vcftools_stats --missing-indv

vcftools --gzvcf ~/scratch/temp/filtervariants/narw_masked_filtered.vcf.gz --out narw_vcftools_stats --missing-site

vcftools --gzvcf ~/scratch/temp/filtervariants/narw_masked_filtered.vcf.gz --out narw_vcftools_stats --TajimaD 100000
```


maff filtering


PLot summary stats of sites that are masked


repeats 110-131 HiC_scaffold_1