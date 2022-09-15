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

### July 9, 2022 ###

I wrote code in python for using scikit allel to calculate ROH, but it was taking a long time to run in th jupyter notebook.

```
# load module

module unload mugqic/python/3.9.1
module load scipy-stack

# virtualenv --no-download PY_ENV

source PY_ENV/bin/activate

# source deactivate


# Inside job reactivate environment
# I will install modules now and they can be imported within scripts

narw_vcf = allel.read_vcf(~/projects/def-frasiert/RW_WGS/vcf/NARW/narw_masked_filtered.vcf.gz, fields=['variants/*','calldata/*'])

pi = []
snp_per_window = []
chrom = []
for SCAFFOLD in np.unique(narw_vcf['variants/CHROM']):
    pos = narw_vcf['variants/POS'][narw_vcf['variants/CHROM'] == SCAFFOLD]
    gt = allel.GenotypeArray(narw_vcf['calldata/GT'][narw_vcf['variants/CHROM'] == SCAFFOLD])
    diversity = allel.windowed_diversity(pos, gt.count_alleles(), size = 10000, start=1, stop=None, is_accessible=None)
    snp_per_window.extend(diversity[3])
    pi.extend(diversity[0])
    chrom.extend(np.repeat(SCAFFOLD, len(diversity[0])))

# Save diversity per window to a .csv file
narw_diversity_wide = pd.DataFrame((chrom, pi, snp_per_window))
narw_sequence_diversity = narw_diversity_wide.transpose()
narw_sequence_diversity.columns = ['chromosome', 'pi', 'SNP_per_window']
pd.DataFrame(narw_sequence_diversity).to_csv('narw_sequence_diversity_10kb.csv') 
```

I want to convert it to a script and run it on Cedar.

First I need to activate a new python environment where I will be able to load modules. 

A lot of troubleshooting was had followed by a cluster slow down and delays.

### July 11, 2022 ###

**Run ROH**

To successfully run the python scripts, I did the following:
```
# Run from inside temp/filtervariants/

module load python
module unload mugqic/python/3.9.1
module load scipy-stack

source ../PY_ENV/bin/activate

sbatch ../../run_python.sh
```

Here is the code for these scripts:

```
### run_python.sh ###

#!/bin/bash
#SBATCH --job-name=roh
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=24:00:00

python roh_python.py

# ------------------------------------------------ #
### roh_python.py ###

# Import Libraries
import numpy as np
import scipy
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import allel; print('scikit-allel', allel.__version__)
import hmmlearn

narw_vcf = allel.read_vcf('narw_masked_filtered.vcf.gz', fields=['variants/*','calldata/*'])

sample_no = []
scaff = []
narw_percent_roh = []
narw_num_roh = []
for SCAFFOLD in np.unique(narw_vcf['variants/CHROM']):
    narw_gt = allel.GenotypeArray(narw_vcf['calldata/GT'][narw_vcf['variants/CHROM'] == SCAFFOLD])
    for IND in range(0,12):
		gv = narw_gt[:,IND]
        pos = narw_vcf['variants/POS'][narw_vcf['variants/CHROM'] == SCAFFOLD]
        scaffold_size = pos[len(pos)-1]
        roh = allel.roh_mhmm(gv, pos, phet_roh=0.001, phet_nonroh=(0.0025, 0.01), transition=1e-06, min_roh=0, contig_size=scaffold_size)
        sample_no.append(IND)
        scaff.append(SCAFFOLD)
        narw_percent_roh.append(roh[1])        narw_num_roh.append(len(roh[0]))
		
# Save ROH data per sample per scaffold to a .csv file
narw_roh_wide = pd.DataFrame((sample_no, scaff, narw_num_roh, narw_percent_roh))
narw_roh_persamp = narw_roh_wide.transpose()
narw_roh_persamp.columns = ['Sample_No', 'Scaffold', 'Number_ROHs', 'Percent_Scaffold_in_ROH']
pd.DataFrame(narw_roh_persamp).to_csv('narw_roh_persamp.csv')
```

The data seemed wrong to me, and while I used the defaults, I realized thee main arguement that was a problem was the min_roh. I reran the command abve on both species but used the following:

```
# Import Libraries
import numpy as np
import scipy
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import allel; print('scikit-allel', allel.__version__)
import hmmlearn

narw_vcf = allel.read_vcf('narw_masked_filtered.vcf.gz', fields=['variants/*','calldata/*'])

sample_no = []
scaff = []
narw_percent_roh = []
narw_num_roh = []
for SCAFFOLD in np.unique(narw_vcf['variants/CHROM']):
    narw_gt = allel.GenotypeArray(narw_vcf['calldata/GT'][narw_vcf['variants/CHROM'] == SCAFFOLD])
    for IND in range(0,12):
        gv = narw_gt[:,IND]
        pos = narw_vcf['variants/POS'][narw_vcf['variants/CHROM'] == SCAFFOLD]
        scaffold_size = pos[len(pos)-1]
        roh = allel.roh_mhmm(gv, pos,min_roh=100000, contig_size=scaffold_size)
        sample_no.append(IND)
        scaff.append(SCAFFOLD)
        narw_percent_roh.append(roh[1])
        narw_num_roh.append(len(roh[0]))

# Save ROH data per sample per scaffold to a .csv file
narw_roh_wide = pd.DataFrame((sample_no, scaff, narw_num_roh, narw_percent_roh))
narw_roh_persamp = narw_roh_wide.transpose()
narw_roh_persamp.columns = ['Sample_No', 'Scaffold', 'Number_ROHs', 'Percent_Scaffold_in_ROH']
pd.DataFrame(narw_roh_persamp).to_csv('narw_roh_persamp.csv')
```

** Determining Sex Scaffolds **

From qualimap, it was quite difficult to identify depth of single individuals at each chromosome. There should have been a report generated for each individual, but I didn't see one.

I am running qualimap on 2 indivivuals to assess differences in read depth across scaffolds in hope of identifying the sex scaffolds. In the male, EGL272-1, read depth on the sex scaffolds (X and Y) should be roughly half that of the other scaffolds.

In the female, there should be one scaffold with very low depth (ideally zero) - this will be the Y chromosome and should be the same as one of the two identified in the male.


### July 14, 2022 ###

** Mapping SRW to NARW **

I have trimmed SRW reads, but I want to map them to the NARW genome too.
Started with mapping.
```
#!/bin/bash
#SBATCH --job-name=srw_map_to_narw
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-30%12
#SBATCH --cpus-per-task=28
#SBATCH --mem=15G
#SBATCH --time=05:00:00

SAMPLELANE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/projects/def-frasiert/RW_WGS/srw_fastq_list)
SAMPLE=$(grep -oP '(?<=\b).*?(?=_2)' <<< "$SAMPLELANE")
LANE=${SAMPLELANE: -4}
READGROUP="@RG\tID:${LANE}\tSM:${SAMPLE}\tLB:RW\tCN:MCGILL\tPL:ILLUMINA"

bwa mem -M -t 28 \
  -R $READGROUP \
  ~/projects/def-frasiert/RW_WGS/reference/Eubalaena_glacialis_HiC_min1Mb.fasta \
  ~/scratch/trimmed_reads/$SAMPLELANE-FP.fq.gz \
  ~/scratch/trimmed_reads/$SAMPLELANE-RP.fq.gz > \
  ~/scratch/temp_mapping/$SAMPLELANE-aln_min1Mb_SRWonNARW.sam \
  2> ~/projects/def-frasiert/RW_WGS/QC/std_out/$SAMPLELANE-SRWonNARW-bwa.err
```



### July 13, 2022 ###

Yesterday I finished mapping SRW to NARW and finishing sorting those bam file. 

Today I am going to delete the .sam files first, then I am going to run a script to merge the bam files and mark duplicates.

### July 14, 2022 ###

Running HaplotypeCaller on SRW mapped to NARW.

**KING**

I installed KING in our project folder. It is a precompiled binary, so I saved it and unzipped with tar in a new directory called programs/.

I needed to run a script supplied by the digital alliance to change some library configurations to make it work. This is all the code I used to set it up.

```
## Inside programs/

wget https://www.kingrelatedness.com/Linux-king.tar.gz
tar -xzvf Linux-king.tar.gz

setrpaths.sh --path ~/projects/def-frasiert/RW_WGS/programs/
```

### July 15, 2022 ###

**Checking HaplotypeCaller on SRW mapped to NARW**

Even with 8 hours, a few files were timing out. It has only been about 24 hours so this file is not expected to be finished yet, but as jobs are ending, I am chacking the tail of the .out file to confirm they finished in full. 

The following is a list I will continue to add to of files that will need to be run with more time and/or resources:
	- Eau017 HiC_scaffold_20
	- Eau018 HiC_scaffold_13
	- Eau018 HiC_scaffold_20
	- Eau019 HiC_scaffold_6
	- Eau019 HiC_scaffold_13
	- Eau019 HiC_scaffold_14
	- Eau019 HiC_scaffold_16
	- Eau019 HiC_scaffold_20	


**Running KING**

I first tested running KING on a single scaffold data to confirm it was running as desired.

Now I am going to run PLINK without linkage pruning (as recommended by KING) on the entire NARW filtered vcf file. I will then need to edit the .bim file to remove the 'HiC_scaffold_' prefix from the chromosome names and then run KING --kinship and --ibdseq specifying the --sexchr 21.

```
## Inside an interactive node

module load plink/1.9b_6.21-x86_64

## From inside ~/scratch/relatedness/
## Run PLINK to generate .bed, .fam and .bim files

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/NARW/narw_masked_filtered.vcf.gz --allow-extra-chr --out narw

## Edit the .bim file to remove the HiC_scaffold_ prefix from column 1

sed 's,HiC_scaffold_\(.*\),\1,' narw.bim > narw_edit.bim

## Run KING

~/projects/def-frasiert/RW_WGS/programs/king -b ~/scratch/relatedness/narw.bed --fam ~/scratch/relatedness/narw.fam --bim ~/scratch/relatedness/narw_edit.bim --kinship --degree 3 --prefix narw_kinship --sexchr 21
```

The results are saved in relatedness/narw_kinship.kin0

Two pairs of samples have a kinship coefficient of >0.25 which is the equivalent of parent-offspring or full siblings. 

Pairs:  EGL140-1 and EGL312-1a
		EGL272-1 and EGL276-1
		
One additional pair also had a kinship coefficient value which falls in KING's bin for 1st degree relationship (0.177-0.354).

Pair:   EGL00252-1 and EGL276-1 : 0.1773 
		
Since EGL140-1 and EGL272-1 are both males, I think we could exlude these two and keep the remaining 10 female samples. There was also no real difference in mean depth between the related pairs to prompt keeping one over the other either.

Although, if we remove EGL276-1 instead, we also remove the potential relationship with EGL00252-1.

Going into the admixture analysis, I think I remove the sex scaffold and the 2 sample: EGL140-1 and EGL276-1 from the dataset.


I am going to run the same analyses on Southern Right Whales.

```
## Inside an interactive node

module load plink/1.9b_6.21-x86_64

## From inside ~/scratch/relatedness/
## Run PLINK to generate .bed, .fam and .bim files

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/SRW/srw_masked_filtered.vcf.gz --allow-extra-chr --out srw

## Edit the .bim file to remove the HiC_scaffold_ prefix from column 1 and change a few chromosome names so as to clearly specify sex scaffolds.

sed 's,HiC_scaffold_\(.*\),\1,' srw.bim > srw_edit.bim
sed -i 's,21\(\t.*\),1\1,' srw_edit.bim
sed -i 's,8\(\t.*\),21\1,' srw_edit.bim

## Run KING

~/projects/def-frasiert/RW_WGS/programs/king -b ~/scratch/relatedness/srw.bed --fam ~/scratch/relatedness/srw.fam --bim ~/scratch/relatedness/srw_edit.bim --kinship --degree 5 --prefix srw_kinship --sexchr 8
```

There was a strong relationship between Eau10b and Eau019 (kinship coeefficient >0.31). We will removed the individual with poorer quality data for admixture analyses (Eau019).

**Removing Individuals and Sex Scaffolds and HiC_scaffold_111**

```
module load mugqic/vcftools/0.1.14

vcftools --remove-indv EGL140-1 --remove-indv EGL276-1 --not-chr HiC_scaffold_21 --not-chr HiC_scaffold_111 --gzvcf ~/projects/def-frasiert/RW_WGS/vcf/NARW/narw_masked_filtered.vcf.gz --recode --out narw_unrelated_filtered.vcf.gz

vcftools --remove-indv Eau019 --not-chr HiC_scaffold_8 --gzvcf ~/projects/def-frasiert/RW_WGS/vcf/SRW/srw_masked_filtered.vcf.gz --recode --out srw_unrelated_filtered.vcf.gz

```

**LD pruning**

```
module load plink/1.9b_6.21-x86_64

#SRW
plink --vcf srw_unrelated_filtered.vcf.gz.recode.vcf --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 \
--out srw_plink

plink --vcf srw_unrelated_filtered.vcf.gz.recode.vcf --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--extract srw_plink.prune.in \
--make-bed --pca --out srw_plink

# NARW
plink --vcf narw_unrelated_filtered.vcf.gz.recode.vcf --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 \
--out narw_plink

plink --vcf narw_unrelated_filtered.vcf.gz.recode.vcf --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--extract narw_plink.prune.in \
--make-bed --pca --out narw_plink
```

The number of SNPs retained after LD pruning was:
	- NARW: 40426
	- SRW: 221557
	
This is quite low for NARW, typically would want to see around 100,000 for ADMIXTURE type analyses. This could be achieved by increasing the r2 threshold from 0.1.

**LD Decay**

The program PopLDdecay can be used to understand linkage decay across the genome. This program can be loaded as a module.

```
module load poplddecay/3.41

PopLDdecay    -InVCF  srw_unrelated_filtered.vcf.gz.recode.vcf  -OutStat LDdecay

```

The program is supposed to be able to plot the data, but I'm not sure this will work on the cluster. I'm sure I can plot myself in R if the code doesn't work.


## July 16, 2022 ##

Downloaded and tried to run the Plot_OnePop.pl from PopLDdecay. It threw an error. Will try to figure it out Monday. Brief overview:ran wget to save file to programs/ and then tried to run. Error was:
```
Bareword found where operator expected at /home/crossman/projects/def-frasiert/RW_WGS/programs/Plot_OnePop.pl line 9, near ""en" data"
        (Missing operator before data?)
Bareword found where operator expected at /home/crossman/projects/def-frasiert/RW_WGS/programs/Plot_OnePop.pl line 9, near ""auto" data"
        (Missing operator before data?)
Bareword found where operator expected at /home/crossman/projects/def-frasiert/RW_WGS/programs/Plot_OnePop.pl line 9, near ""light" data"
        (Missing operator before data?)
Bareword found where operator expected at /home/crossman/projects/def-frasiert/RW_WGS/programs/Plot_OnePop.pl line 9, near ""dark" data"
        (Missing operator before data?)
Can't modify numeric lt (<) in scalar assignment at /home/crossman/projects/def-frasiert/RW_WGS/programs/Plot_OnePop.pl line 9, near ""en" data"
syntax error at /home/crossman/projects/def-frasiert/RW_WGS/programs/Plot_OnePop.pl line 9, near ""en" data"
Excessively long <> operator at /home/crossman/projects/def-frasiert/RW_WGS/programs/Plot_OnePop.pl line 21.
```


**PLINK with larger r2 threshold for NARW**

```
module load plink/1.9b_6.21-x86_64


# NARW
plink --vcf narw_unrelated_filtered.vcf.gz.recode.vcf --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.2 \
--out narw_plink_r.2

plink --vcf narw_unrelated_filtered.vcf.gz.recode.vcf --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--extract narw_plink_r.2.prune.in \
--make-bed --pca --out narw_plink_r.2
```

**ADMIXTURE**

I want to run 10 iterations of admixture for each K of 1-5.

Had to change scaff names to numbers only
```
module load admixture/1.3.0

#!/bin/bash
#SBATCH --job-name=narw_admixture
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-10
#SBATCH --cpus-per-task=22
#SBATCH --mem=10G
#SBATCH --time=24:00:00

r=${SLURM_ARRAY_TASK_ID} 

prefix=narw_plink_r.2

for K in {1..5};
do
	admixture --cv -s ${RANDOM} ${prefix}.bed $K -j20
	mv ${prefix}.${K}.Q out/${prefix}.K${K}r${r}.Q
done
```


### July 17, 2022 ###

**GenomicsDBImport**

I confirmed that all the haplotypecaller files had run to completion.

I then started running genomicsDBimport to create the genomics databases.

### July 18, 2022 ###

**GenomicsDBImport**

After speaking with Michael, we decided it would be best to generate the genomics database files for both NARW and SRW mapped to NARW together as a single dataset. That way, alleles fixed in each population would be included in the final merged vcf. 

To do this, I copied the NARW g.vcf.gz files from their NARW/ to SRWonNARW/ and ran the genomicsDBimport scripts from with of the scaffold directories. I gave each script twice as much time to run and will check the status this evening.

I will add the bowhead to this database this afternoon.

**Visualizing Admixture**

I ran pong locally in the windows command prompt using the following:
```
pong -m narw_k1-5_Qfilemap --greedy -s .95
```
It generated a series of plots, but nothing really interesting - as anticipated.

The cross validation of error values will help determine the most appropriate K value.

```
#Generate file with cv error rates from log file

grep "^CV error" ../narw_admixture-39683972.out | sed 's,CV error (K=\(.*\)): \(.*\),\1\t\2,' > narw_plink_r.2.CV.txt
```
The data best fit the models assuming a single population as indicated by the lack of convergence of the admixture plots as well as the lowest CV error at K=1.


**Southern Right Whale Admixture**

Running Southern Right Whale admixture alone. 

ADMIXTURE SCRIPT:
```
#!/bin/bash
#SBATCH --job-name=srw_admixture
#SBATCH --output=srw_admixture_log.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=22
#SBATCH --mem=10G
#SBATCH --time=1:00:00

prefix=srw_plink

for r in {1..10}; do for K in {1..5};
do
        admixture --cv -s ${RANDOM} ${prefix}.bed $K -j20
        mv ${prefix}.${K}.Q out/${prefix}.K${K}r${r}.Q
done; done
```

Now I need to plot the data. First I will plot the admixture plots with pong. I will create a file map and then transfer the files locally to use pong.

```
# Create file map
createQmap(){
local r=$1
local K=$2
awk -v K=$K -v r=$r -v file=${prefix}.K${K}r${r} 'BEGIN{ \
printf("K%dr%d\t%d\t%s.Q\n",K,r,K,file)
}' >> ${prefix}_k1-5_Qfilemap
}
export -f createQmap
for K in {1..5}; do for r in {1..10}; do createQmap $r $K; \
done; done


# Run pong 
pong -m srw_k1-5_Qfilemap --greedy -s .95


#Generate file with cv error rates from log file
grep "^CV error" ../srw_admixture_log.out | sed 's,CV error (K=\(.*\)): \(.*\),\1\t\2,' > srw_plink_r.2.CV.txt
```

Just as with the NARW, there was clearly 1 population as indicated by the lack of convergence of the admixture plots as well as the lowest CV error at K=1.

**Add to GenomicsDatabases**

If I use a slightly newer version of GATK (v4.1.8 instead of v4.1.0), I can add samples to a genomics database and then genotype GVCFs across all of them.

First I will unload the older verison of GATK and load the newer one.

```
module unload mugqic/GenomeAnalysisTK/4.1.0
module load mugqic/GenomeAnalysisTK/4.1.8
```

Tested on scaff111:

```
java -Xms7G -Xmx7G -jar ${GATK_JAR} GenomicsDBImport \
--genomicsdb-update-workspace-path ~/scratch/genomicsDB/NARW_workspace/scaff111/ \
-V ~/scratch/gvcf/bowhead/HiC_scaffold_111/SRR1685383-HiC_scaffold_111.g.vcf.gz \
-L HiC_scaffold_111 \
--tmp-dir ~/scratch/temp
```

Writing a script for running through all scaffolds:

```
#!/bin/bash
#SBATCH --job-name=add_bowhead_array
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-22
#SBATCH --cpus-per-task=6
#SBATCH --mem=25G
#SBATCH --time=12:00:00

SCAFF=${SLURM_ARRAY_TASK_ID}

FILE=~/scratch/gvcf/bowhead/HiC_scaffold_${SCAFF}/SRR1685383-HiC_scaffold_${SCAFF}.g.vcf.gz

java -Xms20G -Xmx20G -jar ${GATK_JAR} GenomicsDBImport \
        --genomicsdb-update-workspace-path ~/scratch/genomicsDB/NARW_workspace/scaff${SCAFF}/ \
        -V ${FILE} \
        --tmp-dir ~/scratch/temp
```

I also ran this updating the SRW_on_NARW databases too.

**Genotype GVCFs**

Ran this script to call genotypes on all files mapped to NARW.
```
#!/bin/bash
#SBATCH --job-name=genotype_gvcf
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-23
#SBATCH --cpus-per-task=6
#SBATCH --mem=25G
#SBATCH --time=24:00:00

FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" srw_on_narw_gendb_inputs)

FOLDER=$(echo $FILE | awk '{print $1}')
SCAFFOLD=$(echo $FILE | awk '{print $2}')

java -Xms20G -Xmx20G -jar ${GATK_JAR} GenotypeGVCFs \
   -R ~/projects/def-frasiert/RW_WGS/reference/Eubalaena_glacialis_HiC_min1Mb.fasta \
   -V gendb://$FOLDER/ \
   -O ~/scratch/vcf/all_on_narw_HiC_scaffold_$SCAFFOLD.vcf.gz \
   2> ~/projects/def-frasiert/RW_WGS/QC/gvcfQC/all_on_narw_HiC_scaffold_$SCAFFOLD-gvcf.out
```

### July 19, 2022 ###

**Genotyping gvcfs**

I checked on the status of the genotype gvcfs this morning. Scaff 16 and scaff 20 were unlikely to finish. after 12 hours of running, I killed those jobs and restarted with 48 hours instead of 24. This could mean they won't be done until Thursday morning! 

**Preparing Code**

Once I have called genotypes, I need to merge my vcf files using *bcftools concat*
```
module load bcftools/1.11

bcftools concat all_on_narw_HiC_scaffold_1.vcf.gz all_on_narw_HiC_scaffold_2.vcf.gz all_on_narw_HiC_scaffold_3.vcf.gz all_on_narw_HiC_scaffold_4.vcf.gz all_on_narw_HiC_scaffold_5.vcf.gz all_on_narw_HiC_scaffold_6.vcf.gz all_on_narw_HiC_scaffold_7.vcf.gz all_on_narw_HiC_scaffold_8.vcf.gz all_on_narw_HiC_scaffold_9.vcf.gz all_on_narw_HiC_scaffold_10.vcf.gz all_on_narw_HiC_scaffold_11.vcf.gz all_on_narw_HiC_scaffold_12.vcf.gz all_on_narw_HiC_scaffold_13.vcf.gz all_on_narw_HiC_scaffold_14.vcf.gz all_on_narw_HiC_scaffold_15.vcf.gz all_on_narw_HiC_scaffold_16.vcf.gz all_on_narw_HiC_scaffold_17.vcf.gz all_on_narw_HiC_scaffold_18.vcf.gz all_on_narw_HiC_scaffold_19.vcf.gz all_on_narw_HiC_scaffold_20.vcf.gz all_on_narw_HiC_scaffold_21.vcf.gz all_on_narw_HiC_scaffold_22.vcf.gz all_on_narw_HiC_scaffold_111.vcf.gz -Oz -o ~/projects/def-frasiert/RW_WGS/vcf/ALL/merged_all_on_narw.vcf.gz
```

With merged VCF files, I need to perform some filtering. First I will mask repeat regions defined in the .gff file. Then, I will filters the variants sites.

```
module load mugqic/bedtools/2.30.0 bcftools/1.11 mugqic/vcftools/0.1.14


bedtools intersect -v -a ~/projects/def-frasiert/RW_WGS/vcf/ALL/merged_all_on_narw.vcf.gz -b ~/scratch/Eubalaena_glacialis_HiC.repeatmasker.trf.windowmasker.gff.gz -wa -header > ~/scratch/temp/filtervariants/all_on_narw_masked.vcf

gzip ~/scratch/temp/filtervariants/all_on_narw_masked.vcf

vcftools --gzvcf ~/scratch/temp/filtervariants/all_on_narw_masked.vcf.gz --out ~/scratch/temp/filtervariants/all_on_narw_masked_formatfilter.vcf.gz --minDP 10 --minGQ 30
```

Finish the filtering with appropriate mapping quality and max depth.
2x median value seems to work again. 2x depth = 1370
```
bcftools filter -i 'COUNT(FORMAT/GT="mis")<6 && INFO/MQ>30 && INFO/DP<1370' -Oz -o ~/scratch/temp/filtervariants/all_on_narw_masked_filtered.vcf.gz  ~/scratch/temp/filtervariants/all_on_narw_masked.vcf.gz
```

Removing Related Individuals:

Now I am going to run PLINK without linkage pruning (as recommended by KING) on the entire NARW filtered vcf file. I will then need to edit the .bim file to remove the 'HiC_scaffold_' prefix from the chromosome names and then run KING --kinship and --ibdseq specifying the --sexchr 21.

```
## Inside an interactive node

module load plink/1.9b_6.21-x86_64

## From inside ~/scratch/relatedness/
## Run PLINK to generate .bed, .fam and .bim files

plink --vcf ~/scratch/temp/filtervariants/all_on_narw_masked_filtered.vcf.gz  --allow-extra-chr --out all_on_narw

## Edit the .bim file to remove the HiC_scaffold_ prefix from column 1

sed -i 's,HiC_scaffold_\(.*\),\1,' all_on_narw.bim

## Run KING

~/projects/def-frasiert/RW_WGS/programs/king -b ~/scratch/relatedness/all_on_narw.bed --fam ~/scratch/relatedness/all_on_narw.fam --bim ~/scratch/relatedness/all_on_narw.bim --kinship --degree 3 --prefix all_kinship --sexchr 21
```

Inspecting the results of KING, we confirmed the same relationships we identified before:

Pairs:  EGL140-1 and EGL312-1a
		EGL272-1 and EGL276-1
		Eau019 and Eau10b
		
One additional pair also had a kinship coefficient value which falls just outside KING's bin for 1st degree relationship (0.177-0.354).

Pair:   EGL00252-1 and EGL276-1 : 0.1724 

We will remove the same individuals as with the inslge species analyses.

Parsing data to remove sex scaffolds, related individuals and Linkage Pruning

```
module load mugqic/vcftools/0.1.14 plink/1.9b_6.21-x86_64

vcftools --remove-indv EGL140-1 --remove-indv EGL276-1 --remove-indv Eau019 --not-chr HiC_scaffold_21 --not-chr HiC_scaffold_111 --gzvcf ~/scratch/temp/filtervariants/all_on_narw_masked_filtered.vcf.gz --recode --out all_on_narw_unrelated_filtered


## Linkage pruning
gzip ~/scratch/relatedness/all_on_narw_unrelated_filtered.recode.vcf

plink --vcf ~/scratch/relatedness/all_on_narw_unrelated_filtered.recode.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 \
--out all_on_narw_unrelated_filtered_pruned

plink --vcf ~/scratch/relatedness/all_on_narw_unrelated_filtered.recode.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--extract ~/scratch/relatedness/all_on_narw_unrelated_filtered_pruned.prune.in \
--make-bed --pca --out all_on_narw_unrelated_filtered_pruned

```

Save eigenval and eigenvec files locally for plotting.

### July 20,2022 ###

**Genotype GVCFs**

All scaffolds have finished with the exception of HiC_scaffold_16. It should finish this evening, unless the node speeds up. I will monitor throughout the day. The node appears to be very slow at the moment.

**Stairway Plot**

In the meantime, I can run stairway plot on the folded (unpolarized) VCF for NARW alone (and again for SRW alone.

Stairway plot takes an input file called a 'blueprint' file. They provide an example two-epoch_fold.blueprint, that I have modified for NARW below.

I'm not sure what L should be here. It could be the total length of scaffolds included - 2169635744, it could be the number of variant sites used (even if they are monomorphic, 1234365). I think it is the later due to the wording 'after filtering'. This is what I used below.


```
#NARW blueprint file
#input setting
popid: NARW_alone_folded # id of the population (no white space)
nseq: 22 # number of sequences
L: 1234365 # total number of observed nucleic sites, including polymorphic and monomorphic
whether_folded: true # whether the SFS is folded (true or false)
SFS:    34947 264108 173133 135052 112356 100878  93553  86746	82456  83237  67899 # snp frequency spectrum: number of singleton, number of doubleton, etc. (separated by white space)
#smallest_size_of_SFS_bin_used_for_estimation: 1 # default is 1; to ignore singletons, uncomment this line and change this number to 2
#largest_size_of_SFS_bin_used_for_estimation: 15 # default is nseq/2 for folded SFS
pct_training: 0.67 # percentage of sites for training
nrand: 7        15      22      28 # number of random break points for each try (separated by white space)
project_dir: narw_two-epoch_fold # project directory
stairway_plot_dir: stairway_plot_es # directory to the stairway plot files
ninput: 200 # number of input files to be created for each estimation
#random_seed: 6
#output setting
mu: 1.54e-9 # assumed mutation rate per site per generation
year_per_generation: 35 # assumed generation time (in years) per Taylor 2007 pre whaling estimates
#plot setting
plot_title: NARW_two-epoch_fold # title of the plot
xrange: 0.1,10000 # Time (1k year) range; format: xmin,xmax; "0,0" for default
yrange: 0,100 # Ne (1k individual) range; format: xmin,xmax; "0,0" for default
xspacing: 2 # X axis spacing
yspacing: 2 # Y axis spacing
fontsize: 12 # Font size
```

First I will need to calculate site frequency spectrum. I will do this with scikit allel, using a vcf that contains only bi-allelic SNPs.

I will use bcftools to filter my vcf to include only bi-allelic SNPs using the following code and then generate the SFS in python scikit allel locally.

```
module load bcftools/1.11

bcftools view -m2 -M2 -v snps -Oz -o narw_unrelated_filtered_biallelic.vcf.gz.recode.vcf narw_unrelated_filtered.vcf
```

To run stairway plot, I did the following:

Created the batch file:
```
java -cp stairway_plot_es Stairbuilder narw_alone.blueprint
```

Create a file to submit a scheduled job for the batch:
```
#!/bin/bash
#SBATCH --job-name=stariway_plot_narw
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=6
#SBATCH --mem=25G
#SBATCH --time=12:00:00

bash narw_alone.blueprint.sh
```

Ran this job ```sbatch run_stairwayplot.sh```.

The plots were empty....

I found the raw data and the population sizes were WAYYY off. Spoke with Josquin and the L is a bit more complicated.

It is in fact the number of polymorphic and monomorphic sites in the genome used after filtering. So this means, we can idenitfy how many sites were masked or removed by filtering.

To start, 292469 SNPs were removed by filtering after masking. This includes removing sites on HiC_scaffold_111. Only 2 HiC_scaffold_111 sites were removed by filters. Only 8 were included before masking. 

Using lengths calculated in R, I wrote the following as a .bed file. I will then intersect with the gff in bedtools and keep only the fasta info. Then I should be able to calculate remainging sites.

```
HiC_scaffold_1	0	104503784
HiC_scaffold_2	0	 82659560
HiC_scaffold_3	0	 73844109
HiC_scaffold_4	0	 32574810
HiC_scaffold_5	0	 84139403
HiC_scaffold_6	0	 83526999
HiC_scaffold_7	0	 56280534
HiC_scaffold_8	0	132874210
HiC_scaffold_9	0	129918987
HiC_scaffold_10	0	 96444604
HiC_scaffold_11	0	 97336298
HiC_scaffold_12	0	 96411063
HiC_scaffold_13	0	161762203
HiC_scaffold_14	0	157489382
HiC_scaffold_15	0	 56491294
HiC_scaffold_16	0	169365006
HiC_scaffold_17	0	101413572
HiC_scaffold_18	0	 75544364
HiC_scaffold_19	0	 84528899
HiC_scaffold_20	0	185427052
HiC_scaffold_21	0	104246847
HiC_scaffold_22	0	  2852764
```
Using bedtools intersect to find overlapping regions
```
module load mugqic/bedtools/2.30.0

bedtools maskfasta -fi ~/projects/def-frasiert/RW_WGS/reference/Eubalaena_glacialis_HiC_min1Mb.fasta -bed narw_masked_gff.bed -fo Eubalaena_glacialis_HiC_min1Mb_masked.fasta

grep -v "^>" Eubalaena_glacialis_HiC_min1Mb_masked.fasta | tr -cd N | wc -c
```

Total sum of bases in repeat regions for NARW is: 1126573292
Total bases filtered = # variants in masked VCF - number of bases in the biallelic, no missing VCF = 1526834 - 1194017 = 332817 

For Stairway plot, the number of bases used will be calculated as:

Total # bases in reference sequence for scaff1-22 - bases in those scaffolds masked - bases filtered (+2 to account for bases filtered on scaff111)

2169635744 - 1126573292 - 332817 + 2

Total observed nucleic sites: 1042729637

I ran stairway plot making a few changes:
	- Substituting the number above for L
	- Corrected the nrand numbers per the manual suggestions 
	- Changed the nseq to 20 as the number of sequences (aka genotypes)
	- SFS then had 1 too many values, removed the first which represented monomorphic sites (34947)

The new blueprint file was:
```
#NARW blueprint file
#input setting
popid: NARW_alone_folded # id of the population (no white space)
nseq: 20 # number of sequences
L: 1042729637 # total number of observed nucleic sites, including polymorphic and monomorphic
whether_folded: true # whether the SFS is folded (true or false)
SFS:   256470	167531	130382	108229	96836	89794	83115	78978	81093	67899  # snp frequenc
#smallest_size_of_SFS_bin_used_for_estimation: 1 # default is 1; to ignore singletons, uncomment t>
#largest_size_of_SFS_bin_used_for_estimation: 15 # default is nseq/2 for folded SFS
pct_training: 0.67 # percentage of sites for training
nrand: 4        9      13      18 # number of random break points for each try (separated by white>project_dir: narw_two-epoch_fold # project directory
stairway_plot_dir: stairway_plot_es # directory to the stairway plot files
ninput: 200 # number of input files to be created for each estimation
#random_seed: 6
#output setting
mu: 1.54e-9 # assumed mutation rate per site per generation
year_per_generation: 35 # assumed generation time (in years) per Taylor 2007 pre whaling estimates
#plot setting
plot_title: NARW_two-epoch_fold # title of the plot
xrange: 0.1,10000 # Time (1k year) range; format: xmin,xmax; "0,0" for default
yrange: 0,100 # Ne (1k individual) range; format: xmin,xmax; "0,0" for default
xspacing: 2 # X axis spacing
yspacing: 2 # Y axis spacing
fontsize: 12 # Font size

```
	
**GenotypeGVCFs and filtering**

The final scaffold finished running. Using the code putlined above, I concatenated the vcf files and overnight I am gzipping them and doing initial filtering.

Tomorrow I will plot the depth to determine max depth to use and finish filtering with bcftools.

### July 21, 2022 ###

**Stairwayplot on SRW alone**

Using lengths calculated in R, I wrote the following as a .bed file. I will then intersect with the gff in bedtools and keep only the fasta info. Then I should be able to calculate remaining sites.

```
HiC_scaffold_1	0	178961223
HiC_scaffold_2	0	166458207
HiC_scaffold_3	0	58059695
HiC_scaffold_4	0	88227396
HiC_scaffold_5	0	103990946
HiC_scaffold_6	0	195279244
HiC_scaffold_7	0	136626744
HiC_scaffold_8	0	122997816
HiC_scaffold_9	0	170687504
HiC_scaffold_10	0	87991585
HiC_scaffold_11	0	77160057
HiC_scaffold_12	0	78958011
HiC_scaffold_13	0	101526481
HiC_scaffold_14	0	87763565
HiC_scaffold_15	0	139896917
HiC_scaffold_16	0	59443939
HiC_scaffold_17	0	34393929
HiC_scaffold_18	0	86839656
HiC_scaffold_19	0	107160415
HiC_scaffold_20	0	112042483
HiC_scaffold_21	0	101845965
```

Using bedtools intersect to find overlapping regions
```
module load mugqic/bedtools/2.30.0

bedtools maskfasta -fi ~/projects/def-frasiert/RW_WGS/reference/RWref_HiC_min1Mb.fasta -bed srw_masked_gff.bed -fo RWref_HiC_min1Mb_masked.fasta

grep -v "^>" RWref_HiC_min1Mb_masked.fasta | tr -cd N | wc -c
```

Total sum of bases in repeat regions for SRW is: 1231666157

I will also need to calculate site frequency spectrum. I will do this with scikit allel, using a vcf that contains only bi-allelic SNPs with no missing data.

I will use bcftools to filter my vcf to include only bi-allelic SNPs using the following code and then generate the SFS in python scikit allel locally.

```
module load bcftools/1.11

bcftools view -m2 -M2 -v snps -g ^miss -Oz -o srw_unrelated_filtered_biallelic_nomiss.vcf.gz srw_unrelated_filtered.vcf.gz.recode.vcf.gz
```

Total bases filtered = # variants in masked VCF - number of sites in the biallelic, no missing VCF = 8320997 - 7126131 = 1194866 

For Stairway plot, the number of bases used will be calculated as:

Total # bases in reference sequence for scaff1-21 - bases in those scaffolds masked - bases filtered

2296311778 - 1231666157 - 1194866 

Total observed nucleic sites: 1063450755

The blueprint files looks as follows:
```
#SRW blueprint file
#input setting
popid: SRW_alone_folded # id of the population (no white space)
nseq: 18 # number of sequences
L: 1063450755 # total number of observed nucleic sites, including polymorphic and monomorphic
whether_folded: true # whether the SFS is folded (true or false)
SFS: 2369932	1214357	797685	613145	512186	454147	420778	404189	213611 # snp frequency spectrum: number of singleton, number of doubleton, etc. (separated by white space)
#smallest_size_of_SFS_bin_used_for_estimation: 1 # default is 1; to ignore singletons, uncommen>#largest_size_of_SFS_bin_used_for_estimation: 15 # default is nseq/2 for folded SFS
pct_training: 0.67 # percentage of sites for training
nrand: 4      8      12     16 # number of random break points for each try (separated by white space)
project_dir: srw_two-epoch_fold # project directory
stairway_plot_dir: stairway_plot_es # directory to the stairway plot files
ninput: 200 # number of input files to be created for each estimation
#random_seed: 6
#output setting
mu: 1.54e-9 # assumed mutation rate per site per generation
year_per_generation: 29 # assumed generation time (in years) per Taylor 2007 pre whaling estimates)
#plot setting
plot_title: SRW_two-epoch_fold # title of the plot
xrange: 0.1,10000 # Time (1k year) range; format: xmin,xmax; "0,0" for default
yrange: 0,100 # Ne (1k individual) range; format: xmin,xmax; "0,0" for default
xspacing: 2 # X axis spacing
yspacing: 2 # Y axis spacing
fontsize: 12 # Font size

```

**ADMIXTURE on 2 species mapped to NARW**

First I don't need to include bowhead in the admixture analysis.

```
nano bowhead # SRR1685383

plink --vcf ~/scratch/relatedness/all_on_narw_unrelated_filtered.recode.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--remove bowhead \
--extract ~/scratch/relatedness/all_on_narw_unrelated_filtered_pruned.prune.in \
--make-bed --pca --out rw_unrelated_filtered_pruned
```


I want to run 10 iterations of admixture for each K of 1-6.

Had to change scaff names to numbers only
```
module load admixture/1.3.0

#!/bin/bash
#SBATCH --job-name=rw_admixture
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=22
#SBATCH --mem=10G
#SBATCH --time=24:00:00

prefix=rw_unrelated_filtered_pruned

for r in {1..10}; do for K in {1..6};
do
        admixture --cv -s ${RANDOM} ${prefix}.bed $K -j20
        mv ${prefix}.${K}.Q out/${prefix}.K${K}r${r}.Q
done; done
```

**Polarize VCF**

I am using the polarize vcf python script provided by Josquin to polarize the variant site calls on the vcf of NARW and SRW mapped to NARW.

```
~/scratch/PY_ENV/polarizeVCF.py --vcf ~/projects/def-frasiert/RW_WGS/vcf/locked/all_on_narw_unrelated_filtered.recode.vcf.gz --keep outgroup --miss 0.8 -r > all_on_narw_polarized_body.vcf
```
This kept throwing an error...

### July 22, 2022 ###

**Stairwayplot**

I finishing calculating SFS and the appropriate L. Details were editted inline above. 

I started jobs at ~10:30am to run stairway plot for both NARW and SRW.
stariway_plot_srw-4001

**Visualizing Admixture**

I copied all of the admixture output files locally.

I ran pong locally in the windows command prompt using the following:
```
prefix=rw_unrelated_filtered_pruned
# Create file map
createQmap(){
local r=$1
local K=$2
awk -v K=$K -v r=$r -v file=${prefix}.K${K}r${r} 'BEGIN{ \
printf("K%dr%d\t%d\t%s.Q\n",K,r,K,file)
}' >> ${prefix}_k1-6_Qfilemap
}
export -f createQmap
for K in {1..6}; do for r in {1..10}; do createQmap $r $K; \
done; done
```

Locally:
```
scp crossman@cedar.computecanada.ca:scratch/admixture/out/rw_* .

pong -m rw_unrelated_filtered_pruned_k1-6_Qfilemap --greedy -s .95
```
It generated a series of plots, but nothing really interesting - as anticipated.

The cross validation of error values will help determine the most appropriate K value.

```
#Generate file with cv error rates from log file

grep "^CV error" ../rw_admixture-39985587.out | sed 's,CV error (K=\(.*\)): \(.*\),\1\t\2,' > rw_admixture.CV.txt
```
The data best fit the models assuming a single population as indicated by the lack of convergence of the admixture plots as well as the lowest CV error at K=2.

**Calculating SFS**

SFS can't handle missing genotypes. There are two options: remove missing data or do projections. For stairway plot, I will remove missing genotypes as these only account for a small portion of the total sites.

Therefore, I will take the biallelic vcf and remove sites with missing genotype calls.

```
module load bcftools/1.11

bcftools view -m2 -M2 -v snps -g ^mis -Oz -o narw_unrelated_filtered_biallelic_nomiss.vcf.gz ~/scratch/admixture/narw_unrelated_filtered.vcf.gz.recode.vcf
```

**Polarize VCF**

The called vcf had 913 sites where the genotpye call was "." - haploid. I will figure out how this happened later (and include it as a filter), but now I will filter this vcf to exclude sites where GT="."

```
bcftools view -e 'FORMAT/GT="."' -Oz -o all_on_narw_unrelated_filtered_nohaps.vcf.gz ~/projects/def-frasiert/RW_WGS/vcf/locked/all_on_narw_unrelated_filtered.recode.vcf.gz
```

Running polarizeVCF.py script
```
source PY_ENV/bin/activate

cd polarize/

~/scratch/PY_ENV/polarizeVCF.py --vcf all_on_narw_unrelated_filtered_nohaps.vcf.gz --keep outgroup --miss 0.8 -r > all_on_narw_polarized_body.vcf

# Need to add header
bcftools view -h all_on_narw_unrelated_filtered_nohaps.vcf.gz > header

bcftools reheader -h header
```

### July 23, 2022 ###

After some problem with python libraries, I finally got things working again.

**MSMC**

I am going to start compiling the input data required for MSMC.

MSMC requires phased vcf data created from the .bam files directly. I am going to start writing the script here and try it out on a few samples for HiC_scaffold_22.

```
# Start python environment and minimize library issues
# 	1. login
# 	2. don't load rw_mods
# 	3. load modules

module unload mugqic/python/3.9.1
module load python
module load scipy-stack
module load samtools/1.12
module load bcftools/1.11

source PY_ENV/bin/activate

# 	4. go to msmc directory

cd msmc

```

Run this from within a python environment:

```
#!/bin/bash
#SBATCH --job-name=bamcaller
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-23
#SBATCH --cpus-per-task=10
#SBATCH --mem=25G
#SBATCH --time=36:00:00

REF=~/projects/def-frasiert/RW_WGS/reference/Eubalaena_glacialis_HiC_min1Mb.fasta

SAMPLEBAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samplebam_list)
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samplelist)

for CHR in {1..22}; do 
	DEPTH=$(samtools depth -r HiC_scaffold_$CHR ~/projects/def-frasiert/RW_WGS/merged_bam/$SAMPLEBAM | awk '{sum += $3} END {print sum / NR}')
	bcftools mpileup -Ou -q 20 -Q 20 -C 50 -d 1000 -r HiC_scaffold_$CHR -f $REF ~/projects/def-frasiert/RW_WGS/merged_bam/$SAMPLEBAM | bcftools call -c -V indels | ~/scratch/PY_ENV/bamCaller.py $DEPTH $SAMPLE.HiC_scaffold_$CHR.bed.gz | gzip -c > $SAMPLE.HiC_scaffold_$CHR.vcf.gz
done
```

This appears to be running properly :).


**Stairway plot**

I want to run stairway plot on the unfolded data of both NARW and SRW. It will be hard to compare the folded and unfolded SRW at this time because the SRW here are mapped to NARW. I could and will polarize SRW alone so I can compare the results, but that will be for another day.

Again, the two main parameters I will need to calculate are L (number of observed sites) and SFS (can be calculated in scikit-allel).

First, in order to calculate SFS, I will use bcftools to subset the vcf to include only NARW (or SRW) samples, only biallelic sites, and site with no missing data.

```
module load bcftools/1.11

bcftools view -S NARW ~/projects/def-frasiert/RW_WGS/vcf/locked/all_on_narw_unrelated_filtered.recode.vcf.gz | bcftools view -m2 -M2 -v snps -g ^miss -Oz -o narw_on_narw_unrelated_filtered_biallelic_nomiss.vcf.gz 

bcftools view -S SRW ~/projects/def-frasiert/RW_WGS/vcf/locked/all_on_narw_unrelated_filtered.recode.vcf.gz | bcftools view -m2 -M2 -v snps -g ^miss -Oz -o srw_on_narw_unrelated_filtered_biallelic_nomiss.vcf.gz

```

*Calculating L*

This is done the dame way as before:
Total number of bases: 2169635744
Total sum of bases in repeat regions for NARW is: 1126573292
Total number of filtered sites: Number of sites in vcf after masking (15538544) - number of sites in final vcf ( NARW: 13523440 , SRW: 13417990 )
	NARW on NARW: 2015104
	SRW on NARW: 2120554


For Stairway plot, the number of bases used will be calculated as:

Total # bases in reference sequence for scaff1-22 - bases in those scaffolds masked - bases filtered (+ 3 to account for bases filtered on scaff111)

Total observed nucleic sites: 
	NARW on NARW: 1041047351
	SRW on NARW: 1040941901

Tomorrow I will calculate SFS in scikit allel from the files: narw_on_narw_unrelated_filtered_biallelic_nomiss.vcf.gz and srw_on_narw_unrelated_filtered_biallelic_nomiss.vcf.gz

### July 24, 2022 ###


**Stairway plot**

Calculated SFS in scikit allel on unfolded vcf files. I copied them below into the blueprint files, but it says that the SFS has too many values. The program expects nseq-1.

Created blueprint files.

```
#NARW blueprint file
#input setting
popid: NARW # id of the population (no white space)
nseq: 20 # number of sequences
L: 1041047351 # total number of observed nucleic sites, including polymorphic and monomorphic
whether_folded: false # whethr the SFS is folded (true or false)
SFS:  247889   151323	111304	88963	74622	64934	55653	51227	47221	62232	31816	26071	25821	23092	20359	17215	16633	13684	13715	8918 # snp frequency spectrum: number of singleton, number of doubleton, etc. (separated by white space)
#smallest_size_of_SFS_bin_used_for_estimation: 1 # default is 1; to ignore singletons, uncomment this line and change this number to 2
#largest_size_of_SFS_bin_used_for_estimation: 29 # default is n-1; to ignore singletons, uncomment this line and change this number to nseq-2
pct_training: 0.67 # percentage of sites for training
nrand: 4    9     13     18 # number of random break points for each try (separated by white space)
project_dir: NARW_unfolded # project directory
stairway_plot_dir: stairway_plot_es # directory to the stairway plot files
ninput: 200 # number of input files to be created for each estimation
#random_seed: 6
#output setting
mu: 1.54e-9 # assumed mutation rate per site per generation
year_per_generation: 35 # assumed generation time (in years)
#plot setting
plot_title: NARW_unfolded_two-epoch # title of the plot
xrange: 0.1,10000 # Time (1k year) range; format: xmin,xmax; "0,0" for default
yrange: 0,0 # Ne (1k individual) range; format: xmin,xmax; "0,0" for default
xspacing: 2 # X axis spacing
yspacing: 2 # Y axis spacing
fontsize: 12 # Font size
```

```
#SRW blueprint file
#input setting
popid: SRW # id of the population (no white space)
nseq: 18 # number of sequences
L: 1040941901 # total number of observed nucleic sites, including polymorphic and monomorphic
whether_folded: false # whethr the SFS is folded (true or false)
SFS: 2210056	1071748	663746	479469	373678	308830	264730	234941	222682	177420	164521	153477	144275	138502	135219	136205	143626	761575 # snp frequency spectrum: number of singleton, number of doubleton, etc. (separated by white space)
#smallest_size_of_SFS_bin_used_for_estimation: 1 # default is 1; to ignore singletons, uncomment this line and change this number to 2
#largest_size_of_SFS_bin_used_for_estimation: 29 # default is n-1; to ignore singletons, uncomment this line and change this number to nseq-2
pct_training: 0.67 # percentage of sites for training
nrand: 4    8     12     16 # number of random break points for each try (separated by white space)
project_dir: SRW_unfolded # project directory
stairway_plot_dir: stairway_plot_es # directory to the stairway plot files
ninput: 200 # number of input files to be created for each estimation
#random_seed: 6
#output setting
mu: 1.54e-9 # assumed mutation rate per site per generation
year_per_generation: 29 # assumed generation time (in years)
#plot setting
plot_title: SRW_unfolded_two-epoch # title of the plot
xrange: 0.1,10000 # Time (1k year) range; format: xmin,xmax; "0,0" for default
yrange: 0,0 # Ne (1k individual) range; format: xmin,xmax; "0,0" for default
xspacing: 2 # X axis spacing
yspacing: 2 # Y axis spacing
fontsize: 12 # Font size
```

Created batch files:
```
java -cp stairway_plot_es Stairbuilder narw_unfolded.blueprint
java -cp stairway_plot_es Stairbuilder srw_unfolded.blueprint
```


### July 25, 2022 ###

**Continue data prep for MSMC**

Per the Statistical population genomics book, and the script provided by Josquin.

	- Removes multi-allelic sites in the VCF
	- Makes a list of sites to be excluded in the main run for phasing by shapeit
	- Runs shapeit with --exclude-snp and -nomcmc generating two output files 
	- These two files can be converted into VCF format by shapeit -convert. 
	- Merge the phased VCF sample1.chr$CHR.onlyPhased.vcf.gz and the unphased (original) VCF sample1.chr$CHR.fixedformat.vcf.gz, keeping all unphased sites from the original VCF, but replacing the phased ones.

```
module load bcftools/1.11

#!/bin/bash
#SBATCH --job-name=phasing_vcf
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-23
#SBATCH --cpus-per-task=21
#SBATCH --mem=50G
#SBATCH --time=48:00:00


SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samplelist)

for CHR in {1..22}; do
    UNPHASED_VCF=$SAMPLE.HiC_scaffold_$CHR.vcf.gz
    UNPHASED_VCF_NOMAS=$SAMPLE.HiC_scaffold_$CHR.noMultiAllelicSites.vcf.gz
    LOG_ALIGN=$SAMPLE.HiC_scaffold_$CHR.alignments
    LOG_MAIN=$SAMPLE.HiC_scaffold_$CHR.main
    
    PHASED_HAPS=$SAMPLE.HiC_scaffold_$CHR.phased.haps.gz
    PHASED_SAMPLE=$SAMPLE.HiC_scaffold_$CHR.phased.samples
    PHASED_VCF=$SAMPLE.HiC_scaffold_$CHR.onlyPhased.vcf
    
    LOG_CONVERT=$SAMPLE.HiC_scaffold_$CHR.convert
    FINAL_VCF=$SAMPLE.HiC_scaffold_$CHR.phased.vcf.gz
    
    
	#Preparation
	bcftools view -M 2 -O z $UNPHASED_VCF > $UNPHASED_VCF_NOMAS
	shapeit -check --input-vcf $UNPHASED_VCF_NOMAS --output-log $LOG_ALIGN --thread 20

	#Main run
	shapeit -V $UNPHASED_VCF_NOMAS --output-max $PHASED_HAPS $PHASED_SAMPLE --no-mcmc --output-log $LOG_MAIN --thread 20 --window 0.5
	
    shapeit -convert --input-haps $PHASED_HAPS $PHASED_SAMPLE --output-vcf $PHASED_VCF --output-log $LOG_CONVERT --thread 20

	#Zipping and indexing
	bcftools view -Oz $PHASED_VCF > $PHASED_VCF.gz
	bcftools index -f $PHASED_VCF

	#Merging phased and unphased vcfs, keeping all unphased sites from the original vcf, but replacing the phased ones.
	bcftools merge --force-samples $UNPHASED_VCF $PHASED_VCF | awk 'BEGIN {ofs=”\t”}
        $0 ~ /^#CHROM/ {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}
        $0 !~ /^#/ {
            if(substr($11, 1, 3) != "./.")
                $10 = $11
            print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10
        }' | bcftools view -O z > $FINAL_VCF
done
```

ERROR: You cannot phase less than 10 samples wihout using a reference panel!


Going to merge VCF files across all 23 samples, or actoss each species.
Then try phasing the single vcf file for each chromosome.

I am changing the script above to not run a for loop, but instead run on an array with a job for each chromosome. It will then be easy to test on a CHR22 alone. 

First I will merge the vcf files for CHR22.
```
bcftools merge -Oz -o narwmerged_HiC_scaffold_22_msmc.vcf.gz --file-list chr22
```

There was a problem in that the vcf.gz files don't appear to be compressed properly. I will recompress and reindex them.

```
#!/bin/bash

for sample in EGL00252-1.HiC_scaffold_22 EGL276-1.HiC_scaffold_22 Eau019.HiC_scaffold_22 EGL013-3qa.HiC_scaffold_22 EGL308-1a.HiC_scaffold_22 Eau023.HiC_scaffold_22 EGL140-1.HiC_scaffold_22 EGL312-1a.HiC_scaffold_22 Eau029.HiC_scaffold_22 EGL183-1.HiC_scaffold_22 EGL336_1b.HiC_scaffold_22 Eau034A.HiC_scaffold_22 EGL254-1.HiC_scaffold_22 Eau017.HiC_scaffold_22 Eau10b.HiC_scaffold_22 EGL272-1.HiC_scaffold_22 Eau018.HiC_scaffold_22 Eau283.HiC_scaffold_22 Eau7.HiC_scaffold_22 Eau9c.HiC_scaffold_22 SID179132.HiC_scaffold_22 SID181803.HiC_scaffold_22 SRR1685383.HiC_scaffold_22; do
	mv $sample.vcf.gz $sample.vcf
	bcftools view -Oz -o $sample.vcf.gz $sample.vcf
	bcftools index $sample.vcf.gz
done
```

```
bcftools merge -Oz -o merged_HiC_scaffold22_msmc.vcf.gz --file-list chr22
```




If the shapeit works on this, I will then redo this step for each scaffold.

```
#!/bin/bash
#SBATCH --job-name=phasing_vcf
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-23
#SBATCH --cpus-per-task=21
#SBATCH --mem=50G
#SBATCH --time=48:00:00


CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samplelist)


UNPHASED_VCF=HiC_scaffold_${CHR}_all_on_narw_filtered.vcf.gz
UNPHASED_VCF_NOMAS=HiC_scaffold_${CHR}.noMultiAllelicSites.vcf.gz
LOG_ALIGN=HiC_scaffold_${CHR}.alignments
LOG_MAIN=HiC_scaffold_${CHR}.main
    
PHASED_HAPS=HiC_scaffold_${CHR}.phased.haps.gz
PHASED_SAMPLE=HiC_scaffold_${CHR}.phased.samples
PHASED_VCF=HiC_scaffold_${CHR}.onlyPhased.vcf
    
LOG_CONVERT=HiC_scaffold_${CHR}.convert
FINAL_VCF=HiC_scaffold_${CHR}.phased.vcf.gz
    
    
#Preparation
bcftools view -M 2 -O z $UNPHASED_VCF > $UNPHASED_VCF_NOMAS
shapeit -check --input-vcf $UNPHASED_VCF_NOMAS --output-log $LOG_ALIGN --thread 20

#Main run
shapeit -V $UNPHASED_VCF_NOMAS --output-max $PHASED_HAPS $PHASED_SAMPLE --output-log $LOG_MAIN --thread 20
	
shapeit -convert --input-haps $PHASED_HAPS $PHASED_SAMPLE --output-vcf $PHASED_VCF --output-log $LOG_CONVERT --thread 20

#Zipping and indexing
bcftools view -Oz $PHASED_VCF > $PHASED_VCF.gz
bcftools index -f $PHASED_VCF

#Merging phased and unphased vcfs, keeping all unphased sites from the original vcf, but replacing the phased ones.
bcftools merge --force-samples $UNPHASED_VCF $PHASED_VCF | awk 'BEGIN {ofs=”\t”}
    $0 ~ /^#CHROM/ {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}
    $0 !~ /^#/ {
        if(substr($11, 1, 3) != "./.")
           $10 = $11
        print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10
	}' | bcftools view -O z > $FINAL_VCF

```

ERROR: You cannot phase less than 10 samples wihout using a reference panel!
So merged vcf, but now... TOO MUCH MISSING DATA IN MERGED VCF
	Try merging NARW and SRW southern

Stairway plot age different should be ok, becuase corrected for time and gen on top.

---------------------------------

*Trying vcfAllSiteParser.py*

OK, running the msmc prep pipeline as recommended from the bam files encounters a couple problems. 1) SHAPEIT requires more than 10 samples to do statistica phasing without a map file and the pipeline tries to call a single sample vcf and 2) merging the vcf produced by this pipeline to run in shapeit creates a VCF with a lots of missing data.

MSMC requires a phased vcf and a masking file. There is a python script called vcfAllSiteParser, that will generate the appropriate masking file for MSMC. 

- I think I can use my gvcf files to create part 1 of a relevant masking file and then make a bed file of vcf sites that were and then run SHAPEIT on my nice filtered vcf file. This should provide me with the two file necessary to run MSMC. 

Masking 2 things: 
	1. sites not mapped (could be identified by the gvcf files) (Use vcfAllSiteParser.py to create mask bed file, ignore the vcf it creates)
	2. sites that are filtered (difference between vcf file for each scaffold and final vcf).
			bedtools subtract ~/projects/def-frasiert/RW_WGS/vcf/ALL/merged_all_on_narw.vcf.gz
			~/projects/def-frasiert/RW_WGS/vcf/locked/all_on_narw_unrelated_filtered.vcf.gz
			BED file of sites removed
	3. Create the real mask for each individual with bedtools subtract first vcf filteredvariants.bed

The recommended pipeline involves using the bam files to call variants, then phase them. In this process, a mask file is generated. The Statistical Genomics book recommends doing this for each sample at each chromosome separately, but without a map file, SHAPEIT2 cannot call variants on a single sample. I tried to merge the vcfs that were created from the bam files, both across all samples and across the NARW samples, and in both cases, there was an error noting too much missing data. You can force the program to work and ignore the error or when merging tell the program that missing data should be called as the reference but I don't think we want to do either of these as the important part here is really knowing what is a monomorphic site vs missing data. 

The speciation genomics workshop website goes through phasing, but not MSMC, so they can use any vcf. Before we go to read base phasing, I had another idea.

The generatemultihetsep.py script sets up the input file for MSMC by merging mask.bed files and phased .vcf files across individuals. The mask files help distinguish sites that were unable to be called (either un mapped or filtered), from homozygous sites. I think we can break this down into two parts:
	1. We can use phase our good filtered vcf files using SHAPEIT2
		- We can then split this vcf by individual if require for the generatemultihetsep.py script.
	2. We can generate our masking file
		- As part of the MSMC tools there is a python script called vcfAllSiteParser, that will generate the appropriate masking file for MSMC from a gvcf file. 
		- This however will not account for filtered variants if we are using our filtered vcf. So we would need to also create a bed file for the variants that we filtered out (bedtools subtract VCF from joint Genotyping before filters and VCF after fitlers).
		- We could then combine these two bed files with bedtools to create the appropriate mask file for each individual (the gvcf mask would be unique for each individual and the variant masking would be the same for all individuals).
		
		
**StairwayPlot on unfolded data**

The last value in the SFS represents sites that are fixed for the alternate allele. Therefore, we don't include those in our blueprint file. I am going to rerun the stairway plot code from above, but just remove the final SFS value. 

Going forward, I am going to try a few different combinations:
	- Talk with Tim about the appropriate generation time
	- NARW folded and unfolded
	- SRW folded alone, folded on NARW, and unfolded on SRW and unfolded on NARW



**SHAPEIT**

- split merged/filtered vcf into scaffolds
```
for scaffold in HiC_scaffold_1 HiC_scaffold_2 HiC_scaffold_3 HiC_scaffold_4 HiC_scaffold_5 HiC_scaffold_6 HiC_scaffold_7 HiC_scaffold_8 HiC_scaffold_9 HiC_scaffold_10 HiC_scaffold_11 HiC_scaffold_12 HiC_scaffold_13 HiC_scaffold_14 HiC_scaffold_15 HiC_scaffold_16 HiC_scaffold_17 HiC_scaffold_18 HiC_scaffold_19 HiC_scaffold_20 HiC_scaffold_22; do
	bcftools view -Oz -r $scaffold -o ~/scratch/msmc/shapeit/${scaffold}_all_on_narw_filtered.vcf.gz ~/projects/def-frasiert/RW_WGS/vcf/locked/all_on_narw_unrelated_filtered.rezip.vcf.gz
done
```
At the same time, I should have indexed, so I will run the following:
```
for scaffold in HiC_scaffold_1 HiC_scaffold_2 HiC_scaffold_3 HiC_scaffold_4 HiC_scaffold_5 HiC_scaffold_6 HiC_scaffold_7 HiC_scaffold_8 HiC_scaffold_9 HiC_scaffold_10 HiC_scaffold_11 HiC_scaffold_12 HiC_scaffold_13 HiC_scaffold_14 HiC_scaffold_15 HiC_scaffold_16 HiC_scaffold_17 HiC_scaffold_18 HiC_scaffold_19 HiC_scaffold_20 HiC_scaffold_22; do
	bcftools index ~/scratch/msmc/shapeit/${scaffold}_all_on_narw_filtered.vcf.gz 
done
```


For each multisample scaffold run shapeit script
```
#!/bin/bash
#SBATCH --job-name=phasing_vcf
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-22
#SBATCH --cpus-per-task=21
#SBATCH --mem=50G
#SBATCH --time=12:00:00


CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scaffold_list)

#Set up Input Files
UNPHASED_VCF=HiC_scaffold_${CHR}_all_on_narw_filtered.vcf.gz
UNPHASED_VCF_NOMAS=HiC_scaffold_${CHR}.noMultiAllelicSites.vcf.gz
LOG_ALIGN=HiC_scaffold_${CHR}.alignments
LOG_MAIN=HiC_scaffold_${CHR}.main

PHASED_HAPS=HiC_scaffold_${CHR}.phased.haps.gz
PHASED_SAMPLE=HiC_scaffold_${CHR}.phased.samples
PHASED_VCF=HiC_scaffold_${CHR}.onlyPhased.vcf
NEW_HEAD=HiC_scaffold_${CHR}_newhead.onlyPhased.vcf

LOG_CONVERT=HiC_scaffold_${CHR}.convert
FINAL_VCF=HiC_scaffold_${CHR}.phased.vcf


#Preparation
bcftools view -M 2 -O z $UNPHASED_VCF > $UNPHASED_VCF_NOMAS
shapeit -check --input-vcf $UNPHASED_VCF_NOMAS --output-log $LOG_ALIGN --thread 20

#Main run with force arguement
shapeit -V $UNPHASED_VCF_NOMAS --output-max $PHASED_HAPS $PHASED_SAMPLE --output-log $LOG_MAIN --thread 20 --force

shapeit -convert --input-haps $PHASED_HAPS $PHASED_SAMPLE --output-vcf $PHASED_VCF --output-log $LOG_CONVERT --thread 20

echo "SHAPEIT DONE"

#new header on phased vcf
bcftools view -h $UNPHASED_VCF | bcftools reheader -h - $PHASED_VCF -o ${NEW_HEAD}
bcftools view -Oz -o $NEW_HEAD.gz $NEW_HEAD
bcftools index ${NEW_HEAD}.gz
echo "REHEADER DONE"


#Extracting sites in the unphased vcf, not present in the phased sites (i.e. those sites which could not be phased)

bcftools isec -n -1 -p ${CHR}/ -c all $UNPHASED_VCF $NEW_HEAD.gz
echo "ISEC DONE"

#Concatenate unphased sites with phased sites and sort the vcf.
bcftools concat $NEW_HEAD ${CHR}/0000.vcf -Oz -o $FINAL_VCF
echo "CONCAT DONE"
bcftools sort -Oz -o HiC_scaffold_${CHR}_sorted.vcf.gz $FINAL_VCF

mv HiC_scaffold_${CHR}_sorted.vcf.gz phased/HiC_scaffold_${CHR}_sorted.vcf.gz
```

### July 26, 2022 ###

**SHAPEIT for Phasing**

I finishing the script above and ran it sucessfully on all 22 chromosomes. This generated 22 vcf files. The main thing to note is that I used -force in the main SHAPEIT command to override the missing data warnings.

**Creating bed file for masking for MSMC** 

To generate the MSMC input file, the python script required a bed file "which gives the regions on the chromosome on which the genome of that individual was covered sufficiently"

Creating our own masking file was proving ot be a bit difficult. We can generate our masking file:
	- As part of the MSMC tools there is a python script called vcfAllSiteParser, that will generate the appropriate masking file for MSMC from a gvcf file. 
	- This however will not account for filtered variants if we are using our filtered vcf. So we would need to also create a bed file for the variants that we filtered out (bedtools subtract VCF from joint Genotyping before filters and VCF after filters).
	- We could then combine these two bed files with bedtools to create the appropriate mask file for each individual (the gvcf mask would be unique for each individual and the variant masking would be the same for all individuals).


First I will create the bed files for sites that were able to be mapped. I am doing this using the pythonpackage *gvcf2bed* from here: <https://pypi.org/project/gvcf2bed/>.
```
#!/bin/bash
#SBATCH --job-name=gvcf_bed_narw
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-22
#SBATCH --cpus-per-task=15
#SBATCH --mem=20G
#SBATCH --time=4:00:00

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scaffold_list)

for SAMPLE in EGL00252-1 EGL013-3qa EGL308-1a EGL312-1a EGL183-1 EGL336_1b EGL254-1 EGL272-1 SID179132 SID181803; do \
        gvcf2bed -I ~/scratch/gvcf/NARW/HiC_scaffold_${CHR}/${SAMPLE}-HiC_scaffold_${CHR}.g.vcf.gz -O ${SAMPLE}_${CHR}.bed;
done
```

Next I need to create a bed file that identifies regions of the genome that we masked and/or filtered. I will do this using bedtools subtract to variants not identified in the filtered VCF, but that are listed in the merged vcf. This will output a bed file of variant positions we discarded for a number of reasons. 
```
#!/bin/bash
#SBATCH --job-name=gvcf_bed_narw
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-21
#SBATCH --cpus-per-task=15
#SBATCH --mem=20G
#SBATCH --time=4:00:00

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/scratch/msmc/createBed/scaffold_list)

for SAMPLE in EGL00252-1 EGL013-3qa EGL308-1a EGL312-1a EGL183-1 EGL336_1b EGL254-1 EGL272-1 SID179132 SID181803; do \
	vcftools --gzvcf ~/scratch/msmc/shapeit/phased/HiC_scaffold_${CHR}_sorted.vcf.gz  --indv ${SAMPLE} --recode --out ${SAMPLE}_${CHR}.filtered
	vcftools --gzvcf ~/projects/def-frasiert/RW_WGS/vcf/ALL/merged_all_on_narw.vcf.gz --indv ${SAMPLE} --chr HiC_scaffold_${CHR} --recode --out ${SAMPLE}_${CHR}.originalmerge 
	bedtools subtract -a ${SAMPLE}_${CHR}.originalmerge.recode.vcf -b ${SAMPLE}_${CHR}.filtered.recode.vcf > ${SAMPLE}_${CHR}_filtered_sites.vcf
	rm ${SAMPLE}_${CHR}.originalmerge.recode.vcf;
done
```

Finally, I correct the header on the vcf of filtered sites, and remove these filtered sites from the original bed file. The new bed file represents the
```
#!/bin/bash
#SBATCH --job-name=gvcf_bed_narw
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=3-21
#SBATCH --cpus-per-task=15
#SBATCH --mem=20G
#SBATCH --time=4:00:00

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/scratch/msmc/createBed/scaffold_list)


for SAMPLE in EGL00252-1 EGL013-3qa EGL308-1a EGL312-1a EGL183-1 EGL336_1b EGL254-1 EGL272-1 SID179132 SID181803; do \
	bcftools view -h ~/scratch/msmc/${SAMPLE}.HiC_scaffold_${CHR}.vcf.gz > ~/scratch/msmc/createBed/mask/temp-${SAMPLE}-${CHR}
	
	cat temp-${SAMPLE}-${CHR} ~/scratch/msmc/createBed/${SAMPLE}_${CHR}_filtered_sites.vcf > ${SAMPLE}_${CHR}_filtered_sites_wheader.vcf

	rm temp-${SAMPLE}-${CHR}

# combine bedtools to remove sites filtered from original masking file	
	bedtools subtract -a ~/scratch/msmc/createBed/${SAMPLE}_${CHR}.bed -b ~/scratch/msmc/createBed/mask/${SAMPLE}_${CHR}_filtered_sites_wheader.vcf > ~/scratch/msmc/createBed/mask/${SAMPLE}_${CHR}_msmcmask.bed;
done
```


### July 27, 2022 ###

**MSMC file set up continued**

Now we need to create the input for the generate_multihetsep.py file. The code itself is quite simple, but I will write some short 1 liners that will enable me to copy and paste long lists of file names.
```
# zip bed files
gzip *.bed

# List of the zipped bed files
ls -lf *.bed.gz | sed 's,\(.*\),--mask=createBed/mask/\1 \\,' > zipped_bed

```
We need to have the a separate vcf for each ind at each chromosome.
```

# list of the phased vcf files 
ls -lf *_phased.vcf.gz | sed 's,\(.*\),shapeit/phased/\1 \\,' > phased_vcf
# removed SRW and bowhead form this list 
grep -v "^shapeit/phased/Eau*" phased_vcf | grep -v "shapeit/phased/SRR*" > narw_phased_vcf
```

It recommends creating a mappability mask file. I don't think we need this as the way we generated our mask file accounts for whether or not we could map the sites for each individual.


PHASED VCF in msmc/shapeit/phased/
MASKING BED in createBed/mask/

TO BE RUN FROM MSMC directory
```
#!/bin/bash

./generate_multihetsep.py --mask=covered_sites_sample1_chr1.bed.txt.gz \
                          --mask=covered_sites_sample2_chr1.bed.txt.gz \
                          --mask=mappability_mask_chr1.bed.txt.gz \
                          phased/HiC_scaffold_${CHR}_sorted.vcf.gz \
						  phased/HiC_scaffold_${CHR}_sorted.vcf.gz
```


At this point, I got an error that there can be no missing data in my vcfs (ValueError: invalid literal for int() with base 10: '.').

When I created the list of the filtered sites, I created a vcf per sample and per scaffold (which I deleted, but should have used instead of having to create them again). In creating these vcf with vcftools, I should have specified --max-missing 1 which would indicate no missing data allowed. Then when I generate the mask files, these sites would be masked for the respective individual.

Therefore, I will go back to this step.  

```
#!/bin/bash
#SBATCH --job-name=gvcf_bed_narw
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-21
#SBATCH --cpus-per-task=15
#SBATCH --mem=20G
#SBATCH --time=4:00:00

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/scratch/msmc/createBed/scaffold_list)

for SAMPLE in EGL00252-1 EGL013-3qa EGL308-1a EGL312-1a EGL183-1 EGL336_1b EGL254-1 EGL272-1 SID179132 SID181803; do \
	vcftools --gzvcf ~/scratch/msmc/shapeit/phased/HiC_scaffold_${CHR}_sorted.vcf.gz  --indv ${SAMPLE} --max-missing 1.0 --recode --out ${SAMPLE}_${CHR}.filtered
	vcftools --gzvcf ~/projects/def-frasiert/RW_WGS/vcf/ALL/merged_all_on_narw.vcf.gz --indv ${SAMPLE} --chr HiC_scaffold_${CHR} --recode --out ${SAMPLE}_${CHR}.originalmerge 
	bedtools subtract -a ${SAMPLE}_${CHR}.originalmerge.recode.vcf -b ${SAMPLE}_${CHR}.filtered.recode.vcf > ${SAMPLE}_${CHR}_filtered_sites.vcf
	rm ${SAMPLE}_${CHR}.originalmerge.recode.vcf;

	bcftools view -h ~/scratch/msmc/${SAMPLE}.HiC_scaffold_${CHR}.vcf.gz > ~/scratch/msmc/createBed/mask/temp-${SAMPLE}-${CHR}
	
	cat ~/scratch/msmc/createBed/mask/temp-${SAMPLE}-${CHR} ~/scratch/msmc/createBed/${SAMPLE}_${CHR}_filtered_sites.vcf > ~/scratch/msmc/createBed/mask/${SAMPLE}_${CHR}_filtered_sites_wheader.vcf

	rm ~/scratch/msmc/createBed/mask/temp-${SAMPLE}-${CHR}

#combine bedtools to remove sites filtered from original masking file	
	bedtools subtract -a ~/scratch/msmc/createBed/${SAMPLE}_${CHR}.bed -b ~/scratch/msmc/createBed/mask/${SAMPLE}_${CHR}_filtered_sites_wheader.vcf > ~/scratch/msmc/createBed/mask/${SAMPLE}_${CHR}_msmcmask.bed
	gzip ~/scratch/msmc/createBed/mask/${SAMPLE}_${CHR}_msmcmask.bed;
done
```

This code didn't run in all parts. So I ran the script in chunks to ensure it ran properly.

Creating the script file by the following three one-liners
```
# List of the zipped bed files
ls -lf *.bed.gz | sed 's,\(.*\),--mask=createBed/mask/\1 \\,' > zipped_bed

```
We need to have the a separate vcf for each ind at each chromosome.
```
# I didn't write the code, but I used bcftools to zip the recode.vcf files as the next step wanted them zipped.
# list of the phased vcf files 
ls -lf *filtered.recode.vcf.gz | sed 's,\(.*\),createBed/mask/\1 \\,' > phased_vcf
```

It recommends creating a mappability mask file. I don't think we need this as the way we generated our mask file accounts for whether or not we could map the sites for each individual.


TO BE RUN FROM MSMC directory. The sample set up of the file is below:
```
#!/bin/bash

./generate_multihetsep.py --mask=covered_sites_sample1_chr1.bed.txt.gz \
                          --mask=covered_sites_sample2_chr1.bed.txt.gz \
                          --mask=mappability_mask_chr1.bed.txt.gz \
                          phased/HiC_scaffold_${CHR}_sorted.vcf.gz \
						  phased/HiC_scaffold_${CHR}_sorted.vcf.gz
```

I made mine by adding the #!/bin/bash and the path the the python script manually after the following:
```
cat ~/scratch/msmc/createBed/mask/zipped_bed ~/scratch/msmc/createBed/mask/phased_vcf > generate_msmc_input.sh
```

**IBDSeq**

Another program we will use for historical inference is IBDNe.
IBDNe required a file of identity by descent blocks that can be created with IBDSeq.

IBDSeq is a java program I installed in the RW_WGS/programs/IBDSeq/.

It has two required arguements, gt={vcf file} and out={out prefix}. nthreads={#} can be specified for parallel computing.

```
java -Xmx7G -jar ~/projects/def-frasiert/RW_WGS/programs/IBDseq/ibdseq.r1206.jar \
	gt=~/projects/def-frasiert/RW_WGS/vcf/locked/all_on_narw_unrelated_filtered.rezip.vcf.gz \ out=all_on_narw_ibdseq \
	nthreads=7
```

This code only runs on the first Chromosome list unless a CHROM is specified. Will need to run for each scaffold.

This code also threw errors when there was missing data. It also yielded very little IBDsegments and discarded most sites because of the r2max=0.15 parameter. Increasing this increased the data available. 

Firat I will deal with the missing data. I am going to use the NARW and SRW filtered and unrelated vcf files that were mapped to their own reference genomes. Using VCFtools I will be extracting all sites with missing data. 

```
# Ran in the IBDNe/IBDSeq/

vcftools --gzvcf ~/projects/def-frasiert/RW_WGS/vcf/locked/narw_unrelated_filtered.vcf.gz.recode.vcf --max-missing 1.0 --recode --stdout | gzip -c > narw_all_filters_no_missing.vcf.gz

vcftools --gzvcf ~/projects/def-frasiert/RW_WGS/vcf/locked/srw_unrelated_filtered.vcf.gz.recode.vcf --max-missing 1.0 --recode --stdout | gzip -c > srw_all_filters_no_missing.vcf.gz
```

Trying IBDseq again on NARW with no missing data.

```
java -Xmx7G -jar ~/projects/def-frasiert/RW_WGS/programs/IBDseq/ibdseq.r1206.jar \
	gt=~narw_all_filters_no_missing.vcf.gz out=narw_ibdseq_1 nthreads=7 chrom=HiC_scaffold_1
```

It works, but still generates a very small number of IBD segments with teh default r2max.

I am going to set it up to run on each scaffold with three different r2max values: 0.3, 0.5, 0.7. 

```
for R in 0.1 0.3 0.5 0.7 0.9; do \
	for SCAFFOLD in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 22; do \
		java -Xmx8G -jar ~/projects/def-frasiert/RW_WGS/programs/IBDseq/ibdseq.r1206.jar \
			gt=narw_all_filters_no_missing.vcf.gz out=r${R}/narw_ibdseq_scaff${SCAFFOLD}_${R} nthreads=8 chrom=HiC_scaffold_${SCAFFOLD} r2max=${R};
done; done
```

**IBDNe**

IBDNe requires a map file in PLINK format. For each of my vcf files (NARW and SRW), I will generate a map file in PLINK with the following:
```
module load plink/1.9b_6.21-x86_64
plink --vcf ~/scratch/IBDNe/IBDSeq/narw_all_filters_no_missing.vcf.gz --allow-extra-chr --cm --out narw_plink4ibdne
```

Running IBDNe first for r2max=0.5
```
#!/bin/bash
#SBATCH --job-name=IBDNE_NARW
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=20
#SBATCH --mem=50G
#SBATCH --time=24:00:00

cat ~/scratch/IBDNe/IBDSeq/r0.5/*.ibd | java -jar ~/projects/def-frasiert/RW_WGS/programs/IBDNe/ibdne23Apr20.ae9.jar map=narw_plink4ibdne.map out=narw_ibdne_r.5 nthreads=20
```

### July 28, 2022 ###

**MSMC**

No output file was generated by generate_multihetsep.py because I didn't specify and > out.file arguement at the end. I added this and restarted the script. I worry that it might not work properly as there is something being written to the files while it is running, even though something was written to the screen before and then saved in the log file which is also currently empty. 

This ran, but only had 25 IBD identified. I realized, it needed to be run with each scaffold separately.

I ran the following and saved MSMC input files in input_files/
```
#!/bin/bash
#SBATCH --job-name=generate_MSMC_input
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-21
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G
#SBATCH --time=1:00:00

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scaffold_list)

~/scratch/PY_ENV/generate_multihetsep.py \
--mask=createBed/mask/EGL00252-1_${CHR}_msmcmask.bed.gz \
--mask=createBed/mask/EGL013-3qa_${CHR}_msmcmask.bed.gz \
--mask=createBed/mask/EGL183-1_${CHR}_msmcmask.bed.gz \
--mask=createBed/mask/EGL254-1_${CHR}_msmcmask.bed.gz \
--mask=createBed/mask/EGL272-1_${CHR}_msmcmask.bed.gz \
--mask=createBed/mask/EGL308-1a_${CHR}_msmcmask.bed.gz \
--mask=createBed/mask/EGL312-1a_${CHR}_msmcmask.bed.gz \
--mask=createBed/mask/EGL336_1b_${CHR}_msmcmask.bed.gz \
--mask=createBed/mask/SID179132_${CHR}_msmcmask.bed.gz \
--mask=createBed/mask/SID181803_${CHR}_msmcmask.bed.gz \
createBed/mask/EGL00252-1_${CHR}.filtered.recode.vcf.gz \
createBed/mask/EGL013-3qa_${CHR}.filtered.recode.vcf.gz \
createBed/mask/EGL183-1_${CHR}.filtered.recode.vcf.gz \
createBed/mask/EGL254-1_${CHR}.filtered.recode.vcf.gz \
createBed/mask/EGL272-1_${CHR}.filtered.recode.vcf.gz \
createBed/mask/EGL308-1a_${CHR}.filtered.recode.vcf.gz \
createBed/mask/EGL312-1a_${CHR}.filtered.recode.vcf.gz \
createBed/mask/EGL336_1b_${CHR}.filtered.recode.vcf.gz \
createBed/mask/SID179132_${CHR}.filtered.recode.vcf.gz \
createBed/mask/SID181803_${CHR}.filtered.recode.vcf.gz > input_files/narw_${CHR}_multihetsep.txt
```

This ran successfully and it is finally time to start running MSMC.

```
#!/bin/bash
#SBATCH --job-name=MSMC2
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=20
#SBATCH --mem=124G
#SBATCH --time=24:00:00

msmc2 -t 20 -o results/narw_msmc2 \
	input_files/narw XX
	
	LIST FILES HERE
	XXX
	XXX
```


**IBDNe**

IBDNe requires a linkage map. We obviously don't have this for right whales, so we could use a constant recombination rate and try a few different vlues to test the robustness of our results. First we need to generate a linkage map file.

To do this, we can use PLINK, but then we need to change a few values to make it exactly in the format we are looking for (set constant recombination rate, and change to centimorgans).  

```
#generate basic plink map file

module load plink/1.9b_6.21-x86_64

plink --vcf ~/scratch/IBDNe/IBDSeq/narw_all_filters_no_missing.vcf.gz --recode --allow-extra-chr --out narw_4ibdne

#Run awk to add position to second column, set recombination rate and change the fourth column to be the position in cm along

awk '{$2=$4; $3=0.8; $4=$2/1000000*$3; print $0}' narw_4ibdne.map > r0.8/narw_r0.8.map

awk '{$2=$4; $3=1.0; $4=$2/1000000*$3; print $0}' narw_4ibdne.map > r1.0/narw_r1.0.map

awk '{$2=$4; $3=1.2; $4=$2/1000000*$3; print $0}' narw_4ibdne.map > r1.2/narw_r1.2.map

```

Now I can run 15 different combinations of IBDNe 3 recombination values (represented by 3 different map files) and 5 different r2max values from IDBSeq.

```
#!/bin/bash
#SBATCH --job-name=IBDNE_NARW
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-15
#SBATCH --cpus-per-task=20
#SBATCH --mem=20G
#SBATCH --time=24:00:00

INPUTS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ibdne_options)
R2MAX=$(echo $INPUTS | awk '{print $1}')
RECOM=$(echo $INPUTS | awk '{print $2}')

cat ~/scratch/IBDNe/IBDSeq/r${R2MAX}/*.ibd | java -jar ~/projects/def-frasiert/RW_WGS/programs/IBDNe/ibdne.23Apr20.ae9.jar map=~/scratch/IBDNe/r${RECOM}/fixed_narw_r${RECOM}.map out=narw_r2max${R2MAX}_recom${RECOM} nthreads=20

```

### July 29, 2022 ###

**MSMC2**

I ran msmc2 on all haplotypes and it failed with a lack of memory.

I reran the generate input with 4 samples (and therefore 8 haplotypes). We can rerun this with different samples to compare the results.

Repeated generate_multihetsep.py and reran msmc2 with just:
	EGL00252-1
	EGL183-1
	EGL013-3qa
	SID181803

**IBDNe**

The variability in our results with different recombination rates underscored the importance of finding an appropriate recombination rate and not applying a generic constant value.

We can use LDJump (an R package that employs the C/C++ program LDHat) to generate a recombination map file. 


### August 2, 2022 ###

**Start Haplotype Calling for Bowhead and SRW**
First I need to add the bowhead to the SRW genomics DB

*Add to GenomicsDatabases*

If I use a slightly newer version of GATK (v4.1.8 instead of v4.1.0), I can add samples to a genomics database and then genotype GVCFs across all of them.

First I will unload the older verison of GATK and load the newer one.

```
module unload mugqic/GenomeAnalysisTK/4.1.0.0
module load mugqic/GenomeAnalysisTK/4.1.8.1
```

To run through all scaffolds:

```
#!/bin/bash
#SBATCH --job-name=add_bowhead_array
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=2-21
#SBATCH --cpus-per-task=6
#SBATCH --mem=25G
#SBATCH --time=12:00:00

SCAFF=${SLURM_ARRAY_TASK_ID}

FILE=~/scratch/gvcf/bowhead/map_to_SRW/HiC_scaffold_${SCAFF}/SRR1685383-HiC_scaffold_${SCAFF}-SRW.g.vcf.gz

java -Xms20G -Xmx20G -jar ${GATK_JAR} GenomicsDBImport \
        --genomicsdb-update-workspace-path ~/scratch/genomicsDB/SRW_workspace/scaff${SCAFF}/ \
        -V ${FILE} \
        --tmp-dir ~/scratch/temp
```

*Genotype GVCFs*

Ran this script to call genotypes of SRW with bowhead included.
```
#!/bin/bash
#SBATCH --job-name=genotype_gvcf
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=2-21
#SBATCH --cpus-per-task=6
#SBATCH --mem=25G
#SBATCH --time=36:00:00

FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" srw_gendb_inputs)

FOLDER=$(echo $FILE | awk '{print $1}')
SCAFFOLD=$(echo $FILE | awk '{print $2}')

java -Xms20G -Xmx20G -jar ${GATK_JAR} GenotypeGVCFs \
   -R ~/projects/def-frasiert/RW_WGS/reference/RWref_HiC_min1Mb.fasta \
   -V gendb://$FOLDER/ \
   -O ~/scratch/vcf/SRW_w_bowhead_HiC_scaffold_$SCAFFOLD.vcf.gz \
   2> ~/projects/def-frasiert/RW_WGS/QC/gvcfQC/SRW_w_bowhead_HiC_scaffold_$SCAFFOLD-gvcf.out
```

**MSMC2 different set of samples**

I want to run the MSMC2 analysis on a different set of samples for comparison.

```
#!/bin/bash
#SBATCH --job-name=generate_MSMC_input
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-21
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G
#SBATCH --time=1:00:00

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scaffold_list)

~/scratch/PY_ENV/generate_multihetsep.py \
--mask=createBed/mask/EGL254-1_${CHR}_msmcmask.bed.gz \
--mask=createBed/mask/EGL272-1_${CHR}_msmcmask.bed.gz \
--mask=createBed/mask/EGL308-1a_${CHR}_msmcmask.bed.gz \
--mask=createBed/mask/EGL312-1a_${CHR}_msmcmask.bed.gz \
createBed/mask/EGL254-1_${CHR}.filtered.recode.vcf.gz \
createBed/mask/EGL272-1_${CHR}.filtered.recode.vcf.gz \
createBed/mask/EGL308-1a_${CHR}.filtered.recode.vcf.gz \
createBed/mask/EGL312-1a_${CHR}.filtered.recode.vcf.gz \
> input_files/narw_${CHR}_multihetsep_set2.txt
```

### August 3, 2022 ###

Last night I got an email from the system admins at Compute Canada telling me I should write to $SLURM_TMPDIR for temporary file storage. Where possible, I should put this in my code for scheduled jobs!



**Running MSMC2**

I tried running this yesterday, but it required an incredible amount of memory (even quit with 200G). So I can use the same generate MSMC input files and run on just 6 haplotypes. If this works, I think it will be best to generate the multihet file on all samples and then run in iterations of different sample sets.

```
#!/bin/bash
#SBATCH --job-name=MSMC2
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=20
#SBATCH --mem=124G
#SBATCH --time=24:00:00

msmc2 -t 20 -I 0,1,2,3 -o results/narw_msmc2/set2/ \
	input_files/narw XX
	
	LIST FILES HERE
	XXX
	XXX
```
Still didn't work with this code to compare only a subset of haplotypes. I think this still has to read the whole files and that doesn't work. Will need to generate multihet sep again on fewer samples.

**Genotype gvcfs on scaff1**

Ran this script to do chr 1 as it is kind of different.

### August 4, 2022 ###

**SRW with bowhead filtering**

Once I have called genotypes, I need to merge my vcf files using *bcftools concat*
```
module load bcftools/1.11

bcftools concat SRW_w_bowhead_HiC_scaffold_1.vcf.gz SRW_w_bowhead_HiC_scaffold_2.vcf.gz SRW_w_bowhead_HiC_scaffold_3.vcf.gz SRW_w_bowhead_HiC_scaffold_4.vcf.gz SRW_w_bowhead_HiC_scaffold_5.vcf.gz SRW_w_bowhead_HiC_scaffold_6.vcf.gz SRW_w_bowhead_HiC_scaffold_7.vcf.gz SRW_w_bowhead_HiC_scaffold_8.vcf.gz SRW_w_bowhead_HiC_scaffold_9.vcf.gz SRW_w_bowhead_HiC_scaffold_10.vcf.gz SRW_w_bowhead_HiC_scaffold_11.vcf.gz SRW_w_bowhead_HiC_scaffold_12.vcf.gz SRW_w_bowhead_HiC_scaffold_13.vcf.gz SRW_w_bowhead_HiC_scaffold_14.vcf.gz SRW_w_bowhead_HiC_scaffold_15.vcf.gz SRW_w_bowhead_HiC_scaffold_16.vcf.gz SRW_w_bowhead_HiC_scaffold_17.vcf.gz SRW_w_bowhead_HiC_scaffold_18.vcf.gz SRW_w_bowhead_HiC_scaffold_19.vcf.gz SRW_w_bowhead_HiC_scaffold_20.vcf.gz SRW_w_bowhead_HiC_scaffold_21.vcf.gz -Oz -o ~/projects/def-frasiert/RW_WGS/vcf/SRW/merged_SRW_w_bowhead.vcf.gz
```

With merged VCF files, I need to perform some filtering. First I will mask repeat regions defined in the .gff file. Then, I will filters the variants sites.

```
module load mugqic/bedtools/2.30.0 bcftools/1.11 mugqic/vcftools/0.1.14


bedtools intersect -v -a ~/projects/def-frasiert/RW_WGS/vcf/SRW/merged_SRW_w_bowhead.vcf.gz -b ~/scratch/RWref_HiC.repeatmasker.trf.windowmasker.gff.gz -wa -header > ~/scratch/temp/filtervariants/SRW_w_bowhead_masked.vcf

vcftools --vcf ~/scratch/temp/filtervariants/SRW_w_bowhead_masked.vcf --recode-INFO-all -Oz --out ~/scratch/temp/filtervariants/SRW_w_bowhead_masked_formatfilter --minDP 10 --minGQ 30
```

Finish the filtering with appropriate mapping quality and max depth.
2x median value seems to work again. 2x depth = 636
```
bcftools filter -i 'COUNT(FORMAT/GT="mis")<6 && INFO/MQ>30 && INFO/DP<636' -Oz -o ~/scratch/temp/filtervariants/SRW_w_bowhead_masked_filtered.vcf.gz  ~/scratch/temp/filtervariants/SRW_w_bowhead_masked_formatfilter.recode.vcf.gz
```

** EASYSFS **

Ran EasySFS locally and it yielded the exact same SFS as scikit allel. This rules out at least one potential source of discrepancy between the species.

### August 5, 2022 ###


**CHECK FILTERING FROM PREVIOUS TESTS**

It looks like there may have been a mistake in the filtering where the format filters were not included. I am going to rerun the filtering on the SRW called alone, and compare this new vcf to my previous vcf. 

```
vcftools --vcf ~/scratch/temp/filtervariants/srw_masked.vcf --recode --recode-INFO-all --out ~/scratch/temp/filtervariants/srw_masked_formatfilter_aug5 --minDP 10 --minGQ 30
```
This doesn't remove sites, so the depth cutoff should be the same (checked)
```
bcftools filter -i 'COUNT(FORMAT/GT="mis")<3 && INFO/MQ>30 && INFO/DP<558' -Oz -o ~/scratch/temp/filtervariants/srw_masked_filtered_aug5.vcf.gz   ~/scratch/temp/filtervariants/srw_masked_formatfilter_aug5.vcf.gz
```

Used vcf-compare to check and the final filtered vcfs were not the same as I suspected. I will recreated the NARW and all_on_NARW final filtered vcfs from the masked vcf, and repeat any analyses as needed.

The depth thresholds will be the same.

SRW alone done to filtered vcf. Number of SNPs updated.

NARW alone
```
vcftools --gzvcf ~/scratch/temp/filtervariants/narw_masked.vcf.gz --recode --recode-INFO-all --out ~/scratch/temp/filtervariants/narw_masked_formatfilter_aug5 --minDP 10 --minGQ 30

gzip narw_masked_formatfilter_aug5.recode.vcf

bcftools filter -i 'COUNT(FORMAT/GT="mis")<=3 && INFO/MQ>30 && INFO/DP<790' -Oz -o ~/scratch/temp/filtervariants/narw_masked_filtered_aug5.vcf.gz   ~/scratch/temp/filtervariants/narw_masked_formatfilter_aug5.recode.vcf.gz
```

All on NARW
```
vcftools --gzvcf ~/scratch/temp/filtervariants/all_on_narw_masked.vcf.gz --recode --recode-INFO-all --out ~/scratch/temp/filtervariants/all_on_narw_masked_formatfilter_aug5 --minDP 10 --minGQ 30

gzip all_on_narw_masked_formatfilter_aug5.recode.vcf

bcftools filter -i 'COUNT(FORMAT/GT="mis")<6 && INFO/MQ>30 && INFO/DP<1370' -Oz -o ~/scratch/temp/filtervariants/all_on_narw_masked_filtered_aug5.vcf.gz  ~/scratch/temp/filtervariants/all_on_narw_masked_formatfilter_aug5.recode.vcf.gz
```

** Unrelated filters on all fixed vcf files **

*Rerunning KING*
```
## Inside an interactive node

module load plink/1.9b_6.21-x86_64

## From inside ~/scratch/relatedness/
## Run PLINK to generate .bed, .fam and .bim files

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/NARW/narw_masked_filtered_aug5.vcf.gz --allow-extra-chr --out narw

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/SRW/srw_masked_filtered_aug5.vcf.gz --allow-extra-chr --out srw

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/ALL/all_on_narw_masked_filtered_aug5.vcf.gz --allow-extra-chr --out all

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/SRW/SRW_w_bowhead_masked_filtered.vcf.gz --allow-extra-chr --out srw_w_bowhead

## Edit the .bim file to remove the HiC_scaffold_ prefix from column 1 and rename a few to clearly identify sex chromosome


sed 's,HiC_scaffold_\(.*\),\1,' all.bim > all_edit.bim
sed 's,HiC_scaffold_\(.*\),\1,' narw.bim > narw_edit.bim
sed 's,HiC_scaffold_\(.*\),\1,' srw_w_bowhead.bim > srw_w_bowhead_edit.bim
sed -i 's,21\(\t.*\),1\1,' srw_w_bowhead_edit.bim
sed -i 's,8\(\t.*\),21\1,' srw_w_bowhead_edit.bim
sed 's,HiC_scaffold_\(.*\),\1,' srw.bim > srw_edit.bim
sed -i 's,21\(\t.*\),1\1,' srw_edit.bim
sed -i 's,8\(\t.*\),21\1,' srw_edit.bim

## Run KING

~/projects/def-frasiert/RW_WGS/programs/king -b ~/scratch/relatedness/narw.bed --fam ~/scratch/relatedness/narw.fam --bim ~/scratch/relatedness/narw_edit.bim --kinship --degree 3 --prefix narw_kinship --sexchr 21

~/projects/def-frasiert/RW_WGS/programs/king -b ~/scratch/relatedness/all.bed --fam ~/scratch/relatedness/all.fam --bim ~/scratch/relatedness/all_edit.bim --kinship --degree 3 --prefix all_kinship --sexchr 21

~/projects/def-frasiert/RW_WGS/programs/king -b ~/scratch/relatedness/srw.bed --fam ~/scratch/relatedness/srw.fam --bim ~/scratch/relatedness/srw_edit.bim --kinship --degree 5 --prefix srw_kinship --sexchr 21

~/projects/def-frasiert/RW_WGS/programs/king -b ~/scratch/relatedness/srw_w_bowhead.bed --fam ~/scratch/relatedness/srw_w_bowhead.fam --bim ~/scratch/relatedness/srw_w_bowhead_edit.bim --kinship --degree 5 --prefix srw_w_bowhead_kinship --sexchr 21
```


*Removing Individuals and Sex Scaffolds and HiC_scaffold_111*

```
###WROTE THIS INTO A SCRIPT

module load mugqic/vcftools/0.1.14

vcftools --remove-indv EGL312-1a --remove-indv SID181803 --remove-indv EGL276-1 --not-chr HiC_scaffold_21 --not-chr HiC_scaffold_111 --gzvcf ~/projects/def-frasiert/RW_WGS/vcf/NARW/narw_masked_filtered_aug5.vcf.gz --recode --stdout | gzip -c > ~/projects/def-frasiert/RW_WGS/vcf/locked/narw_unrelated_filtered_aug5.vcf.gz

vcftools --remove-indv Eau019 --not-chr HiC_scaffold_8 --gzvcf ~/projects/def-frasiert/RW_WGS/vcf/SRW/srw_masked_filtered_aug5.vcf.gz --recode --stdout | gzip -c > ~/projects/def-frasiert/RW_WGS/vcf/locked/srw_unrelated_filtered_aug5.vcf.gz

vcftools  --remove-indv EGL312-1a --remove-indv SID181803 --remove-indv EGL276-1 --remove-indv Eau019 --not-chr HiC_scaffold_21 --not-chr HiC_scaffold_111 --gzvcf ~/projects/def-frasiert/RW_WGS/vcf/ALL/all_on_narw_masked_filtered_aug5.vcf.gz --recode --stdout | gzip -c > ~/projects/def-frasiert/RW_WGS/vcf/locked/all_on_narw_unrelated_filtered_aug5.vcf.gz

vcftools --remove-indv Eau019 --not-chr HiC_scaffold_8 --gzvcf ~/projects/def-frasiert/RW_WGS/vcf/SRW/SRW_w_bowhead_masked_filtered.vcf.gz --recode --stdout | gzip -c > ~/projects/def-frasiert/RW_WGS/vcf/locked/SRW_w_bowhead_unrelated_filtered_aug5.vcf.gz
```

*Population Structure PCA*
```
module load plink/1.9b_6.21-x86_64

### NARW ###

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/locked/narw_unrelated_filtered_aug5.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.2 \
--out narw_plink

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/locked/narw_unrelated_filtered_aug5.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--extract narw_plink.prune.in \
--make-bed --pca --out narw_plink

### SRW ###

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/locked/srw_unrelated_filtered_aug5.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 \
--out srw_plink

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/locked/srw_unrelated_filtered_aug5.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--extract srw_plink.prune.in \
--make-bed --pca --out srw_plink

### ALL ###

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/locked/all_on_narw_unrelated_filtered_aug5.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.2 \
--out all_on_narw_plink

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/locked/all_on_narw_unrelated_filtered_aug5.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--extract all_on_narw_plink.prune.in \
--make-bed --pca --out all_on_narw_plink
```

Variants removed with Pruning:
NARW - 1313703 of 1407025 variants removed
SRW -   6097123 of 6252485 variants removed
All - 11985339 of 12793691 variants removed

### August 6, 2022 ###
**LDdecay**
```
module load poplddecay/3.41

PopLDdecay    -InVCF  ~/projects/def-frasiert/RW_WGS/vcf/locked/srw_unrelated_filtered_aug5.vcf.gz  -OutStat LDdecay
PopLDdecay    -InVCF  ~/projects/def-frasiert/RW_WGS/vcf/locked/narw_unrelated_filtered_aug5.vcf.gz  -OutStat NARW_LDdecay
PopLDdecay    -InVCF  ~/projects/def-frasiert/RW_WGS/vcf/locked/all_on_narw_unrelated_filtered_aug5.vcf.gz  -OutStat ALL_LDdecay
```


**Admixture**

*PLINK with larger r2 threshold for NARW, r1 for SRW*
```
module load plink/1.9b_6.21-x86_64

# NARW
plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/locked/narw_unrelated_filtered_aug5.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.2 \
--out narw_plink_r.2

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/locked/narw_unrelated_filtered_aug5.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--extract narw_plink_r.2.prune.in \
--make-bed --pca --out narw_plink_r.2

# SRW
plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/locked/srw_unrelated_filtered_aug5.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 \
--out srw_plink_r.1

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/locked/srw_unrelated_filtered_aug5.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--extract srw_plink_r.1.prune.in \
--make-bed --pca --out srw_plink_r.1

#ALL

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/locked/all_on_narw_unrelated_filtered_aug5.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.2 \
--out all_plink_r.2

plink --vcf ~/projects/def-frasiert/RW_WGS/vcf/locked/all_on_narw_unrelated_filtered_aug5.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--extract all_plink_r.2.prune.in \
--make-bed --pca --out all_plink_r.2
```

Had to change scaff names to numbers only using sed on the .bim file



**ADMIXTURE**

I want to run 10 iterations of admixture for each K of 1-5.

*NARW*
```
module load admixture/1.3.0

#!/bin/bash
#SBATCH --job-name=narw_admixture
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=22
#SBATCH --mem=10G
#SBATCH --time=6:00:00

prefix=narw_plink_r.2

cp ${prefix}* ${SLURM_TMPDIR}
cd ${SLURM_TMPDIR}


for r in {1..10}; do for K in {1..5};
do
        admixture --cv -s ${RANDOM} ${prefix}.bed $K -j20
        mv ${prefix}.${K}.Q ${prefix}.K${K}r${r}.Q
done; done

cp ${prefix}* ${SLURM_SUBMIT_DIR}/out
```

*SRW*


```
#!/bin/bash
#SBATCH --job-name=srw_admixture
#SBATCH --output=srw_admixture_log.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=22
#SBATCH --mem=10G
#SBATCH --time=1:00:00

prefix=srw_plink_r.1

cp ${prefix}* ${SLURM_TMPDIR}
cd ${SLURM_TMPDIR}


for r in {1..10}; do for K in {1..5};
do
        admixture --cv -s ${RANDOM} ${prefix}.bed $K -j20
        mv ${prefix}.${K}.Q ${prefix}.K${K}r${r}.Q
done; done

cp ${prefix}* ${SLURM_SUBMIT_DIR}/out
```

*ALL*

```
#!/bin/bash
#SBATCH --job-name=rw_admixture
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=22
#SBATCH --mem=10G
#SBATCH --time=2:00:00

prefix=all_plink_r.2

cp ${prefix}* ${SLURM_TMPDIR}
cd ${SLURM_TMPDIR}


for r in {1..10}; do for K in {1..6};
do
        admixture --cv -s ${RANDOM} ${prefix}.bed $K -j20
        mv ${prefix}.${K}.Q ${prefix}.K${K}r${r}.Q
done; done

cp ${prefix}* ${SLURM_SUBMIT_DIR}/out
```

**Polarize VCF**

First ensure we are not missing any haplotypes

```
bcftools view -e 'FORMAT/GT="."' -Oz -o all_on_narw_unrelated_filtered_nohaps.vcf.gz ~/projects/def-frasiert/RW_WGS/vcf/locked/all_on_narw_unrelated_filtered_aug5.vcf.gz

bcftools view -e 'FORMAT/GT="."' -Oz -o SRW_w_bowhead_unrelated_filtered_nohaps.vcf.gz ~/projects/def-frasiert/RW_WGS/vcf/locked/SRW_w_bowhead_unrelated_filtered_aug5.vcf.gz
```

Running polarizeVCF.py script
```
source PY_ENV/bin/activate

cd polarize/

~/scratch/PY_ENV/polarizeVCF.py --vcf all_on_narw_unrelated_filtered_nohaps.vcf.gz --keep outgroup --miss 0.8 -r > all_on_narw_polarized_body.vcf

# Need to add header
bcftools view -h all_on_narw_unrelated_filtered_nohaps.vcf.gz > all_header

bcftools view -h SRW_w_bowhead_unrelated_filtered_nohaps.vcf.gz > srw_header

```

PYTHON 3.8.10


### August 7, 2022 ###

**Running pong**

Now I need to plot the data. First I will plot the admixture plots with pong. I will create a file map and then transfer the files locally to use pong.

```
# Create file map
createQmap(){
local r=$1
local K=$2
awk -v K=$K -v r=$r -v file=${prefix}.K${K}r${r} 'BEGIN{ \
printf("K%dr%d\t%d\t%s.Q\n",K,r,K,file)
}' >> ${prefix}_k1-5_Qfilemap
}
export -f createQmap
for K in {1..5}; do for r in {1..10}; do createQmap $r $K; \
done; done

```

I ran pong locally in the windows command prompt using the following:
```
# Run pong 
pong -m srw_plink_r.1_k1-5_Qfilemap --greedy -s .95
pong -m narw_plink_r.2_k1-5_Qfilemap --greedy -s .95
pong -m all_plink_r.2_k1-6_Qfilemap --greedy -s .95
```
It generated a series of plots, but nothing really interesting - as anticipated.

The cross validation of error values will help determine the most appropriate K value.

```
#Generate file with cv error rates from log file

grep "^CV error" ../narw_admixture-41371918.out | sed 's,CV error (K=\(.*\)): \(.*\),\1\t\2,' > narw_plink_r.2.CV.txt


grep "^CV error" ../srw_admixture_log.out | sed 's,CV error (K=\(.*\)): \(.*\),\1\t\2,' > srw_plink_r.1.CV.txt

grep "^CV error" ../rw_admixture-41427188.out | sed 's,CV error (K=\(.*\)): \(.*\),\1\t\2,' > all_plink_r.2.CV.txt
```

**SFS for Stairway Plot**

I will use scikit allel in python to calculate SFS, but I will need to first just have SRW and NARW included.


```
module load bcftools/1.11

bcftools view -S NARW ~/scratch/polarize/all_on_narw_polarized.vcf.gz | bcftools view -m2 -M2 -v snps -g ^miss -Oz -o narw_polarized_biallelic_nomiss.vcf.gz 

bcftools view -S SRW ~/scratch/polarize/srw_w_bowhead_polarized.vcf.gz | bcftools view -m2 -M2 -v snps -g ^miss -Oz -o srw_polarized_biallelic_nomiss.vcf.gz

```

FROM SCIKIT ALLEL:

NARW - 3431049 113590 70218 49824 40587 36325 30455 28591 25640 29375 22821 23079 21706 22868 22065 21690 23493 25637 3233246

SRW - 63606 1043286 487379 297904 213855 168420 138759 119309 105206 96732 85318 79064 74331   71283 69154 69576 71452 78805 1756277


**Calculating L**

*NARW*
This is done the same way as before -
Total number of bases: 2169635744
Total sum of bases in repeat regions for NARW is: 1126573292
Total number of filtered sites: Number of sites in vcf after masking (15538544) - number of sites in final vcf ( NARW: 7272259 ) = 8266285

For Stairway plot, the number of bases used will be calculated as:

Total # bases in reference sequence for scaff1-22 - bases in those scaffolds masked - bases filtered (+ 3 to account for bases filtered on scaff111)

Total observed nucleic sites for NARW: 1034796170
	
*SRW*
This is done the same way as before -
Total number of bases: 2296311778
Total sum of bases in repeat regions for SRW is: 1231666157
Total number of filtered sites: Number of sites in vcf after masking (14387933) - number of sites in final vcf ( SRW: 5089716 ) = 9298217

For Stairway plot, the number of bases used will be calculated as:

Total # bases in reference sequence for scaff1-21 - bases in those scaffolds masked - bases filtered

Total observed nucleic sites for SRW: 1055347404


**Started Stairway plot**


### August 8, 2022 ###

** Admixture **

Finished running admixture and pong on all_on_narw.
CV error lowest at 2, but 3 also low and PONG converged on all 10 runs on K=1, K=2, K=3.

** ROH **

Ran the python script for Runs of Homozygosity on NARW and then SRW.

** Rerun MSMC on newly filtered vcf files **

First, I need to phase the data. I am phasing all_on_NARW and SRW_w_bowhead and I can extract the relevant samples later as phasing cannot be done for less than 10 samples
```
#prepare input vcf files separate for each scaffold
for scaffold in HiC_scaffold_1 HiC_scaffold_2 HiC_scaffold_3 HiC_scaffold_4 HiC_scaffold_5 HiC_scaffold_6 HiC_scaffold_7 HiC_scaffold_8 HiC_scaffold_9 HiC_scaffold_10 HiC_scaffold_11 HiC_scaffold_12 HiC_scaffold_13 HiC_scaffold_14 HiC_scaffold_15 HiC_scaffold_16 HiC_scaffold_17 HiC_scaffold_18 HiC_scaffold_19 HiC_scaffold_20 HiC_scaffold_22; do
	bcftools view -Oz -r $scaffold -o ~/scratch/msmc/shapeit/${scaffold}_all_on_narw_filtered_aug8.vcf.gz all_on_narw_unrelated_filtered_aug5.vcf.gz
	bcftools index ~/scratch/msmc/shapeit/${scaffold}_all_on_narw_filtered_aug8.vcf.gz 
done
```

Phasing VCF for NARW

Per the Statistical population genomics book, and the script provided by Josquin.

	- Removes multi-allelic sites in the VCF
	- Makes a list of sites to be excluded in the main run for phasing by shapeit
	- Runs shapeit with --exclude-snp and -nomcmc generating two output files 
	- These two files can be converted into VCF format by shapeit -convert. 
	- Merge the phased VCF sample1.chr$CHR.onlyPhased.vcf.gz and the unphased (original) VCF sample1.chr$CHR.fixedformat.vcf.gz, keeping all unphased sites from the original VCF, but replacing the phased ones.


```
module load bcftools/1.11 shapeit/2.r904

#!/bin/bash
#SBATCH --job-name=phasing_vcf
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-22
#SBATCH --cpus-per-task=21
#SBATCH --mem=50G
#SBATCH --time=12:00:00


CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scaffold_list)

cp HiC_scaffold_${CHR}_all_on_narw_filtered_aug8.vcf.gz ${SLURM_TMPDIR}
cd ${SLURM_TMPDIR}


#Set up Input Files
UNPHASED_VCF=HiC_scaffold_${CHR}_all_on_narw_filtered.vcf.gz
UNPHASED_VCF_NOMAS=HiC_scaffold_${CHR}.noMultiAllelicSites.vcf.gz
LOG_ALIGN=HiC_scaffold_${CHR}.alignments
LOG_MAIN=HiC_scaffold_${CHR}.main

PHASED_HAPS=HiC_scaffold_${CHR}.phased.haps.gz
PHASED_SAMPLE=HiC_scaffold_${CHR}.phased.samples
PHASED_VCF=HiC_scaffold_${CHR}.onlyPhased.vcf
NEW_HEAD=HiC_scaffold_${CHR}_newhead.onlyPhased.vcf

LOG_CONVERT=HiC_scaffold_${CHR}.convert
FINAL_VCF=HiC_scaffold_${CHR}.phased.vcf


#Preparation
bcftools view -M 2 -O z $UNPHASED_VCF > $UNPHASED_VCF_NOMAS
shapeit -check --input-vcf $UNPHASED_VCF_NOMAS --output-log $LOG_ALIGN --thread 20

#Main run with force arguement
shapeit -V $UNPHASED_VCF_NOMAS --output-max $PHASED_HAPS $PHASED_SAMPLE --output-log $LOG_MAIN --thread 20 --force

shapeit -convert --input-haps $PHASED_HAPS $PHASED_SAMPLE --output-vcf $PHASED_VCF --output-log $LOG_CONVERT --thread 20

echo "SHAPEIT DONE"

#new header on phased vcf
bcftools view -h $UNPHASED_VCF | bcftools reheader -h - $PHASED_VCF -o ${NEW_HEAD}
bcftools view -Oz -o $NEW_HEAD.gz $NEW_HEAD
bcftools index ${NEW_HEAD}.gz
echo "REHEADER DONE"


#Extracting sites in the unphased vcf, not present in the phased sites (i.e. those sites which could not be phased)

bcftools isec -n -1 -p ${CHR}/ -c all $UNPHASED_VCF $NEW_HEAD.gz
echo "ISEC DONE"

#Concatenate unphased sites with phased sites and sort the vcf.
bcftools concat $NEW_HEAD ${CHR}/0000.vcf -Oz -o $FINAL_VCF
echo "CONCAT DONE"
bcftools sort -Oz -o HiC_scaffold_${CHR}_sorted.vcf.gz $FINAL_VCF

mv HiC_scaffold_${CHR}* ${SLURM_SUBMIT_DIR}
```

**Creating masking file for NARW**

```
#!/bin/bash
#SBATCH --job-name=gvcf_bed_narw
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-22
#SBATCH --cpus-per-task=15
#SBATCH --mem=20G
#SBATCH --time=4:00:00

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scaffold_list)
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 22; do \
	for SAMPLE in EGL140-1; do \
			gvcf2bed -I ~/scratch/gvcf/NARW/HiC_scaffold_${CHR}/${SAMPLE}-HiC_scaffold_${CHR}.g.vcf.gz -O ${SAMPLE}_${CHR}.bed;
done; done
```

Can use the masking files I recreated for the first run.

Will also rerun this on SRW before any files are erased and I need to the VCF files.

Next I need to create a bed file that identifies regions of the genome that we masked and/or filtered. I will do this using bedtools subtract to variants not identified in the filtered VCF, but that are listed in the merged vcf. This will output a bed file of variant positions we discarded for a number of reasons. 
```
#!/bin/bash
#SBATCH --job-name=gvcf_bed_narw
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-21
#SBATCH --cpus-per-task=15
#SBATCH --mem=20G
#SBATCH --time=4:00:00

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/scratch/msmc/createBed/scaffold_list)

for SAMPLE in EGL00252-1 EGL013-3qa EGL308-1a EGL140-1 EGL183-1 EGL336_1b EGL254-1 EGL272-1 SID179132; do \
	vcftools --gzvcf ~/scratch/msmc/shapeit/phased/HiC_scaffold_${CHR}_sorted.vcf.gz  --indv 	${SAMPLE} --max-missing 1.0 --recode --out ${SAMPLE}_${CHR}.filtered
	vcftools --gzvcf ~/projects/def-frasiert/RW_WGS/vcf/ALL/merged_all_on_narw.vcf.gz --indv ${SAMPLE} --chr HiC_scaffold_${CHR} --recode --out ${SAMPLE}_${CHR}.originalmerge 
	bedtools subtract -a ${SAMPLE}_${CHR}.originalmerge.recode.vcf -b ${SAMPLE}_${CHR}.filtered.recode.vcf > ${SAMPLE}_${CHR}_filtered_sites.vcf
	rm ${SAMPLE}_${CHR}.originalmerge.recode.vcf;
done
```

Finally, I correct the header on the vcf of filtered sites, and remove these filtered sites from the original bed file. The new bed file represents the
```
#!/bin/bash
#SBATCH --job-name=gvcf_bed_narw
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=3-21
#SBATCH --cpus-per-task=15
#SBATCH --mem=20G
#SBATCH --time=4:00:00

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/scratch/msmc/createBed/scaffold_list)


for SAMPLE in EGL00252-1 EGL013-3qa EGL308-1a EGL183-1 EGL336_1b EGL254-1 EGL272-1 SID179132 EGL140-1; do \
	bcftools view -h ~/scratch/msmc/shapeit/phased/phased_split/${SAMPLE}_${CHR}.recode.vcf > ~/scratch/msmc/createBed/mask/temp-${SAMPLE}-${CHR}
	
	cat temp-${SAMPLE}-${CHR} ~/scratch/msmc/createBed/${SAMPLE}_${CHR}_filtered_sites.vcf > ${SAMPLE}_${CHR}_filtered_sites_wheader.vcf

	rm temp-${SAMPLE}-${CHR}

# combine bedtools to remove sites filtered from original masking file	
	bedtools subtract -a ~/scratch/msmc/createBed/${SAMPLE}_${CHR}.bed -b ~/scratch/msmc/createBed/mask/${SAMPLE}_${CHR}_filtered_sites_wheader.vcf > ~/scratch/msmc/createBed/mask/${SAMPLE}_${CHR}_msmcmask.bed;
done
```

**Phasing SRW**

```
#prepare input vcf files separate for each scaffold
for scaffold in HiC_scaffold_1 HiC_scaffold_2 HiC_scaffold_3 HiC_scaffold_4 HiC_scaffold_5 HiC_scaffold_6 HiC_scaffold_7 HiC_scaffold_9 HiC_scaffold_10 HiC_scaffold_11 HiC_scaffold_12 HiC_scaffold_13 HiC_scaffold_14 HiC_scaffold_15 HiC_scaffold_16 HiC_scaffold_17 HiC_scaffold_18 HiC_scaffold_19 HiC_scaffold_20 HiC_scaffold_21; do
	bcftools view -Oz -r $scaffold -o ~/scratch/msmc/shapeit/${scaffold}_srw_w_bowhead_filtered_aug8.vcf.gz SRW_w_bowhead_unrelated_filtered_aug5.vcf.gz
	bcftools index ~/scratch/msmc/shapeit/${scaffold}_srw_w_bowhead_filtered_aug8.vcf.gz 
done
```

Running SHAPEIT
```
module load bcftools/1.11 nixpkgs/16.09 shapeit/2.r904

#!/bin/bash
#SBATCH --job-name=phasing_vcf
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-22
#SBATCH --cpus-per-task=21
#SBATCH --mem=50G
#SBATCH --time=12:00:00


CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scaffold_list_srw)

cp HiC_scaffold_${CHR}_srw_w_bowhead_filtered_aug8.vcf.gz ${SLURM_TMPDIR}
cd ${SLURM_TMPDIR}


#Set up Input Files
UNPHASED_VCF=HiC_scaffold_${CHR}_srw_w_bowhead_filtered_aug8.vcf.gz
UNPHASED_VCF_NOMAS=HiC_scaffold_${CHR}_srw_w_bowhead.noMultiAllelicSites.vcf.gz
LOG_ALIGN=HiC_scaffold_${CHR}_srw_w_bowhead.alignments
LOG_MAIN=HiC_scaffold_${CHR}_srw_w_bowhead.main

PHASED_HAPS=HiC_scaffold_${CHR}_srw_w_bowhead.phased.haps.gz
PHASED_SAMPLE=HiC_scaffold_${CHR}_srw_w_bowhead.phased.samples
PHASED_VCF=HiC_scaffold_${CHR}_srw_w_bowhead.onlyPhased.vcf
NEW_HEAD=HiC_scaffold_${CHR}_newhead_srw_w_bowhead.onlyPhased.vcf

LOG_CONVERT=HiC_scaffold_${CHR}_srw_w_bowhead.convert
FINAL_VCF=HiC_scaffold_${CHR}_srw_w_bowhead.phased.vcf


#Preparation
bcftools view -M 2 -O z $UNPHASED_VCF > $UNPHASED_VCF_NOMAS
shapeit -check --input-vcf $UNPHASED_VCF_NOMAS --output-log $LOG_ALIGN --thread 20

#Main run with force arguement
shapeit -V $UNPHASED_VCF_NOMAS --output-max $PHASED_HAPS $PHASED_SAMPLE --output-log $LOG_MAIN --thread 20 --force

shapeit -convert --input-haps $PHASED_HAPS $PHASED_SAMPLE --output-vcf $PHASED_VCF --output-log $LOG_CONVERT --thread 20

echo "SHAPEIT DONE"

#new header on phased vcf
bcftools view -h $UNPHASED_VCF | bcftools reheader -h - $PHASED_VCF -o ${NEW_HEAD}
bcftools view -Oz -o $NEW_HEAD.gz $NEW_HEAD
bcftools index ${NEW_HEAD}.gz
echo "REHEADER DONE"


#Extracting sites in the unphased vcf, not present in the phased sites (i.e. those sites which could not be phased)

bcftools isec -n -1 -p ${CHR}/ -c all $UNPHASED_VCF $NEW_HEAD.gz
echo "ISEC DONE"

#Concatenate unphased sites with phased sites and sort the vcf.
bcftools concat $NEW_HEAD ${CHR}/0000.vcf -Oz -o $FINAL_VCF
echo "CONCAT DONE"
bcftools sort -Oz -o HiC_scaffold_${CHR}_srw_w_bowhead_sorted.vcf.gz $FINAL_VCF

mv HiC_scaffold_${CHR}* ${SLURM_SUBMIT_DIR}/output/
```

### August 9, 2022 ###

**Phasing SRW**

Reran the code above using 10 threads as the program doesn't like it when there are too many relative to sample size.

**Gzip masking files**

Running a script that is just gzip on all files.

**MSMC Input**

I started to generate the msmc input files and I realized, I forgot the --max-missing=1 in vcftools. I am reunning that code to split the vcf files and then will recreate the masking files using filter_ned_narw.sh

```
#!/bin/bash
#SBATCH --job-name=filter_bed_narw
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-21
#SBATCH --cpus-per-task=15
#SBATCH --mem=20G
#SBATCH --time=6:00:00

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/scratch/msmc/createBed/scaffold_list)

for SAMPLE in EGL00252-1 EGL013-3qa EGL308-1a EGL312-1a EGL183-1 EGL336_1b EGL254-1 EGL272-1 SID179132 SID181803 EGL140-1 EGL276-1; do \
        vcftools --gzvcf ~/scratch/msmc/shapeit/phased/HiC_scaffold_${CHR}_sorted.vcf.gz --indv ${SAMPLE} --max-missing 1.0 --recode --out ${SAMPLE}_${CHR}.filtered
        vcftools --gzvcf ~/projects/def-frasiert/RW_WGS/vcf/ALL/merged_all_on_narw.vcf.gz --indv ${SAMPLE} --chr HiC_scaffold_${CHR} --recode --out ${SAMPLE}_${CHR}.originalmerge
        bedtools subtract -a ${SAMPLE}_${CHR}.originalmerge.recode.vcf -b ${SAMPLE}_${CHR}.filtered.recode.vcf > ${SAMPLE}_${CHR}_filtered_sites.vcf
        rm ${SAMPLE}_${CHR}.originalmerge.recode.vcf

#correct header
        bcftools view -h ${SAMPLE}_${CHR}.filtered.recode.vcf > ~/scratch/msmc/createBed/mask/temp-${SAMPLE}-${CHR}

        cat ~/scratch/msmc/createBed/mask/temp-${SAMPLE}-${CHR} ~/scratch/msmc/createBed/${SAMPLE}_${CHR}_filtered_sites.vcf > ~/scratch/msmc/createBed/mask/${SAMPLE}_${CHR}_filtered_sites_wheader.vcf

        rm ~/scratch/msmc/createBed/mask/temp-${SAMPLE}-${CHR}

#combine bedtools to remove sites filtered from original masking file
        bedtools subtract -a ~/scratch/msmc/createBed/${SAMPLE}_${CHR}.bed -b ~/scratch/msmc/createBed/mask/${SAMPLE}_${CHR}_filtered_sites_wheader.vcf > ~/scratch/msmc/createBed/mask/${SAMPLE}_${CHR}_msmcmask.bed
        gzip ~/scratch/msmc/createBed/mask/${SAMPLE}_${CHR}_msmcmask.bed;

done
```

**SRW MSMC BED FILES**

SRW w bowhead finished phasing.
```
#!/bin/bash
#SBATCH --job-name=filter_bed_narw
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-20
#SBATCH --cpus-per-task=15
#SBATCH --mem=20G
#SBATCH --time=6:00:00

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/scratch/msmc/createBed/scaffold_list_srw)

for SAMPLE in Eau018 Eau017 Eau019 Eau023 Eau029 Eau034A Eau283 Eau10b Eau7 Eau9c; do \
        vcftools --gzvcf ~/scratch/msmc/shapeit/output/HiC_scaffold_${CHR}_srw_w_bowhead_sorted.vcf.gz  --indv ${SAMPLE} --max-missing 1.0 --recode --out ${SAMPLE}_${CHR}.filtered
        vcftools --gzvcf ~/projects/def-frasiert/RW_WGS/vcf/SRW/merged_SRW_w_bowhead.vcf.gz --indv ${SAMPLE} --chr HiC_scaffold_${CHR} --recode --out ${SAMPLE}_${CHR}.originalmerge
        bedtools subtract -a ${SAMPLE}_${CHR}.originalmerge.recode.vcf -b ${SAMPLE}_${CHR}.filtered.recode.vcf > ${SAMPLE}_${CHR}_filtered_sites.vcf
        rm ${SAMPLE}_${CHR}.originalmerge.recode.vcf

#correct header
        bcftools view -h ${SAMPLE}_${CHR}.filtered.recode.vcf > ~/scratch/msmc/createBed/mask/temp-${SAMPLE}-${CHR}

        cat ~/scratch/msmc/createBed/mask/temp-${SAMPLE}-${CHR} ~/scratch/msmc/createBed/${SAMPLE}_${CHR}_filtered_sites.vcf > ~/scratch/msmc/createBed/mask/${SAMPLE}_${CHR}_filtered_sites_wheader.vcf

        rm ~/scratch/msmc/createBed/mask/temp-${SAMPLE}-${CHR}

#combine bedtools to remove sites filtered from original masking file
        bedtools subtract -a ~/scratch/msmc/createBed/${SAMPLE}_${CHR}.bed -b ~/scratch/msmc/createBed/mask/${SAMPLE}_${CHR}_filtered_sites_wheader.vcf > ~/scratch/msmc/createBed/mask/${SAMPLE}_${CHR}_msmcmask.bed
        gzip ~/scratch/msmc/createBed/mask/${SAMPLE}_${CHR}_msmcmask.bed;

done
```

### August 10, 2022 ###

**Final MSMC prep**

- gzipped all filtered.recode.vcf files.
- Moved all input files (*_msmcmask.bed.gz, *filtered.recode.vcf.gz) to generate_multisephet_input\

Loaded python environment and started create_MSMC_input.sh on narw 3 samples to start.

First run on NARW included: EGL00252-1, EGL013-3qa, EGL183-1

```
#!/bin/bash
#SBATCH --job-name=generate_MSMC_input
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-21
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G
#SBATCH --time=1:00:00

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scaffold_list)

~/scratch/PY_ENV/generate_multihetsep.py \
--mask=EGL00252-1_${CHR}_msmcmask.bed.gz \
--mask=EGL013-3qa_${CHR}_msmcmask.bed.gz \
--mask=EGL183-1_${CHR}_msmcmask.bed.gz \
./EGL00252-1_${CHR}.filtered.recode.vcf.gz \
./EGL013-3qa_${CHR}.filtered.recode.vcf.gz \
./EGL183-1_${CHR}.filtered.recode.vcf.gz > ../input_files/aug10_narw_${CHR}_multihetsep.txt

```

First run on SRW included: Eau017, Eau018, Eau7
```
#!/bin/bash
#SBATCH --job-name=generate_MSMC_input_srw
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-20
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G
#SBATCH --time=1:00:00

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scaffold_list_srw)

~/scratch/PY_ENV/generate_multihetsep.py \
--mask=Eau017_${CHR}_msmcmask.bed.gz \
--mask=Eau018_${CHR}_msmcmask.bed.gz \
--mask=Eau7_${CHR}_msmcmask.bed.gz \
./Eau017_${CHR}.filtered.recode.vcf.gz \
./Eau018_${CHR}.filtered.recode.vcf.gz \
./Eau7_${CHR}.filtered.recode.vcf.gz > ../input_files/aug10_srw_${CHR}_multihetsep.txt
```

MSMC cannot run with indels and while I ran everything with no missing data, there still appeared to be indels screwing things up. Starting with NARW, I am going to rerun the filtering script, added line to the first vcftools command removing indels.

New code:

```
#!/bin/bash
#SBATCH --job-name=filter_bed_narw
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-21
#SBATCH --cpus-per-task=15
#SBATCH --mem=20G
#SBATCH --time=6:00:00

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/scratch/msmc/createBed/scaffold_list)

for SAMPLE in EGL00252-1 EGL013-3qa EGL308-1a EGL312-1a EGL183-1 EGL336_1b EGL254-1 EGL272-1 SID179132 SID181803 EGL140-1 EGL276-1; do \
        vcftools --gzvcf ~/scratch/msmc/shapeit/phased/HiC_scaffold_${CHR}_sorted.vcf.gz  --indv ${SAMPLE} --max-missing 1.0 --remove-indels --recode --out ${SAMPLE}_${CHR}.filtered
        vcftools --gzvcf ~/projects/def-frasiert/RW_WGS/vcf/ALL/merged_all_on_narw.vcf.gz --indv ${SAMPLE} --chr HiC_scaffold_${CHR} --recode --out ${SAMPLE}_${CHR}.originalmerge
        bedtools subtract -a ${SAMPLE}_${CHR}.originalmerge.recode.vcf -b ${SAMPLE}_${CHR}.filtered.recode.vcf > ${SAMPLE}_${CHR}_filtered_sites.vcf
        rm ${SAMPLE}_${CHR}.originalmerge.recode.vcf

#correct header
        bcftools view -h ${SAMPLE}_${CHR}.filtered.recode.vcf > ~/scratch/msmc/createBed/mask/temp-${SAMPLE}-${CHR}

        cat ~/scratch/msmc/createBed/mask/temp-${SAMPLE}-${CHR} ~/scratch/msmc/createBed/${SAMPLE}_${CHR}_filtered_sites.vcf > ~/scratch/msmc/createBed/mask/${SAMPLE}_${CHR}_filtered_sites_wheader.vcf

        rm ~/scratch/msmc/createBed/mask/temp-${SAMPLE}-${CHR}

#combine bedtools to remove sites filtered from original masking file
        bedtools subtract -a ~/scratch/msmc/createBed/${SAMPLE}_${CHR}.bed -b ~/scratch/msmc/createBed/mask/${SAMPLE}_${CHR}_filtered_sites_wheader.vcf > ~/scratch/msmc/createBed/mask/${SAMPLE}_${CHR}_msmcmask.bed
        gzip ~/scratch/msmc/createBed/mask/${SAMPLE}_${CHR}_msmcmask.bed;

done
```

Reran generate_multihetsep.py

Restarting MSMC2.

```
module load nixpkgs/16.09 gcc/7.3.0 msmc2/2.0.0

#!/bin/bash
#SBATCH --job-name=MSMC2_narw
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=20
#SBATCH --mem=124G
#SBATCH --time=6:00:00

cp input_files/aug10_narw* ${SLURM_TMPDIR}
cd ${SLURM_TMPDIR}

msmc2 -t 20 -o results/narw_msmc2/aug10 \
        ./input_files/aug10_narw_1_multihetsep.txt \
        ./input_files/aug10_narw_2_multihetsep.txt \
        ./input_files/aug10_narw_3_multihetsep.txt \
        ./input_files/aug10_narw_4_multihetsep.txt \
        ./input_files/aug10_narw_5_multihetsep.txt \
        ./input_files/aug10_narw_6_multihetsep.txt \
        ./input_files/aug10_narw_7_multihetsep.txt \
        ./input_files/aug10_narw_8_multihetsep.txt \
        ./input_files/aug10_narw_9_multihetsep.txt \
        ./input_files/aug10_narw_10_multihetsep.txt \
        ./input_files/aug10_narw_11_multihetsep.txt \
        ./input_files/aug10_narw_12_multihetsep.txt \
        ./input_files/aug10_narw_13_multihetsep.txt \
        ./input_files/aug10_narw_14_multihetsep.txt \
        ./input_files/aug10_narw_15_multihetsep.txt \
        ./input_files/aug10_narw_16_multihetsep.txt \
        ./input_files/aug10_narw_17_multihetsep.txt \
        ./input_files/aug10_narw_18_multihetsep.txt \
        ./input_files/aug10_narw_19_multihetsep.txt \
        ./input_files/aug10_narw_20_multihetsep.txt \
        ./input_files/aug10_narw_22_multihetsep.txt

cp results/narw_msmc2/aug10* ${SLURM_SUBMIT_DIR}

```

**Setting up R**

Load newer version of R. (3.7 so the bioconductor package still easy to use) install all packages.

 - Loaded r/4.1.0 mugqic/R_Bioconductor/4.1.0_3.13

```
install.packages("adegenet", repos="https://mirror.rcg.sfu.ca/mirror/CRAN/")
install.packages("ape", repos="https://mirror.rcg.sfu.ca/mirror/CRAN/")
install.packages("pegas", repos="https://mirror.rcg.sfu.ca/mirror/CRAN/")
install.packages("mgcv", repos="https://mirror.rcg.sfu.ca/mirror/CRAN/")
install.packages("snow", repos="https://mirror.rcg.sfu.ca/mirror/CRAN/")
install.packages("data.table", repos="https://mirror.rcg.sfu.ca/mirror/CRAN/")
install.packages("vcfR", repos="https://mirror.rcg.sfu.ca/mirror/CRAN/")
install.packages("stepR", repos="https://mirror.rcg.sfu.ca/mirror/CRAN/")
install.packages("genetics", repos="https://mirror.rcg.sfu.ca/mirror/CRAN/")
install.packages("seqinr", repos="https://mirror.rcg.sfu.ca/mirror/CRAN/")


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings", lib="~/R/x86_64-pc-linux-gnu-library/4.1.2/")

install.packages("~/projects/def-frasiert/RW_WGS/programs/LDJump/LDJump_0.3.1.tar.gz", repos=NULL, type="source")
```


Starting a script to run LDJUMP
I will run some of it in an interactive shell to see if it works at first, but then I will write it as a script for submit with sbatch.

```
##TEST IN INTERACTIVE SHELL##
cp  ~/projects/def-frasiert/RW_WGS/vcf/locked/narw_unrelated_filtered_aug5.vcf.gz ~/scratch/LDJump/
gunzip narw_unrelated_filtered_aug5.vcf #unzip vcf file


# The following with be inside the R script that would be called by the bash script.
require(LDJump)
require()

setwd("~/scratch/LDJump")
vcf_file = "narw_unrelated_filtered_aug5.vcf"
ref_seq = "Eubalaena_glacialis_HiC_min1Mb.fasta"

LDJump(vcf_file, alpha = 0.05, quant = 0.35, segLength = 1000, pathLDhat = "~/projects/def-frasiert/RW_WGS/programs/LDhat/LDhat-master/", pathPhi = "~/projects/def-frasiert/RW_WGS/programs/PhiPack/PhiPack/", format = "vcf", refName = ref_seq, start = 0, constant = F, rescale = F, status = T, polyThres = 0, cores = 1, accept = F, demography = F, regMod = "", out = "NARW_LDJump", lengthofseq = NULL, chr = NULL, startofseq = NULL, endofseq = NULL)
```
It started to run without error.

narw_LDJump.R
```
require(LDJump)

setwd("SLURM_TMPDIR")

vcf_file = "narw_unrelated_filtered_aug5.vcf"
ref_seq = "Eubalaena_glacialis_HiC_min1Mb.fasta"

LDJump(vcf_file, alpha = 0.05, quant = 0.35, segLength = 1000, pathLDhat = "~/projects/def-frasiert/RW_WGS/programs/LDhat/LDhat-master/",
pathPhi = "~/projects/def-frasiert/RW_WGS/programs/PhiPack/PhiPack/", format = "vcf", refName = ref_seq, start = 0, constant = F,
rescale = F, status = T, polyThres = 0, cores = 1, accept = F, demography = F, regMod = "", out = "NARW_LDJump", lengthofseq = NULL, chr = NULL, startofseq = NULL, endofseq = NULL)
```

run_LDJump.sh
```
#!/bin/bash
#SBATCH --job-name=LDJump_narw
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=2
#SBATCH --mem=50G
#SBATCH --time=24:00:00

cp narw_unrelated_filtered_aug5.vcf ${SLURM_TMPDIR}
cp Eubalaena_glacialis_HiC_min1Mb.fasta ${SLURM_TMPDIR}

cd ${SLURM_TMPDIR}

# Run R script

Rscript narw_LDJump.R

# Finally, I would save all output files back to the submit directory
cp NARW_LDJump* ${SLURM_SUBMIT_DIR}
```


Many errors trying to run LDJump in vcf
"Error in 1:(segs + 1) : argument of length 0"
"Calls: LDJump -> vcf_statistics" 

Maybe try parsing phased vcf to fasta and make a couple recombination maps and go from there.

### August 11, 2022 ###

**NARW MSMC2**

This finished running and I plotted the results in Python.

**SRW MSMC2**

I need to redo the filtering/bed masking. Zipped recode.vcf and moved these and .bed.gz and then started generate_multihetsep.py

Running MSMC:

```
#!/bin/bash
#SBATCH --job-name=MSMC2_srw
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=20
#SBATCH --mem=124G
#SBATCH --time=6:00:00

cp input_files/aug10_srw* ${SLURM_TMPDIR}
cd ${SLURM_TMPDIR}

msmc2 -t 20 -o srw_msmc2_aug10 \
        aug10_srw_1_multihetsep.txt \
        aug10_srw_2_multihetsep.txt \
        aug10_srw_3_multihetsep.txt \
        aug10_srw_4_multihetsep.txt \
        aug10_srw_5_multihetsep.txt \
        aug10_srw_6_multihetsep.txt \
        aug10_srw_7_multihetsep.txt \
        aug10_srw_9_multihetsep.txt \
        aug10_srw_10_multihetsep.txt \
        aug10_srw_11_multihetsep.txt \
        aug10_srw_12_multihetsep.txt \
        aug10_srw_13_multihetsep.txt \
        aug10_srw_14_multihetsep.txt \
        aug10_srw_15_multihetsep.txt \
        aug10_srw_16_multihetsep.txt \
        aug10_srw_17_multihetsep.txt \
        aug10_srw_18_multihetsep.txt \
        aug10_srw_19_multihetsep.txt \
        aug10_srw_20_multihetsep.txt \
        aug10_srw_21_multihetsep.txt

cp srw_msmc2_aug10* ${SLURM_SUBMIT_DIR}/results/
```

**LDJump**

Running vcfR_to_fasta in R:
```
require(LDJump)
seq="narw_unrelated_filtered_aug5.vcf"
ref="Eubalaena_glacialis_HiC_min1Mb.fasta"

vcfR_to_fasta(seq,refName=ref)
```

### August 12, 2022 ###

**LDJump**

LDJump needs to be run on a single chromosome at a time. I subset just HiC_scaffold_22, and I am trying LDJump on that file alone first.
``` bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' narw_unrelated_filtered_scaff22.recode.vcf > narw_unrelated_filtered_scaff22_uniqueID.vcf ```

LDJUMP
```
require(LDJump)
ref="Eubalaena_glacialis_HiC_min1Mb.fasta"
seq="narw_unrelated_filtered_scaff22_uniqueID.vcf"
vcfR_to_fasta(seq,refName=ref,fa_start=1, fa_end=2852764) #this worked to create VCF

LDJump(seq, alpha = 0.05, quant = 0.35, segLength = 1000, pathLDhat = "~/projects/def-frasiert/RW_WGS/programs/LDhat/LDhat-master/", pathPhi = "~/projects/def-frasiert/RW_WGS/programs/PhiPack/PhiPack/", format = "vcf", refName = ref, start = 0, constant = F, rescale = F, status = T, polyThres = 0, cores = 4, accept = F, demography = F, regMod = "", out = "NARW_LDJump", lengthofseq = NULL, chr = HiC_scaffold_22, startofseq = 1, endofseq = 2852764)

# Error in 1:(segs + 1) : argument of length 0

```
The LDJump was killed so possibly an oom error. So ran it the Rscript inside bash script. Error in 1:(segs + 1) : argument of length 0.

Even though this seems to be a known problem with the vcf, it looks as though lengthofseq is required when working with vcf files. Trying with that equal to the sequence length and if that doesn't work I will use LDJump with a fasta file.
Error:
```
vcftools: /cvmfs/soft.mugqic/lfs/7.6/tools/lib/libc.so.6: version `GLIBC_2.17' not found (required by /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/gcccore/9.3.0/lib64/libstdc++.so.6)
vcftools: /cvmfs/soft.mugqic/lfs/7.6/tools/lib/libc.so.6: version `GLIBC_2.16' not found (required by /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/gcccore/9.3.0/lib64/libstdc++.so.6)
vcftools: /cvmfs/soft.mugqic/lfs/7.6/tools/lib/libc.so.6: version `GLIBC_2.14' not found (required by /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/gcccore/9.3.0/lib64/libstdc++.so.6)
vcftools: /cvmfs/soft.mugqic/lfs/7.6/tools/lib/libc.so.6: version `GLIBC_2.18' not found (required by /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/gcccore/9.3.0/lib64/libstdc++.so.6)
vcftools: /cvmfs/soft.mugqic/lfs/7.6/tools/lib/libc.so.6: version `GLIBC_2.14' not found (required by /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/gcccore/9.3.0/lib64/libgcc_s.so.1)
Error in file(file, "r") : cannot open the connection
Calls: LDJump ... vcf_statistics -> vcfR_to_fasta -> <Anonymous> -> scan -> file
In addition: Warning message:
In file(file, "r") :
  cannot open file 'temp/sel_1_1001.recode.vcf': No such file or directory
Execution halted
cp: cannot stat 'NARW_LDJump*': No such file or directory
```

I had used vcfR_to_fasta and have a fasta, so I will try the LDJump code for a fasta file with that sequence.

```
require(LDJump)
ref="Eubalaena_glacialis_HiC_min1Mb.fasta"
seq="narw_unrelated_filtered_scaff22_uniqueID.vcf.fasta"
vcfR_to_fasta(seq,refName=ref,fa_start=1, fa_end=2852764) #this worked to create VCF

LDJump(seq, alpha = 0.05, quant = 0.35, segLength = 1000, pathLDhat = "~/projects/def-frasiert/RW_WGS/programs/LDhat/LDhat-master/", pathPhi = "~/projects/def-frasiert/RW_WGS/programs/PhiPack/PhiPack/", format = "fasta", refName = ref, start = 0, constant = F, rescale = F, status = T, polyThres = 0, cores = 4, accept = TRUE, demography = F, regMod = "", out = "NARW_LDJump", lengthofseq = NULL)
```

This didn't work as there was no variants in the fasta file. I suspect I would need a multisample fasta file for each scaffold.

I am continuing to work with the vcf version. It required more memory, a different version of vcftools, more time and the option accept=TRUE to estimate recombination when there are no variants in the segment.

*LDJump Scaffold 1*

```
# Separate only a single chromosome
vcftools --vcf narw_unrelated_filtered_aug5.vcf --chr HiC_scaffold_1 --recode --out narw_unfiltered_scaff1

# Correct the ID column
bcftools annotate --set-id +'%CHROM\_%POS' narw_unrelated_filtered_scaff1.recode.vcf > narw_scaff1_uniqueID.vcf 
```


### August 13, 2022 ###

**MSMC2 SRW**

I plotted MSMC for SRW in a Jupyter notebook. I will rerun with a few other individuals to compare. 


**LDJump**

There was a problem with creating the fasta from the vcf. I am going to try anf take the vcf from a single individual phased at a single locus and then separate into two haploid files (one for each sequence) and then run bcftools consensus with the reference fasta to write a fasta sequence for each.


### August 15, 2022 ###

**LDJump**

I am going to try running the vcfR to fasta again.

```
require(LDJump)
seq="EGL00252-1_1.filtered.recode_uniqueID.vcf"
ref="Eubalaena_glacialis_HiC_min1Mb.fasta"
scaff="HiC_scaffold_1"
end=104503784

vcfR_to_fasta(seq, refName=ref, attr_name="HiC_scaffold_1", fa_start=1, fa_end=end, start=1, ref=ref)
```

```
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' EGL00252-1_1.filtered.recode.vcf > EGL00252-1_1.filtered.recode_uniqueID.vcf 
```
I identified the problem in vcfR. In the names=names(ref) called by seqinr, there it lists all scaffolds and then only writes to the last one instead of the identified one. I am going to try an recreate this function on my own here. 

```
vcfR_to_fasta_edit = function(seqName, refName = NULL, ext.ind = T, cons = F, ext.haps = T, start = NULL, ref = NULL, fa_start= NULL, fa_end = NULL, attr_name = NULL) {

  pop = vcfR::read.vcfR(file = seqName, limit = 1e+07, nrows = -1, skip = 0, cols = NULL, convertNA = F, verbose = T)
  
  seqinr::write.fasta(sequences = ref[["HiC_scaffold_1"]][1:104503784], names = attr_name, nbchar = 80, file.out = paste0("temp/", fa_start, "_", fa_end, "_out.fa"))
  
  pop.dnabin = vcfR::vcfR2DNAbin(pop, extract.indels = T, consensus = F,
                                 extract.haps = T, ref.seq = read.dna(paste0("temp/", fa_start, "_", fa_end, "_out.fa"), format = "fasta"),
                                 start.pos = start, verbose = TRUE)

  ape::write.dna(x = pop.dnabin, file = paste0(seqName, ".fasta", sep = ""), format = "fasta", colsep = "")
  print(paste(seqName, " converted to fasta file: ", seqName, ".fasta", sep = ""))
}
```

**SRW MSMC2**

I want to re run MSMC2 on 3 different SRW individuals.

Using: Eau023 Eau029 and Eau283
```
#!/bin/bash
#SBATCH --job-name=generate_MSMC_input_srw
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-20
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G
#SBATCH --time=1:00:00

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scaffold_list_srw)

~/scratch/PY_ENV/generate_multihetsep.py \
--mask=Eau017_${CHR}_msmcmask_aug11.bed.gz \
--mask=Eau018_${CHR}_msmcmask_aug11.bed.gz \
--mask=Eau7_${CHR}_msmcmask_aug11.bed.gz \
./Eau017_${CHR}_aug11.filtered.recode.vcf.gz \
./Eau018_${CHR}_aug11.filtered.recode.vcf.gz \
./Eau7_${CHR}_aug11.filtered.recode.vcf.gz > ../input_files/aug10_srw_${CHR}_multihetsep.txt
```

I am going to try an run vcf2fasta locally and then import the fasta.

### August 16, 2022 ###

**ROH**

Fixed ROH to be 1Mb for both SRW and NARW (SRW before was run with 100kb), and now removed related individuals and sex chromosomes.

Also wrote a script to save the df generated by scikit allel to save the length of each ROH over 10kb.
```
# Import Libraries
import numpy as np
import scipy
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import allel; print('scikit-allel', allel.__version__)
import hmmlearn

narw_vcf = allel.read_vcf('narw_unrelated_filtered_aug5.vcf.gz', fields=['variants/*','calldata/*'])


narw_roh = pd.DataFrame()
for SCAFFOLD in np.unique(narw_vcf['variants/CHROM']):
    narw_gt = allel.GenotypeArray(narw_vcf['calldata/GT'][narw_vcf['variants/CHROM'] == SCAFFOLD])
    for IND in range(0,9):
        gv = narw_gt[:,IND]
        pos = narw_vcf['variants/POS'][narw_vcf['variants/CHROM'] == SCAFFOLD]
        scaffold_size = pos[len(pos)-1]
        roh = allel.roh_mhmm(gv, pos,min_roh=10000, contig_size=scaffold_size)
        roh[0]['SAMPLE']=IND
        roh[0]['CHROMOSOME']=SCAFFOLD
        narw_roh=narw_roh.append(roh[0])

# Save ROH data per sample per scaffold to a .csv file
narw_roh.to_csv('narw_roh_details.csv')
```

### August 17, 2022 ###
 
**LDJump**

I am going to try and create the fasta file from bcftools consensus again. If this works, I could run this on all individuals for each chromosome and concatenate the files to create one fasta file per scaffold. 

Originally this looked too long, but that is because wc -c counts new lines as characters, but if you count the number of lines and subtract that from the total number of characters you get the correct number. There was an issue with non-biallelic sites. So I subset biallelic sites then created the consensus sequences, renamed them, then joined them into a single vcf.

```
bcftools view -m2 -M2 -v snps -Oz -o EGL00252-1_1_test.vcf.gz EGL00252-1_1.filtered.recode.vcf.gz

bcftools index EGL00252-1_1_test.vcf.gz

bcftools consensus -H 1pIu -f Eubalaena_glacialis_HiC_min1Mb_scaffold1.fasta -o EGL00252-1_scaff1_1.fasta EGL00252-1_1_test.vcf.gz
sed -i 's,>HiC_scaffold_1\(.*\),>EGL00252-1_1\1,' EGL00252-1_scaff1_1.fasta

bcftools consensus -H 2pIu -f Eubalaena_glacialis_HiC_min1Mb_scaffold1.fasta -o EGL00252-1_scaff1_2.fasta EGL00252-1_1_test.vcf.gz
sed -i 's,>HiC_scaffold_1\(.*\),>EGL00252-1_2\1,' EGL00252-1_scaff1_2.fasta

cat EGL00252-1_scaff1_1.fasta EGL00252-1_scaff1_2.fasta > EGL00252-1_scaffold1.fasta
```
Now I can run the code above in a script to run over all samples at all scaffolds, so we have a fasta file for each scaffold. I will also need a fasta reference sequence for each scaffold. 

### August 18, 2022 ###

**Stairway Plot**

Changed the generation time for both NARW and SRW to be 32 (mean of pre-whaling estimates from Taylor 2007).

To run stairway plot, I did the following:

Created the batch file:
```
java -cp stairway_plot_es Stairbuilder NARW_unfolded_Aug18.blueprint
java -cp stairway_plot_es Stairbuilder SRW_unfolded_Aug18.blueprint
```

Create a file to submit a scheduled job for the batch:
```
#!/bin/bash
#SBATCH --job-name=stariway_plot_narw
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=6
#SBATCH --mem=25G
#SBATCH --time=12:00:00

bash NARW_unfolded_Aug18.blueprint.sh
```

Ran this job ```sbatch run_stairwayplot.sh```.


### August 19, 2022 ###

**Create input for LDJump**

Testing scaffold22 on one sample to see if LDJump will work and just needs more time/resources.

```
bcftools view -m2 -M2 -v snps -Oz -o EGL00252-1_22_biallelicsnps.vcf.gz ~/scratch/msmc/generate_multisephet_input/EGL00252-1_22.filtered.recode.vcf.gz

bcftools index EGL00252-1_22_biallelicsnps.vcf.gz

bcftools consensus -H 1pIu -f NARW_ref_scaff22.fasta -o EGL00252-1_scaff22_1.fasta EGL00252-1_22_biallelicsnps.vcf.gz
sed -i 's,>HiC_scaffold_22\(.*\),>EGL00252-1_22_1\1,' EGL00252-1_scaff22_1.fasta

bcftools consensus -H 2pIu -f NARW_ref_scaff22.fasta -o EGL00252-1_scaff22_2.fasta EGL00252-1_22_biallelicsnps.vcf.gz
sed -i 's,>HiC_scaffold_22\(.*\),>EGL00252-1_22_2\1,' EGL00252-1_scaff22_2.fasta

bcftools view -m2 -M2 -v snps -Oz -o EGL254-1_22_biallelicsnps.vcf.gz ~/scratch/msmc/generate_multisephet_input/EGL254-1_22.filtered.recode.vcf.gz

bcftools index EGL254-1_22_biallelicsnps.vcf.gz

bcftools consensus -H 1pIu -f NARW_ref_scaff22.fasta -o EGL254-1_scaff22_1.fasta EGL254-1_22_biallelicsnps.vcf.gz
sed -i 's,>HiC_scaffold_22\(.*\),>EGL254-1_22_1\1,' EGL254-1_scaff22_1.fasta

bcftools consensus -H 2pIu -f NARW_ref_scaff22.fasta -o EGL254-1_scaff22_2.fasta EGL254-1_22_biallelicsnps.vcf.gz
sed -i 's,>HiC_scaffold_22\(.*\),>EGL254-1_22_2\1,' EGL254-1_scaff22_2.fasta


bcftools view -m2 -M2 -v snps -Oz -o EGL272-1_22_biallelicsnps.vcf.gz ~/scratch/msmc/generate_multisephet_input/EGL272-1_22.filtered.recode.vcf.gz

bcftools index EGL272-1_22_biallelicsnps.vcf.gz

bcftools consensus -H 1pIu -f NARW_ref_scaff22.fasta -o EGL272-1_scaff22_1.fasta EGL272-1_22_biallelicsnps.vcf.gz
sed -i 's,>HiC_scaffold_22\(.*\),>EGL272-1_22_1\1,' EGL272-1_scaff22_1.fasta

bcftools consensus -H 2pIu -f NARW_ref_scaff22.fasta -o EGL272-1_scaff22_2.fasta EGL272-1_22_biallelicsnps.vcf.gz
sed -i 's,>HiC_scaffold_22\(.*\),>EGL272-1_22_2\1,' EGL272-1_scaff22_2.fasta


cat EGL272-1_scaff22_1.fasta EGL272-1_scaff22_2.fasta EGL00252-1_scaff22_1.fasta EGL00252-1_scaff22_2.fasta EGL254-1_scaff22_1.fasta EGL254-1_scaff22_2.fasta> EGL00252-1_scaffold22.fasta
```

This seems to work. Now I will create a short script that will iterate over all SRW sample and all SRW chromosomes to generate the fasta files.

```
for CHR in 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 17 18 19 20 21; do
	for IND in Eau017 Eau018 Eau7 Eau9c Eau023 Eau029 Eau034A Eau10b Eau283; do
		bcftools view -m2 -M2 -v snps -Oz -o ${IND}_${CHR}_biallelicsnps.vcf.gz ~/scratch/msmc/generate_multisephet_input/${IND}_${CHR}.filtered.recode.vcf.gz

		bcftools index ${IND}_${CHR}_biallelicsnps.vcf.gz

		bcftools consensus -H 1pIu -f SRW_ref_scaff${CHR}.fasta -o ${IND}_scaff${CHR}_1.fasta ${IND}_${CHR}_biallelicsnps.vcf.gz
		sed -i "s,>HiC_scaffold_${CHR},>${IND}_${CHR}_1," ${IND}_scaff${CHR}_1.fasta

		bcftools consensus -H 2pIu -f SRW_ref_scaff${CHR}.fasta -o ${IND}_scaff${CHR}_2.fasta ${IND}_${CHR}_biallelicsnps.vcf.gz
		sed -i "s,>HiC_scaffold_${CHR},>${IND}_${CHR}_2," ${IND}_scaff${CHR}_2.fasta
		
		cat ${IND}_scaff${CHR}_1.fasta ${IND}_scaff${CHR}_2.fasta > ${IND}_scaff${CHR}_haplotypes.fasta
		
		rm ${IND}_${CHR}_biallelicsnps.vcf.gz
		rm ${IND}_${CHR}_biallelicsnps.vcf.gz.csi
		rm ${IND}_scaff${CHR}_1.fasta
		rm ${IND}_scaff${CHR}_2.fasta
		
	done;

	cat Eau017_scaff${CHR}_haplotypes.fasta Eau018_scaff${CHR}_haplotypes.fasta Eau7_scaff${CHR}_haplotypes.fasta Eau9c_scaff${CHR}_haplotypes.fasta Eau023_scaff${CHR}_haplotypes.fasta Eau029_scaff${CHR}_haplotypes.fasta Eau034A_scaff${CHR}_haplotypes.fasta Eau10b_scaff${CHR}_haplotypes.fasta Eau283_scaff${CHR}_haplotypes.fasta > SRW_scaffold${CHR}.fasta
done	
```

### August 22, 2022 ###

**IBDSeq**

There are two components of IBDNe (IBDSeq and IBDNe). IBDNe requires the recombination rate. This is what I am trying to figure out in LDJump. 

In IBDSeq, we calculate the identity by descent blocks with different r2max values. The default is 0.15, but that yielded very low numbers. I may try r2max=0.3. First, I need to create vcf files with no missing data.

```
# Ran in the IBDNe/IBDSeq/

vcftools --gzvcf ~/projects/def-frasiert/RW_WGS/vcf/locked/narw_unrelated_filtered_aug5.vcf.gz --max-missing 1.0 --recode --stdout | gzip -c > narw_all_filters_no_missing_aug22.vcf.gz

vcftools --gzvcf ~/projects/def-frasiert/RW_WGS/vcf/locked/srw_unrelated_filtered_aug5.vcf.gz --max-missing 1.0 --recode --stdout | gzip -c > srw_all_filters_no_missing_aug22.vcf.gz
```

I am going to set it up to run on each scaffold with three different r2max values: 0.2, 0.3, 0.4

```
for R in 0.2 0.3 0.4; do \
	for SCAFFOLD in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 22; do \
		java -Xmx8G -jar ~/projects/def-frasiert/RW_WGS/programs/IBDseq/ibdseq.r1206.jar \
			gt=narw_all_filters_no_missing_aug22.vcf.gz out=r${R}/narw_ibdseq_scaff${SCAFFOLD}_${R} nthreads=8 chrom=HiC_scaffold_${SCAFFOLD} r2max=${R};
done; done

for R in 0.2 0.3 0.4; do \
	for SCAFFOLD in 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 17 18 19 20 21; do \
		java -Xmx8G -jar ~/projects/def-frasiert/RW_WGS/programs/IBDseq/ibdseq.r1206.jar \
			gt=srw_all_filters_no_missing_aug22.vcf.gz out=r${R}/srw_ibdseq_scaff${SCAFFOLD}_${R} nthreads=8 chrom=HiC_scaffold_${SCAFFOLD} r2max=${R};
done; done
```

### August 23, 2022 ###

** PLINK MAPS **

Creating PLINK Map files. I am going to use the recombination rate we created as constant for each scaffold.

```
#generate basic plink map file

module load plink/1.9b_6.21-x86_64

plink --vcf ~/scratch/IBDNe/IBDSeq/narw_all_filters_no_missing_aug22.vcf.gz --recode --allow-extra-chr --out narw_ibdne

plink --vcf ~/scratch/IBDNe/IBDSeq/srw_all_filters_no_missing_aug22.vcf.gz --recode --allow-extra-chr --out srw_ibdne

sed -i 's,\(HiC_scaffold_.*\)\t\.\(\t0\t\)\(.*\),\1\t\1_\3\2\3,' narw_ibdne.map
```

### August 24, 2022 ###

**LDJump**

LDJump didn't really work. I let it run on 1 scaffold on SRW for 24 hours. It quit with an out of memory error and another out of bounds error. 

Based on the findings yesterday that the genetic map should have relative positions and not absolute genetic map distance (ie. the genetic morgans should be recorded as the distance in centimorgans to the previous base).

Need to adjust the map files based on recombination rate set differently for each scaffold.
```
#Adjust the genetic distance



for (s in c("HiC_scaffold_1", "HiC_scaffold_2", "HiC_scaffold_3", "HiC_scaffold_4", "HiC_scaffold_5", "HiC_scaffold_6", "HiC_scaffold_7","HiC_scaffold_9", "HiC_scaffold_10", "HiC_scaffold_11", "HiC_scaffold_12", "HiC_scaffold_13", "HiC_scaffold_14", "HiC_scaffold_15", "HiC_scaffold_16", "HiC_scaffold_17", "HiC_scaffold_18", "HiC_scaffold_19", "HiC_scaffold_20", "HiC_scaffold_21")){
sites <- which(srw_ibdne_test[,1]==s)
sites <- sites[-1]
for (i in sites){
srw_ibdne_test[i,3] = (srw_ibdne_test[i,4] - srw_ibdne_test[i-1,4]) * 1.1 / 1000000
}}


#Will need to change the last and first line of each scaffold manually.
```

This took too long to run on the interactive node. Saved it to a job called fix_genetic_distance_srw.R and ran as a scheduled job.

Still took too long. Wrote a python script and it ran in a minute or less. Note to self -> use the right tool for the job, even if it is less intuitive or harder to write. 

Python script:
```
import pandas as pd
import numpy as np

df = pd.read_table(r'srw_ibdne_test.map', sep='\s+|\t+',  header=None)

arr = np.array(df)

scaffolds = ['HiC_scaffold_1','HiC_scaffold_2','HiC_scaffold_3','HiC_scaffold_4','HiC_scaffold_5','HiC_scaffold_6','HiC_scaffold_7','HiC_scaffold_9','HiC_scaffold_10','HiC_scaffold_11','HiC_scaffold_12','HiC_scaffold_13','HiC_scaffold_14','HiC_scaffold_15','HiC_scaffold_16','HiC_scaffold_17','HiC_scaffold_18','HiC_scaffold_19','HiC_scaffold_20','HiC_scaffold_21']

for s in scaffolds:
    print(s)
    rows=np.where(arr[:,0]==s)
    sites=np.delete(rows,0)
    for pos in sites:
        arr[pos,2] = (arr[pos,3] - arr[pos-1,3]) * 1.1 / 1000000

dataframe = pd.DataFrame(arr)
pd.DataFrame(dataframe).to_csv('srw_ibdne_gendist.map', sep='\t', header=False, index=False)
```
Second script to fix first snp:
```
import pandas as pd
import numpy as np

df = pd.read_table(r'srw_ibdne_gendist.map', sep='\s+|\t+',  header=None)

arr = np.array(df)

scaffolds = ['HiC_scaffold_1','HiC_scaffold_2','HiC_scaffold_3','HiC_scaffold_4','HiC_scaffold_5','HiC_scaffold_6','HiC_scaffold_7','HiC_scaffold_9','HiC_scaffold_10','HiC_scaffold_11','HiC_scaffold_12','HiC_scaffold_13','HiC_scaffold_14','HiC_scaffold_15','HiC_scaffold_16','HiC_scaffold_17','HiC_scaffold_18','HiC_scaffold_19','HiC_scaffold_20','HiC_scaffold_21']

for s in scaffolds:
    print(s)
    rows=np.where(arr[:,0]==s)
    pos=rows[0][0]
    arr[pos,2] = (arr[pos,3]) * 1.1 / 1000000

dataframe = pd.DataFrame(arr)
pd.DataFrame(dataframe).to_csv('srw_ibdne_gendist_1.1.map', sep='\t', header=False, index=False)
```

Started running IBDNe.

### August 25, 2022 ###

**IBDNe**

This didn't run properly yesterday because the distance was not in order. I checked the hapmap files and where I thought the position should be relative to the previous position, this is not the case and this matches the error I received. 

Recreated the map files with the following.
```
import pandas as pd
import numpy as np

df = pd.read_table(r'srw_ibdne_test.map', sep='\s+|\t+',  header=None)

arr = np.array(df)

scaffolds = ['HiC_scaffold_1','HiC_scaffold_2','HiC_scaffold_3','HiC_scaffold_4','HiC_scaffold_5','HiC_scaffold_6','HiC_scaffold_7','HiC_scaffold_9','HiC_scaffold_10','HiC_scaffold_11','HiC_scaffold_12','HiC_scaffold_13','HiC_scaffold_14','HiC_scaffold_15','HiC_scaffold_16','HiC_scaffold_17','HiC_scaffold_18','HiC_scaffold_19','HiC_scaffold_20','HiC_scaffold_21']

for s in scaffolds:
    print(s)
    rows=np.where(arr[:,0]==s)
    for pos in rows:
        arr[pos,2] = (arr[pos,3]) * 1.1 / 1000000

dataframe = pd.DataFrame(arr)
pd.DataFrame(dataframe).to_csv('srw_ibdne_gendist.map', sep='\t', header=False, index=False)
```

**MSMC2**

Setting up different iterations to run 10 simulations of 6 haplotypes each.

Creating MSMC input files first:

SRW

1. Eau017, Eau018, Eau7 (Aug10)
2. Eau023, Eau029, Eau283 (Aug15)
3. Eau9c, Eau10b, Eau034A
4. Eau017, Eau023, Eau9c
5. Eau018, Eau029, Eau10b
6. Eau7, Eau283, Eau034A
7. Eau017, Eau018, Eau029
8. Eau7, Eau023, Eau9c
9. Eau283, Eau10b, Eau034A
10. Eau018, Eau023, Eau10b



### August 26, 2022 ###

**MSMC2**

Starting running msmc2 as an array job on SRW sets 3-10 (listed above).

**IBDSeq**

The number of IBD segments was extremely low and that was in large part due to so many variants being filtered out due to the r2max filter. This makes sense due to our LDDecay plots. I am going to rerun IBDSeq on NARW with r2max of 0.8. 

```
for SCAFFOLD in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 22; do \
	java -Xmx7G -jar ~/projects/def-frasiert/RW_WGS/programs/IBDseq/ibdseq.r1206.jar \
		gt=narw_all_filters_no_missing_aug22.vcf.gz out=r0.8/narw_ibdseq_scaff${SCAFFOLD}_0.8 nthreads=2 chrom=HiC_scaffold_${SCAFFOLD} r2max=0.8;
done
```

**IBDNe**

```
#!/bin/bash
#SBATCH --job-name=IBDNE_NARW
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=20
#SBATCH --mem=20G
#SBATCH --time=1:00:00

cat ~/scratch/IBDNe/IBDSeq/r0.8/narw*.ibd | java -jar ~/projects/def-frasiert/RW_WGS/programs/IBDNe/ibdne.23Apr20.ae9.jar map=~/scratch/IBDNe/narw_ibdne_gendist.map out=narw_r2max0.8_recom1.1 nthreads=20 mincm=3
```

### August 29, 2022 ###

HapIBD didn't run in 24 hours.
LDJump on scaff4 of SRW yielded and OOM error.

For IBDNe, I am going to use r2max of 0.8 (to include more variants, to keep sample size high for number of IBD blocks, performance shouldn't deteriorate as evidenced by diCal paper). Then running IBDNe on 5 recombination rates 0.8-1.2 to present variation in likely range as reported by other mammal species. Present 1.0cM/Mbp in paper.

**MSMC2**

In order to finish MSMC2, I need to run 9 more iterations for NARW.

Similar to what I did above on SRW, I will generate input and run MSMC 9 more times for new combinations of samples.

NARW

1. EGL00252-1, EGL013-3qa, EGL183-1 (Aug10)
2. EGL254-1, EGL272-1, EGL140-1
3. EGL336_1b, EGL308-1a, SID179132
4. EGL00252-1, EGL254-1, EGL336_1b
5. EGL013-3qa, EGL272-1, EGL308-1a
6. EGL183-1, EGL140-1, SID179132
7. EGL00252-1, EGL013-1qa, EGL272-1
8. EGL183-1 EGL254-1, EGL336_1b
9. EGL140-1, EGL308-1a, SID179132
10. EGL00252-1, EGL254-1, EGL308-1a


Starting running msmc2 as an array job on NARW sets 2-10 (listed above).

** IBDSeq **

Running IBDSeq with r2max=0.8 for each species.

```
for SCAFFOLD in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 22; do \
	java -Xmx7G -jar ~/projects/def-frasiert/RW_WGS/programs/IBDseq/ibdseq.r1206.jar \
		gt=narw_all_filters_no_missing_aug22.vcf.gz out=r0.8/narw_ibdseq_scaff${SCAFFOLD}_0.8 nthreads=2 chrom=HiC_scaffold_${SCAFFOLD} r2max=0.8;
done
```

```
for SCAFFOLD in 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 17 18 19 20 21; do \
	java -Xmx7G -jar ~/projects/def-frasiert/RW_WGS/programs/IBDseq/ibdseq.r1206.jar \
		gt=srw_all_filters_no_missing_aug22.vcf.gz out=r0.8/srw_ibdseq_scaff${SCAFFOLD}_0.8 nthreads=2 chrom=HiC_scaffold_${SCAFFOLD} r2max=0.8;
done
```

**PLINK MAPS**

First I need to recreate the plink map files.
```
#generate basic plink map file

module load plink/1.9b_6.21-x86_64

plink --vcf ~/scratch/IBDNe/IBDSeq/narw_all_filters_no_missing_aug22.vcf.gz --recode --allow-extra-chr --out narw_ibdne_aug29

plink --vcf ~/scratch/IBDNe/IBDSeq/srw_all_filters_no_missing_aug22.vcf.gz --recode --allow-extra-chr --out srw_ibdne_aug29

sed -i 's,\(HiC_scaffold_.*\)\t\.\(\t0\t\)\(.*\),\1\t\1_\3\2\3,' narw_ibdne_aug29.map
sed -i 's,\(HiC_scaffold_.*\)\t\.\(\t0\t\)\(.*\),\1\t\1_\3\2\3,' srw_ibdne_aug29.map
```

Then I need to correct the genetic distance column using code similar to below, except changing the recombination rate each run to be one of 0.8, 0.9, 1.0, 1.1, 1.2

```
import pandas as pd
import numpy as np

df = pd.read_table(r'srw_ibdne_aug29.map', sep='\s+|\t+',  header=None)

arr = np.array(df)

scaffolds = ['HiC_scaffold_1','HiC_scaffold_2','HiC_scaffold_3','HiC_scaffold_4','HiC_scaffold_5','HiC_scaffold_6','HiC_scaffold_7','HiC_scaffold_9','HiC_scaffold_10','HiC_scaffold_11','HiC_scaffold_12','HiC_scaffold_13','HiC_scaffold_14','HiC_scaffold_15','HiC_scaffold_16','HiC_scaffold_17','HiC_scaffold_18','HiC_scaffold_19','HiC_scaffold_20','HiC_scaffold_21']

for s in scaffolds:
    print(s)
    rows=np.where(arr[:,0]==s)
    for pos in rows:
        arr[pos,2] = (arr[pos,3]) * 1.1 / 1000000

dataframe = pd.DataFrame(arr)
pd.DataFrame(dataframe).to_csv('srw_ibdne_gendist_1.1_aug29.map', sep='\t', header=False, index=False)
```


**IBDNe**
Running IBDNe 5 times for each species (5 different recombination rates and therefore different map files).
```
#!/bin/bash
#SBATCH --job-name=IBDNE_NARW
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=20
#SBATCH --mem=20G
#SBATCH --time=5:00:00

cat ~/scratch/IBDNe/IBDSeq/r0.8/narw*.ibd | java -jar ~/projects/def-frasiert/RW_WGS/programs/IBDNe/ibdne.23Apr20.ae9.jar map=~/scratch/IBDNe/narw_ibdne_gendist_0.8_aug29.map out=narw_r2max0.8_recom0.8 nthreads=20 mincm=4

cat ~/scratch/IBDNe/IBDSeq/r0.8/narw*.ibd | java -jar ~/projects/def-frasiert/RW_WGS/programs/IBDNe/ibdne.23Apr20.ae9.jar map=~/scratch/IBDNe/narw_ibdne_gendist_0.9_aug29.map out=narw_r2max0.8_recom0.9 nthreads=20 mincm=4

cat ~/scratch/IBDNe/IBDSeq/r0.8/narw*.ibd | java -jar ~/projects/def-frasiert/RW_WGS/programs/IBDNe/ibdne.23Apr20.ae9.jar map=~/scratch/IBDNe/narw_ibdne_gendist_1.0_aug29.map out=narw_r2max0.8_recom1.0 nthreads=20 mincm=4

cat ~/scratch/IBDNe/IBDSeq/r0.8/narw*.ibd | java -jar ~/projects/def-frasiert/RW_WGS/programs/IBDNe/ibdne.23Apr20.ae9.jar map=~/scratch/IBDNe/narw_ibdne_gendist_1.1_aug29.map out=narw_r2max0.8_recom1.1 nthreads=20 mincm=4

cat ~/scratch/IBDNe/IBDSeq/r0.8/narw*.ibd | java -jar ~/projects/def-frasiert/RW_WGS/programs/IBDNe/ibdne.23Apr20.ae9.jar map=~/scratch/IBDNe/narw_ibdne_gendist_1.2_aug29.map out=narw_r2max0.8_recom1.2 nthreads=20 mincm=4
```
And for SRW:
```
#!/bin/bash
#SBATCH --job-name=IBDNE_SRW
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=20
#SBATCH --mem=20G
#SBATCH --time=5:00:00

cat ~/scratch/IBDNe/IBDSeq/r0.8/srw*.ibd | java -jar ~/projects/def-frasiert/RW_WGS/programs/IBDNe/ibdne.23Apr20.ae9.jar map=~/scratch/IBDNe/srw_ibdne_gendist_0.8_aug29.map out=srw_r2max0.8_recom0.8 nthreads=20 mincm=4

cat ~/scratch/IBDNe/IBDSeq/r0.8/srw*.ibd | java -jar ~/projects/def-frasiert/RW_WGS/programs/IBDNe/ibdne.23Apr20.ae9.jar map=~/scratch/IBDNe/srw_ibdne_gendist_0.9_aug29.map out=srw_r2max0.8_recom0.9 nthreads=20 mincm=4

cat ~/scratch/IBDNe/IBDSeq/r0.8/srw*.ibd | java -jar ~/projects/def-frasiert/RW_WGS/programs/IBDNe/ibdne.23Apr20.ae9.jar map=~/scratch/IBDNe/srw_ibdne_gendist_1.0_aug29.map out=srw_r2max0.8_recom1.0 nthreads=20 mincm=4

cat ~/scratch/IBDNe/IBDSeq/r0.8/srw*.ibd | java -jar ~/projects/def-frasiert/RW_WGS/programs/IBDNe/ibdne.23Apr20.ae9.jar map=~/scratch/IBDNe/srw_ibdne_gendist_1.1_aug29.map out=srw_r2max0.8_recom1.1 nthreads=20 mincm=4

cat ~/scratch/IBDNe/IBDSeq/r0.8/srw*.ibd | java -jar ~/projects/def-frasiert/RW_WGS/programs/IBDNe/ibdne.23Apr20.ae9.jar map=~/scratch/IBDNe/srw_ibdne_gendist_1.2_aug29.map out=srw_r2max0.8_recom1.2 nthreads=20 mincm=4
```

**IBDSeq**
Also ran the below for r2max=0.3
```
for SCAFFOLD in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 22; do \
	java -Xmx7G -jar ~/projects/def-frasiert/RW_WGS/programs/IBDseq/ibdseq.r1206.jar \
		gt=narw_all_filters_no_missing_aug22.vcf.gz out=r0.5/narw_ibdseq_scaff${SCAFFOLD}_0.5 nthreads=2 chrom=HiC_scaffold_${SCAFFOLD} r2max=0.5;
done
```

```
for SCAFFOLD in 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 17 18 19 20 21; do \
	java -Xmx7G -jar ~/projects/def-frasiert/RW_WGS/programs/IBDseq/ibdseq.r1206.jar \
		gt=srw_all_filters_no_missing_aug22.vcf.gz out=r0.5/srw_ibdseq_scaff${SCAFFOLD}_0.5 nthreads=2 chrom=HiC_scaffold_${SCAFFOLD} r2max=0.5;
done
```

IBDNe results seem to be robust to mincm (at least in lower generations), and the recombination rate doesn't seem to change much when the plots are overlaid at least for the past 25-50 generations

I am going to test the 3 parameters for r2max and see how they work. 

Overall, results pretty robust to parameter selection. Interested to see how r2max works out.

### August 30, 2022 ###

**MSMC2**

Two of the simulations did not run to completion for unknown reasons.

I am rerunning just sets 5 and 10 for NARW and will include those in the plots. 

**MSMC-IM**

Inorder to run MSMC-IM, we need to run MSMC on within population and between population pairs. 

I am going to create the bed files and filtered variant files for 4 individuals (2 NARW and 2 SRW) from the all_on_narw variant files.

*Phasing with ShapeIT2*

The all_on_narw vcf from Aug8 was already split into scaffolds and indexed. 
Each multi sample scaffold was already phased and found in the ~/scratch/msmc/shapeit/phased/ directory.


```
#!/bin/bash
#SBATCH --job-name=gvcf_bed_narw
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-22
#SBATCH --cpus-per-task=15
#SBATCH --mem=20G
#SBATCH --time=4:00:00

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scaffold_list)

for SAMPLE in EGL00252-1 EGL336_1b; do \
        gvcf2bed -I ~/scratch/gvcf/NARW/HiC_scaffold_${CHR}/${SAMPLE}-HiC_scaffold_${CHR}.g.vcf.gz -O ${SAMPLE}_${CHR}.bed;
done

for SAMPLE in Eau017 Eau283; do \
        gvcf2bed -I ~/scratch/gvcf/SRWonNARW/HiC_scaffold_${CHR}/${SAMPLE}-HiC_scaffold_${CHR}_SRWonNARW.g.vcf.gz -O ${SAMPLE}_${CHR}.bed;
done
```


```
#!/bin/bash
#SBATCH --job-name=filter_bed_all
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-21
#SBATCH --cpus-per-task=15
#SBATCH --mem=20G
#SBATCH --time=6:00:00

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/scratch/msmc/createBed/scaffold_list)

for SAMPLE in EGL00252-1 EGL336_1b Eau017 Eau283; do \
        vcftools --gzvcf ~/scratch/msmc/shapeit/phased/HiC_scaffold_${CHR}_sorted.vcf.gz  --indv ${SAMPLE} --max-missing 1.0 --remove-indels --recode --out ${SAMPLE}_${CHR}_all.filtered
        vcftools --gzvcf ~/projects/def-frasiert/RW_WGS/vcf/ALL/merged_all_on_narw.vcf.gz --indv ${SAMPLE} --chr HiC_scaffold_${CHR} --recode --out ${SAMPLE}_${CHR}_all.originalmerge
        bedtools subtract -a ${SAMPLE}_${CHR}_all.originalmerge.recode.vcf -b ${SAMPLE}_${CHR}_all.filtered.recode.vcf > ${SAMPLE}_${CHR}_all_filtered_sites.vcf
        rm ${SAMPLE}_${CHR}_all.originalmerge.recode.vcf

#correct header
        bcftools view -h ${SAMPLE}_${CHR}_all.filtered.recode.vcf > ~/scratch/msmc/createBed/mask/temp-${SAMPLE}-${CHR}_all

        cat ~/scratch/msmc/createBed/mask/temp-${SAMPLE}-${CHR}_all ~/scratch/msmc/createBed/${SAMPLE}_${CHR}_all_filtered_sites.vcf > ~/scratch/msmc/createBed/mask/${SAMPLE}_${CHR}_all_filtered_sites_wheader.vcf

        rm ~/scratch/msmc/createBed/mask/temp-${SAMPLE}-${CHR}_all

#combine bedtools to remove sites filtered from original masking file
        bedtools subtract -a ~/scratch/msmc/createBed/all_on_narw/${SAMPLE}_${CHR}.bed -b ~/scratch/msmc/createBed/mask/${SAMPLE}_${CHR}_all_filtered_sites_wheader.vcf > ~/scratch/msmc/createBed/mask/${SAMPLE}_${CHR}_all_msmcmask.bed
        gzip ~/scratch/msmc/createBed/mask/${SAMPLE}_${CHR}_all_msmcmask.bed;

done
```

Ran the following to generate MSMC inputs:
```
#!/bin/bash
#SBATCH --job-name=generate_MSMC_input
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-21
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G
#SBATCH --time=1:00:00

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scaffold_list)

~/scratch/PY_ENV/generate_multihetsep.py \
--mask=EGL00252-1_${CHR}_all_msmcmask.bed.gz \
--mask=EGL336_1b_${CHR}_all_msmcmask.bed.gz \
--mask=Eau017_${CHR}_all_msmcmask.bed.gz \
--mask=Eau283_${CHR}_all_msmcmask.bed.gz \
./EGL00252-1_${CHR}_all.filtered.recode.vcf.gz \
./EGL336_1b_${CHR}_all.filtered.recode.vcf.gz \
./Eau017_${CHR}_all.filtered.recode.vcf.gz \
./Eau283_${CHR}_all.filtered.recode.vcf.gz > ../input_files/all_${CHR}_multihetsep.txt
```

Started Running MSMC in 3 ways:

On the 4 NARW haplotypes:
```
#!/bin/bash
#SBATCH --job-name=MSMC2_narw
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=20
#SBATCH --mem=124G
#SBATCH --time=24:00:00

cp input_files/all* ${SLURM_TMPDIR}
cd ${SLURM_TMPDIR}

msmc2 -t 20 -I 0,1,2,3 -o ALL_NARW_msmc2 \
        all_1_multihetsep.txt \
        all_2_multihetsep.txt \
        all_3_multihetsep.txt \
        all_4_multihetsep.txt \
        all_5_multihetsep.txt \
        all_6_multihetsep.txt \
        all_7_multihetsep.txt \
        all_8_multihetsep.txt \
        all_9_multihetsep.txt \
        all_10_multihetsep.txt \
        all_11_multihetsep.txt \
        all_12_multihetsep.txt \
        all_13_multihetsep.txt \
        all_14_multihetsep.txt \
        all_15_multihetsep.txt \
        all_16_multihetsep.txt \
        all_17_multihetsep.txt \
        all_18_multihetsep.txt \
        all_19_multihetsep.txt \
        all_20_multihetsep.txt \
        all_22_multihetsep.txt

cp ALL_NARW_msmc2* ${SLURM_SUBMIT_DIR}/results/
```

On the 4 SRW haplotypes:

Same as above except -I 4,5,6,7

On the between species haplotypes:

Same as above except -I 0-4,0-5,0-6,0-7,1-4,1-5,1-6,1-7,2-4,2-5,2-6,2-7,3-4,3-5,3-6,3-7



**IBDNE**

I am also going to rerun IBDNE with different r2max parameters: 0.3 and 0.5.

Started running IBDNe on r2max=0.3
Started running IBDNe on r2max=0.5
Also ran IBDSeq and IBDNe for r2max=1.0
Also ran the below for r2max=0.3
```
for SCAFFOLD in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 22; do \
	java -Xmx7G -jar ~/projects/def-frasiert/RW_WGS/programs/IBDseq/ibdseq.r1206.jar \
		gt=narw_all_filters_no_missing_aug22.vcf.gz out=r1.0/narw_ibdseq_scaff${SCAFFOLD}_1.0 nthreads=2 chrom=HiC_scaffold_${SCAFFOLD} r2max=1.0;
done
```

```
for SCAFFOLD in 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 17 18 19 20 21; do \
	java -Xmx7G -jar ~/projects/def-frasiert/RW_WGS/programs/IBDseq/ibdseq.r1206.jar \
		gt=srw_all_filters_no_missing_aug22.vcf.gz out=r1.0/srw_ibdseq_scaff${SCAFFOLD}_1.0 nthreads=2 chrom=HiC_scaffold_${SCAFFOLD} r2max=1.0;
done
```

### August 31, 2022 ###

**MSMC / MSMC-IM**

The joint run of MSMC with 2 NARW samples and 2 SRW samples ran successfully on each species independently. There is a known issue with version msmc2/2.0.0 for running between sample comparisions. This was corrected in future updates. I am requesting a newer version be added to the cluster. 

Until this is resolved, I am going to hold on on other analyses (As this is kind of the last one to run) :)


### September 2, 2022 ###

**MSMC-IM**
I got MSMC-IM to work. The new version of MSMC2 could run the between script from above. 

First, combine runs:
```
#!/bin/bash
#SBATCH --job-name=generate_MSMC_input
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=00:30:00


~/scratch/PY_ENV/msmc_tools/combineCrossCoal.py ALL_BETWEEN_msmc2.final.txt ALL_NARW_msmc2.final.txt ALL_SRW_msmc2.final.txt > RW_combined_msmc.final.txt
```

The results from this can be used to plot relative CCR.

Second Run MSMC-IM:
```
#!/bin/bash
#SBATCH --job-name=MSMC-IM
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=00:30:00

MU="1.5e-8"
OUT="./msmc-im/rw_msmc_im_sept2"

~/scratch/PY_ENV/msmc_tools/MSMC_IM.py -mu $MU -o $OUT --xlog --ylog RW_combined_msmc.final.txt
```

I want to redo this with 10iterations of different samples.

Set 1 - Eau017 Eau283 EGL00252-1 EGL336_1b
Set 2 - Eau018 Eau7 EGL013-3qa SID179132
Set 3 - Eau023 Eau9c EGL140-1 EGL272-1
Set 4 - Eau029 Eau017 EGL183-1 EGL308-1a
Set 5 - Eau034A Eau018 EGL254-1 EGL00252-1
Set 6 - Eau10b Eau023 EGL272-1 EGL013-3qa
Set 7 - Eau283 Eau029 EGL308-1 EGL140-1
Set 8 - Eau7 Eau034A EGL336_1b EGL183-1
Set 9 - Eau9c Eau10b SID179132 EGL254-1
Set 10 - Eau017 Eau10b SID179132 EGL140-1

**MSMC**
I am rerunning the msmc 10 iterations with the new version of MSMC2.

**Stairway Plot**

First to create blueprint files:
```
java -cp stairway_plot_es Stairbuilder NARW_unfolded_Aug18.blueprint
java -cp stairway_plot_es Stairbuilder SRW_unfolded_Aug18.blueprint
```
Running Stairway plot now on both species.


### September 6, 2022 ###

I need to create a nice script (as it is pretty clever and might be called upon by other later) to generate the masking files and vcf file I need for MSMC inputs. More specifically, that I need for the generate_multihetsep.py script.



```
#!/bin/bash
#SBATCH --job-name=create_bed_vcf_inputs
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-378%20
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --time=1:00:00

################################
### The purpose of this script is to create the bed files required for generate_multihetsep.py for MSMC
###
### Requirements:
###		1) A file with two columns with each row being unique combination of sample name and scaffold 
###		2) A gvcf file for each individual, at each scaffold
###		3) Python script installed gvcf2bed
###		4) A phased vcf for each chromosome
###		5) A vcf of raw variants created after joint variant calling before filtering
###		6) vcftools/0.1.14
###		7) bcftools/1.11
###		8) bedtools/2.30.0
###
################################


###---INPUT FILES AND VARIABLES---###

# Individual and Scaffold/Chromosome List
FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" IND_SCAFFOLDS)
CHR=$(echo $FILE | awk '{print $1}')
SAMPLE=$(echo $FILE | awk '{print $2}')

# Parent directory for gvcf files (Mine are subset with species and scaffold, adjust below as needed)
GVCF_FILE=~/scratch/gvcf/SRWonNARW/HiC_scaffold_${CHR}/${SAMPLE}-HiC_scaffold_${CHR}.g.vcf.gz

# VCF file after joint variant calling, before filtering variants
RAW_VCF=~/projects/def-frasiert/RW_WGS/vcf/ALL/merged_all_on_narw.vcf.gz

# VCF file after all filtering and phasing - the VCF file you deemed 'locked in'
PHASED_VCF=~/scratch/msmc/shapeit/phased/HiC_scaffold_${CHR}_sorted.vcf.gz

# OUTPUT DIRECTORY
OUT_DIR=~/scratch/msmc/generate_multisephet_input/MSMC_IM_runs/

###---RUN CODE IN AN ARRAY WITH EACH JOB A DIFFERENT COMBINATION OF SCAFFOLD AND CHROMOSOME ---###

#create bed file from gvcf file to mark all sites that were able to be mapped
gvcf2bed -I ${GVCF_FILE} -O ${SAMPLE}_${CHR}.bed

#split the phased vcf file into individuals, remove indels and missing data.
vcftools --gzvcf ${PHASED_VCF} --indv ${SAMPLE} --max-missing 1.0 --remove-indels --recode --out ${SAMPLE}_${CHR}_all.filtered

#split the raw vcf file by individual and chromosome
vcftools --gzvcf ${RAW_VCF} --indv ${SAMPLE} --chr HiC_scaffold_${CHR} --recode --out ${SAMPLE}_${CHR}_all.originalmerge

#create a vcf file of sites where we called variants, but that were filtered out for any reason
bedtools subtract -a ${SAMPLE}_${CHR}_all.originalmerge.recode.vcf -b ${SAMPLE}_${CHR}_all.filtered.recode.vcf > ${SAMPLE}_${CHR}_all_filtered_sites.vcf

#correct header
bcftools view -h ${SAMPLE}_${CHR}_all.filtered.recode.vcf > temp-${SAMPLE}-${CHR}_all
cat temp-${SAMPLE}-${CHR}_all ${SAMPLE}_${CHR}_all_filtered_sites.vcf > ${SAMPLE}_${CHR}_all_filtered_sites_wheader.vcf

#use bedtools to create a bed file that removed the sites that were filtered out from the original bed file that listed all sites we could map as we are not confident these sites are not variants and we don't want to make this assumption
bedtools subtract -a ${SAMPLE}_${CHR}.bed -b ${SAMPLE}_${CHR}_all_filtered_sites_wheader.vcf > ${SAMPLE}_${CHR}_all_msmcmask.bed
gzip ${SAMPLE}_${CHR}_all_msmcmask.bed

#remove vcf file with raw variant information as it is no longer needed
rm ${SAMPLE}_${CHR}_all.originalmerge.recode.vcf
#remove tempory header and temporary vcf without correct header
rm temp-${SAMPLE}-${CHR}_all
rm ${SAMPLE}_${CHR}_all_filtered_sites.vcf

#Move good output files
mv ${SAMPLE}_${CHR}_all.filtered.recode.vcf ${OUT_DIR}
mv ${SAMPLE}_${CHR}_all_msmcmask.bed.gz ${OUT_DIR}
```

### September 7, 2022 ###

**Prepping for MSMC-IM**
This morning I zipped the ${SAMPLE}_${CHR}_all.filtered.recode.vcf files and started running the generate_multihetsep.py on all 10 sets as described above.

When these 210 jobs finish, I will run MSMC on the four NARW haplotypes, the four SRW haplotypes and the 16 pairwise comparisons for each set (30 jobs total set up in 3 arrays by set as the array number).

MSMC running with 30 jobs.


### September 8, 2022 ###

**MSMC-IM**
I moved these results (ALL*set*) into the results/msmc-im/.

I am going to join the combine msmc and msmc-im into a single script because they don't take long to run.

The script below is called combine_msmc-im.sh and is saved in the msmc directory and should be run from the results/msmc-im/.

```
#!/bin/bash
#SBATCH --job-name=Combine-MSMC-IM
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-10
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=01:00:00

SET=${SLURM_ARRAY_TASK_ID}

MU="0.9664e-8"
OUT=./msmc-im/rw_msmc_im_set${SET}


~/scratch/PY_ENV/msmc_tools/combineCrossCoal.py ALL_BETWEEN_msmc2_set${SET}.final.txt ALL_NARW_msmc2_set${SET}.final.txt ALL_SRW_msmc2_${SET}.final.txt > RW_combined_msmc_set${SET}.final.txt
~/scratch/PY_ENV/msmc_tools/MSMC_IM.py -mu $MU -o $OUT --xlog --ylog RW_combined_msmc_set${SET}.final.txt

```

I can plot the rCCR from RW_combined_msmc_set${SET}.final.txt and plot the M(t) curve with time from rw_msmc_in_set${SET}*
