## **Read Mapping and Variant Calling- RNA version**

### 1. Read Quality Assessment
#### [FastQC v 0.11.9](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
 With any sequencing project, it is always a good idea to first assess the quality of the data we generated. We can use FastQC for this purpose.

Let look at the _fastq_qual.sh_ script:
```
# load modules
module load FastQC/0.11.9-Java-11

# change me!
DIR=/scratch/group/kitchen-group/MARB_689_Molecular_Ecology/class_working_directories
USER=kitchens
DIR2=Project

#create a new directory for the fastqc report
mkdir ${DIR}/${USER}/${DIR2}/fastqc_report

#run fastqc
for SAMPLE in $(cat SRA.list); do
        fastqc -o ${DIR}/${USER}/${DIR2}/fastqc_report -f fastq -t 4 \
        ${DIR}/${USER}/${DIR2}/${SAMPLE}_1.fastq \
        ${DIR}/${USER}/${DIR2}/${SAMPLE}_2.fastq
done
```
-o = output directory
-f = fastq reads (can be fastq, bam or sam)
*Approximate run time for all data from a single sample = 2-5 mins*

Then copy the html output to your local computer (Mac version below):
```
scp username@login.hprc.tamu.edu:/scratch/group/kitchen-group/MARB_689_Molecular_Ecology/class_working_directories/kitchens/fastqc_report/example_1_fastqc.html ./
```

The html provides the summary statistics and various graphical representations of the data (Fig 2).

<p align="center"> <img src="./fastqc_ATT03.png" width="500" />
Figure 2. Sequence quality per base.

### 2. Read filtering and trimming
Due to the adapter contamination found in the raw reads we need to clean them up. We can remove low quality base pairs at the tails of the reads as well as the Illumina adapters using [cutadapt v5.0](https://cutadapt.readthedocs.io/en/stable/index.html).

```
#load modules
module load GCCcore/13.2.0
module load cutadapt/5.0

# Change me!
DIR=/scratch/group/kitchen-group/MARB_689_Molecular_Ecology/class_working_directories
USER=kitchens
DIR2=Project

# run cutadapt v5.0
for SAMPLE in $(cat SRA.list); do
        cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 50 -q 15 -j 0 \
        -o ${DIR}/${USER}/${DIR2}/${SAMPLE}_filtered_1.fastq \
        -p ${DIR}/${USER}/${DIR2}/${SAMPLE}_filtered_2.fastq \
        ${DIR}/${USER}/${DIR2}/${SAMPLE}_1.fastq \
        ${DIR}/${USER}/${DIR2}/${SAMPLE}_2.fastq
done
```
-a/-A = adapter sequence, if both are present this is paired end data
-m = minimum length, reads shorter than defined length will be discarded
-q = quality trim ends of the reads at 3' end
-j = number of cores to use, if set to 0 automatically detect available num. of threads
-o = filename of the 1st read
-p = filename of the 2nd read

Recall that _**fastq**_ reads provide the sequence data from the Illumina cluster generation and quality scores for each read.
```
@HISEQ:1050:HYJCFBCX2:1:1101:2228:2472 1:N:0:ACTGAT
TTTGAAAATCGTAAAAATTGAAATATACAATTATGAACTATTTCGATATGCACTGTTAAAAATGCTACCAGATTTTTAAAAAAATATTTTAAAGCCAAAGACTTATCTTTTAAACAATGCATTATTATATTGCAAAAGTTTATCTTTTAA
+
DDDDDHIIIIIHIIIIIIIIIIIIIIIIIIIIIIHIIIIHIIHHIIIIHIIIIIIHIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIHIIIHGHIHIIIIIIIIIIIIIIIIIHIIIIIIIHHIIIIIIHHIFHIIIIIIIIII
@HISEQ:1050:HYJCFBCX2:1:1101:4107:2475 1:N:0:ACTGAT
TGTCATTTCATGTCATTTCATCTCATGTGCTATGTCCTGTCCTGTCCTGTCCTGTCATGTCATGTCATGTCATGTTATTTCATGTCATTTCATGTCATATGCTATGTCCTGTCCTGTCCTGTCATGTCATGTCATGTCATGTCCTGTCCT
+
DDDDDIIIIIIIIIIIIIIIIIIIIIIIIIIIHIIIIIIIIIIIIIIIHIIIHIII1CGHI1FHIH1DGHI?HHI11<1DG<FHEE1F1FH1DGHH1C<F11<1<1<<<F1<1<C1<1<D1D<1<<C11<<<<<D1D1<D<C1<D1<<<G
```

### 3. Align reads to the reference genome or transcriptome
Mapping of RNAseq reads against a reference (genome or transcriptome) is one of the key steps of the analysis. This is where we deviate from the DNA protocol because RNAseq reads only align to exons of the genome, and therefore we need to consider reads that span splice junctions. There are also tools like *kallisto* and *salmon* that "pseudoalign" to a reference transcriptome without actually aligning the reads. These are great for differential expression but not great for calling SNPs.  Below is a list of common mapping tools available that you might come across in the literature:

- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)- not splice aware
- [STAR](https://github.com/alexdobin/STAR)
- [HISAT2](https://daehwankimlab.github.io/hisat2/)

These read mappers compare the reads against a reference sequence that has been transformed into a highly accessible data structure called a genome/transcriptome index. Such indexes should be generated before mapping begins. In our case, we will use STAR.

We have already created a genome index as follows (`nano star_index.sh`):
```
#load modules
module load GCC/12.2.0
module load STAR/2.7.10b

# Change me!
REFERENCE=Cassiopea

STAR --runMode genomeGenerate --genomeDir /scratch/group/kitchen-group/MARB_689_Molecular_Ecology/map_reference/ \
--genomeFastaFiles ${REFERENCE}_genome.fa --genomeSAindexNbases 13 \
--outFileNamePrefix ${REFERENCE}
```

Now, we will run STAR tool on one sample like:
```
#load modules
module load GCC/12.2.0
module load STAR/2.7.10b
module load SAMtools/1.17

# Change me!
DIR=/scratch/group/kitchen-group/MARB_689_Molecular_Ecology/class_working_directories
DIR2=Project
USER=kitchens
REFERENCE=Cassiopea
DIR3=/scratch/user/${USER}

# command
cd ${DIR3}

for SAMPLE in $(cat ${DIR}/${USER}/${DIR2}/SRA.list); do
        STAR --genomeDir /scratch/group/kitchen-group/map_reference/${REFERENCE}/ \
        --runThreadN 4 --readFilesIn ${DIR}/${USER}/${DIR2}/${SAMPLE}_filtered_1.fastq ${DIR}/${USER}/${DIR2}/${SAMPLE}_filtered_2.fastq \
        --outFileNamePrefix ${SAMPLE}_  --twopassMode Basic

        # how good was the mapping
        samtools flagstat ${DIR3}/${SAMPLE}_Aligned.out.sam
done
```

The number of threads to use is provided by "--runThreads". To increase sensitivity of the alignments to the exon-exon junctions, we will run the "--twopassMode". In this mode, the genome indices are re-generated from splice junctions obtained from a 1-pass mode with the usual parameters and then run the mapping step (2-pass mapping).

At the end, we run [samtools flagstat tool](https://www.htslib.org/doc/samtools.html) to summarize the number of mapped reads to each category: primary, secondary, supplementary, duplicates, properly paired, singletons.

__HOMEWORK__: Make a table of number of all reads, filtered reads, mapped reads, and properly paired mapped reads for your final report.

### 4. Remove PCR Duplicates
After the tool completes, the alignments will be stored in the Sequence Alignment Map (SAM) format, a standard for storing aligned reads. The SAM format is a generic nucleotide alignment format that describes the alignment of sequencing reads (or query sequences) to a reference. We will convert the SAM file into the binary form of the format (BAM), which is compact and can be rapidly searched (if indexed).

We will use [*Picard Tools*](https://broadinstitute.github.io/picard/command-line-overview.html) to sort the SAM file by coordinates and convert it into the BAM format.

As a way to track the samples, we will assign a read group identifier using the "AddOrReplaceReadGroups" tool in Picard Tools. The read group (RG) assignment for each sample requires the following information:

  - ID= read group identifier, generally takes on the flowcell + lane name and number
  - SM= sample, the name of the sample sequenced in the read group
  - PL= platform, valid values include ILLUMINA, SOLID, LS454, HELICOS and PACBIO
  - LB= library identifier, MarkDuplicates uses the LB field to determine which read groups might contain molecular duplicates, in case the same DNA library was sequenced on multiple lanes

Preparation of sequencing libraries for technologies such as Illumina (used in your study) typically involves PCR amplification. Duplicates can be identified based on their outer alignment coordinates or using sequence-based clustering. One of the common ways for identification of duplicate reads is `MarkDuplicates` utility from *Picard Tools* package. It is designed to identify both PCR and optical duplicates. It is not recommended that we actually remove the duplicate sequences from the file, but simply to mark the flags appropriately in the BAM file, so that those sequences are ignored downstream. At this step you can merge together any BAM files from the same sample that are currently separated.

Run MarkDuplicates tool on the BAM files. This tool with have two outputs, a new BAM file with duplicates marked and a metric file about the original BAM. Once again we will sort and index the resulting BAM file.

```
#load modules
module load picard/2.25.1-Java-11

# change me!
DIR=/scratch/group/kitchen-group/MARB_689_Molecular_Ecology/class_working_directories
DIR2=Project
USER=kitchens
DIR3=/scratch/user/${USER}

# add read group and sort the file
for SAMPLE in $(cat SRA.list); do
        java -Xmx16g -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
        I=${DIR3}/${SAMPLE}_Aligned.out.sam O=${DIR3}/${SAMPLE}_rg_added_sorted.bam \
        SO=coordinate \
        RGID=MARB689 RGLB=${SAMPLE} RGPL=ILLUMINA RGPM=HISEQ RGSM=${SAMPLE}

        # index the resulting bam file
        java -Xmx16g -jar $EBROOTPICARD/picard.jar BuildBamIndex \
        I=${DIR3}/${SAMPLE}_rg_added_sorted.bam

        # remove duplicates
        java -Xmx16g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        TMP_DIR=${DIR3}/ \
        I=${DIR3}/${SAMPLE}_rg_added_sorted.bam \
        O=${DIR3}/${SAMPLE}_dedup.bam \
        METRICS_FILE=${DIR3}/${SAMPLE}_dedup.metrics_test.txt \
        REMOVE_DUPLICATES=false \
        TAGGING_POLICY=All \
        CREATE_INDEX=true

        # remove prior alignment files
        rm ${DIR3}/${SAMPLE}_Aligned.out.sam
        rm ${DIR3}/${SAMPLE}_rg_added_sorted.bam
done
```
