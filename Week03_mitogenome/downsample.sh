#!/bin/bash
#SBATCH --time=00:30:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=10G   # memory per CPU core
#SBATCH --output "DS.log"

# module loading
module load GCC/11.2.0
module load seqtk/1.3

# Change the information below
USER=kitchens
SPECIMEN=A
sample=5000000

# Don't touch

# downsample reads
seqtk sample -s 10 /scratch/group/kitchen-group/02_mitogenome/${SPECIMEN}_R1.fastq.gz ${sample} \
> /scratch/group/kitchen-group/class_working_directories/${USER}/${SPECIMEN}_sub_R1.fastq

seqtk sample -s 10 /scratch/group/kitchen-group/02_mitogenome/${SPECIMEN}_R2.fastq.gz ${sample} \
> /scratch/group/kitchen-group/class_working_directories/${USER}/${SPECIMEN}_sub_R2.fastq

# module change
module unload GCC/11.2.0
module load GCCcore/11.2.0
module load cutadapt/3.5

# run cutadapt 
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 50 -q 15 -j 0 \
-o /scratch/group/kitchen-group/class_working_directories/${USER}/${SPECIMEN}_subset_filtered_R1.fastq \
-p /scratch/group/kitchen-group/class_working_directories/${USER}/${SPECIMEN}_subset_filtered_R2.fastq \
/scratch/group/kitchen-group/class_working_directories/${USER}/${SPECIMEN}_sub_R1.fastq \
/scratch/group/kitchen-group/class_working_directories/${USER}/${SPECIMEN}_sub_R2.fastq

