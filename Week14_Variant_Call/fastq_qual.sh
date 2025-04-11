#!/bin/bash
#SBATCH --time=12:00:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=10G   # memory per CPU core
#SBATCH -J "fqc" # job name
#SBATCH --output "fqc"

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
	${DIR}/${USER}/${DIR2}/${SAMPLE}_1.fastq.gz \
	${DIR}/${USER}/${DIR2}/${SAMPLE}_2.fastq.gz
done

