#!/bin/bash
#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=2   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=20G   # memory per CPU core
#SBATCH -J "ca"   # job name
#SBATCH --output "ca"

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

