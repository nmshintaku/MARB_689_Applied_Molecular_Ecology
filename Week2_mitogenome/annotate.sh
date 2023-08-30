#!/bin/bash
#SBATCH --time=00:30:00   # walltime
#SBATCH --ntasks=8   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=32G   # memory per CPU core
#SBATCH --output "annoA.log"

# load modules
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/2.7.15
module load Java/17.0.2
module load BLAST+/2.8.1

# Change the information below
USER=kitchens
SPECIMEN=A

# Don't touch
/scratch/group/kitchen-group/02_mitogenome/MitoFinder/mitofinder -j ${SPECIMEN}_unknown \
-a /scratch/group/kitchen-group/02_mitogenome/specimen${SPECIMEN}_mitogenome.fasta \
-r /scratch/group/kitchen-group/02_mitogenome/reference.gb -o 5 -p 8 -m 32


