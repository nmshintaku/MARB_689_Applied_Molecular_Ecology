#!/bin/bash
#SBATCH --time=02:00:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=48G   # memory per CPU core
#SBATCH --job-name "mtA" # job name
#SBATCH --output "mtA.log" # log file name

# Change the information below
USER=kitchens
SPECIMEN=A

# Don't touch

perl /scratch/group/kitchen-group/02_mitogenome/NOVOPlasty4.3.3/NOVOPlasty4.3.3.pl -c config_specimen${SPECIMEN}.txt

