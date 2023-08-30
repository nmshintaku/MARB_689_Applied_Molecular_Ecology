#!/bin/bash
#SBATCH --time=00:30:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=32G   # memory per CPU core
#SBATCH --output "nt.log"

module load GCC/11.2.0
module load OpenMPI/4.1.1
module load BLAST+/2.12.0

# Change the information below
USER=kitchens
SPECIMEN=A
sequence=Seed.fasta

# Don't touch
# run blast
blastn -query ${sequence} -db /scratch/data/bio/blast/nt -max_target_seqs 10 -num_threads 4 \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle staxids" > report_specimen${SPECIMEN}.txt

# sort blast results for top hit by bitscore and e-value
sort -k1,1 -k12,12nr -k11,11n report_specimen${SPECIMEN}.txt | sort -u -k1,1 --merge > tophit_specimen${SPECIMEN}.txt


