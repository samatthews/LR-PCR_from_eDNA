#!/bin/bash
#SBATCH --job-name=lrpcr_blast
#SBATCH --account=coenv
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=lrpcr_blast_%j.out        # %j becomes the job ID
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=samatt@uw.edu



apptainer exec --bind /gscratch/ \
  --env DB=/gscratch/coenv/rpkelly/vert_mtDNA/vert_mtDNA,WORD_SIZE=64,THREADS=${SLURM_CPUS_PER_TASK} \
  /gscratch/coenv/rpkelly/apptainers/kellylab.sif \
  bash collab_blast.sh /mmfs1/home/samatt/gscratch/lrpcr_reads /gscratch/coenv/samatt/blast_lrpcr_out_newdb
