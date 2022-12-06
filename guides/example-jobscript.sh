#!/bin/bash

# Request a specific partition:
# The partitions available are: batch, gpu, and bigmem
#SBATCH --partition=batch

# Configure runtime, memory usage, and the number of CPUs:
#SBATCH --time=48:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=32

# Notify the user details about the job:
#SBATCH --mail-user=example@mail.com
#SBATCH --mail-type=ALL

# Paths to bind to container
project_resources=heome
fastq_dir=fastq
wrangler_dir=wrangler

# Wrangler options
probe_sets_used='HeOME96'
sample_sets_used='JJJ'
experiment_id='example_id'
sample_list='sample_list.tsv'
min_capture_length=30

singularity run \
  -B ${project_resources}:/opt/project_resources \
  -B ${fastq_dir}:/opt/data \
  -B ${wrangler_dir}:/opt/analysis \
  --app wrangler miptools.sif \
  -p ${probe_sets_used} -s ${sample_sets_used} -e ${experiment_id} \
  -l ${sample_list} -c ${SLURM_CPUS_PER_TASK} -m ${min_capture_length}