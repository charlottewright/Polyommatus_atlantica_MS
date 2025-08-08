#!/bin/bash

GENOME=$1
BAM=$2
PREFIX=$3
NUM_CPU=$4

# run deepvariant

# example command:
#/software/singularity-v3.6.4/bin/singularity exec --bind /usr/lib/locale/ --bind `pwd -P` docker://google/deepvariant:1.4.0 /opt/deepvariant/bin/run_deepvariant --model_type PACBIO -ref Polyommatus_icarus.fasta --reads Polyommatus_icarus.sorted_dups.bam --output_vcf deepvariant_output/Polyommatus_icarus.vcf.gz --num_shards 20
singularity exec --bind /usr/lib/locale/ \
--bind /lustre/scratch122/tol/teams/blaxter/projects/lepidoptera/cw22/polyommatus_atlantica/Analysis/mapped_pacbio/ \
--bind /lustre/scratch122/tol/teams/blaxter/projects/lepidoptera/cw22/polyommatus_atlantica/Raw_data/genomes/core_dataset/ \
--bind `pwd -P` docker://google/deepvariant:1.4.0 /opt/deepvariant/bin/run_deepvariant --model_type PACBIO -ref $GENOME --reads $BAM --output_vcf ${PREFIX}.vcf.gz --num_shards $NUM_CPU
