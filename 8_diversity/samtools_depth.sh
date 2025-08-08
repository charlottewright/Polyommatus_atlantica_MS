#!/bin/bash

# Author: @cb46

bam=$1

samtools depth $bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > $bam".cov"
