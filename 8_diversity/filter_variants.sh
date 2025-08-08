#!/bin/bash

VCF=$1
REPEAT_GFF=$2
PREFIX=$3
NUM_THREADS=$4

# filter to remove RefCall variants
grep -v 'RefCall' $VCF > ${PREFIX}.PASS.vcf

# filter to only keep biallelic SNPs
bcftools view -m2 -M2 -v snps ${PREFIX}.PASS.vcf --threads $NUM_THREADS > ${PREFIX}.PASS.biallelic.vcf

# filter to remove SNPs in repetitive regions
bedtools subtract -header -a ${PREFIX}.PASS.biallelic.vcf -b $REPEAT_GFF > ${PREFIX}.deepvariant.PASS.biallelic.repeat_filtered.vcf
