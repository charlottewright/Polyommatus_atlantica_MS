 #!/usr/bin/env bash

PATH=$PATH:/lustre/scratch123/tol/teams/blaxter/users/cw22/software/bedtools/bedtools2/bin/
PATH=$PATH:/lustre/scratch123/tol/teams/blaxter/users/cw22/software/bedops/

BREAKPOINT_FILE=$1 # ilCyaSemi1.1_vs_ilLysCori1.1_breakpoints_refined_231123.tsv
INDEX_FILE=$2 # Cyaniris_semiargus.fasta.fai
PREFIX=$3 # e.g. Cyaniris_semiargus_allRepeats
REPEAT_GFF=$4
QUERY_SPECIES=$5
# First generate a set of random regions per chr
# Specify:
# breakpoint_regions file
# Index file of species
# output_prefix
python3 ./random_region_generator.py -b $BREAKPOINT_FILE -i $INDEX_FILE -o $PREFIX -n 10000

# Next calculate the density of repeats per random region
# Specify:
# locations of random regions
# GFF of repeats
bedtools coverage -a locations_of_random_regions_10000_per_chr_${PREFIX}.tsv -b ${REPEAT_GFF} > repeat_density_of_random_regions_10000_${REPEAT_GFF}.tsv

# Next calculate repeat density per breakpoint region
# Need to swap column 4 and 5 if the end is greater than the start value (happens occasionally in ilPolAtal1)
cat $BREAKPOINT_FILE | awk -F'\t' '{if ($2 > $3) {temp=$2; $2=$3; $3=temp} print}' OFS='\t' | awk -v OFS='\t' '{print $1, $2, $3}' >${BREAKPOINT_FILE}.bed
#tail -n +2 $BREAKPOINT_FILE | awk -F'\t' '{if ($4 > $5) {temp=$4; $4=$5; $5=temp} print}' OFS='\t' | awk -v OFS='\t' '{print $1, $4, $5}' >${BREAKPOINT_FILE}.bed
bedtools coverage -a ${BREAKPOINT_FILE}.bed -b $REPEAT_GFF > repeat_density_of_${BREAKPOINT_FILE}_${REPEAT_GFF}.tsv

# Now we can compare the repeat density in random regions to those of breakpoint regions
# Specify:
# repeat_density_of_random_regions_Cyaniris_semiargus.tsv
# repeat_density_of_breakpoints_refined_231123_Cyaniris_semiargus.tsv
python3 ./compare_feature_density_in_breakpoint_regions_to_random.py -b repeat_density_of_${BREAKPOINT_FILE}_${REPEAT_GFF}.tsv -r repeat_density_of_random_regions_10000_${REPEAT_GFF}.tsv -o ${PREFIX}_${QUERY_SPECIES}
