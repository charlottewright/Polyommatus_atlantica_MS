# Run command on one file  $FILE
#!/bin/bash 

# Setup environment
# conda activate busco_env
# The busco cut command is needed to make sure the busco dir names for each genome don't contain '.fasta' as this causes problems with the busco script

FASTA=$1
PREFIX=$2
CPU=$3

echo "Processing ${FASTA}"
busco -i $FASTA \
-l /software/busco_downloads/lineages/lepidoptera_odb10 \
-o busco_lepidoptera.augustus.${PREFIX} \
-c $CPU \
--mode genome \
--augustus \
--offline
