# This script is adapted from: https://github.com/lstevens17/heligmosomoides_MS/blob/main/4_divergent_haplotypes/go_term_enrichment/run_topgo.R

args <- commandArgs(trailingOnly = TRUE)
universeFile = "gene2go.tsv"
interestingGenesFile = args[1]
output_file = args[2]
ont = args[3]

# set the output file
sink(output_file)

# load topGO
library("topGO")

# read in the 'gene universe' file
geneID2GO <- readMappings(file = universeFile)
print('yay file1!')
geneUniverse <- names(geneID2GO)

# read in the genes of interest
genesOfInterest <- read.table(interestingGenesFile,header=FALSE)
genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

# build the GOdata object in topGO
myGOdata <- new("topGOdata", description="My project", ontology=ont, allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
myGOdata

# run the Fisher's exact tests
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")
resultElim <- runTest(myGOdata, algorithm="elim", statistic="fisher")
resultTopgo <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
resultParentchild <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")

# see how many results we get where weight01 gives a P-value <= 0.001:
mysummary <- summary(attributes(resultTopgo)$score <= 0.05)
numsignif <- as.integer(mysummary[[3]]) 
numsignif

# print out the top 'numsignif' results:
allRes <- GenTable(myGOdata, classicFisher = resultClassic, elimFisher = resultElim, topgoFisher = resultTopgo, parentchildFisher = resultParentchild, orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif)
allRes

write.table(allRes, file=paste(args[1], "_", args[3], ".tsv", sep =""), row.names=FALSE, sep="\t", quote = FALSE)

# print out the genes that are annotated with the significantly enriched GO terms:
myterms <- allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
for (i in 1:length(myterms))
{
  myterm <- myterms[i]
  mygenesforterm <- mygenes[myterm][[1]]
  myfactor <- mygenesforterm %in% genesOfInterest # find the genes that are in the list of genes of interest
  mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
  mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
  print(paste("Term",myterm,"genes:",mygenesforterm2))
}
# close the output file
#sink() 
 
# print out the genes that are annotated with the significantly enriched GO terms:
myterms <- allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
for (i in 1:length(myterms))
{
   myterm <- myterms[i]
   mygenesforterm <- mygenes[myterm][[1]]
   myfactor <- mygenesforterm %in% genesOfInterest # find the genes that are in the list of genes of interest
   mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
   mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
   print(paste("Term",myterm,"genes:",mygenesforterm2))
}
# close the output file
sink() 
