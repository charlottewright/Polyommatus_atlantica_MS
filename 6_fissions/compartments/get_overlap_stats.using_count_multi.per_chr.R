#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressMessages(library(regioneR))
suppressMessages(library(rtracklayer))

set.seed(12345)

nuwts<-import(args[1])
genes<-import(args[2])

print('Number of overlaps:')
numOverlaps(nuwts, genes, count.once=TRUE)
print('Number of overlaps in total:')
numOverlaps(nuwts, genes)

genome<-read.table(args[3])
#test<-randomizeRegions(genes, genome=genome, per.chromosome=TRUE, allow.overlaps=FALSE)
print(args[6])
pt<-permTest(A=nuwts, B=genes, ntimes=as.numeric(args[6]), genome=genome, randomize.function = randomizeRegions, evaluate.function = numOverlaps,per.chromosome=TRUE)
#pt<-overlapPermTest(A=nuwts, B=genes, ntimes=as.numeric(args[6]), genome=genome, non.overlapping=FALSE, per.chromosome=TRUE,mc.set.seed=FALSE, mc.cores=10,count.once=TRUE)
pdfname=sprintf("%s.%s.pdf", args[4],args[5])
pdf(pdfname)
plot(pt)
dev.off()

line=sprintf("%s,%s,%s,%s", args[4], args[5], pt$numOverlaps$pval,pt$numOverlaps$zscore)
filename_out=sprintf("%s.%s.overlap_stats.txt", args[4],args[5])
write(line,file=filename_out)
