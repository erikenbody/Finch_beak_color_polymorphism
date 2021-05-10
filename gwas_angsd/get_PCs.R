#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

outname <- gsub("pcangsd.cov","PC1_PC2.txt", args[1])
pca.sample <- read.table(args[1])
eigen.sample <- eigen(pca.sample, symmetric = T)
eigenvectors.sample <- as.data.frame(eigen.sample$vectors)
write.table(eigenvectors.sample[,c(1,2)], outname, col.names = F, quote = F, row.names = F, sep = "\t")
