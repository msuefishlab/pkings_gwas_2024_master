library(GenomicFeatures)

txdb <- makeTxDbFromGFF("paramormyrops_ncbi_0.4.gff")
k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k, columns = "TXNAME", keytype = "GENEID")
tx2gene <- df[, 2:1]
head(tx2gene)
write.csv(tx2gene, "paramormyrops_ncbi_0.4_tx2gene.csv")
