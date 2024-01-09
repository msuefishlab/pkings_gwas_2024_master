library(getopt)
getwd()
source('code/05_PopGen/gdsSubset.R')

h <- function(x) {
    cat("Usage: Rscript --vanilla convert_GDS_to_PLINK.R -i VCF_FILE -o PLINK_FILE_PREFIX [-h]\n\n")
    quit(save="no", status=x)
}

opt <- getopt(matrix(c(
    'help',   'h', 0, "logical",
    'output', 'o', 1, "character",
    'input', 'i', 1, "character"
), ncol=4, byrow=TRUE));

if (! is.null(opt$help)) { h(0) }

file.ouptut <- ifelse(is.null(opt$input), "output", opt$output)

if (! is.null(opt$input)) {
    gds.file <- opt$input
    if (! file.exists(gds.file)) { cat(sprintf("GDS file (%s) was not found!\n", gds.file)); h(1) }
} else {
    h(1)
}

require(SNPRelate)

gdsfile_new<-snpgdsOpen(gds.file)
snpset<-snpgdsLDpruning(gdsfile_new,ld.threshold=0.1,maf=0.1,missing.rate=0.1,num.thread=2)
snpset.id<-unlist(snpset)
snpgdsGDS2BED(gdsfile_new, opt$output, sample.id=NULL, snp.id=snpset.id, snpfirstdim=NULL,verbose=TRUE)
snpgdsClose(gdsfile_new)
