# 
# source("https://bioconductor.org/biocLite.R")
# biocLite("dada2")

library(dada2); 
packageVersion("dada2")
## [1] '1.0.3'
library(ShortRead); 
packageVersion("ShortRead")
## [1] '1.31.1'
library(ggplot2); 
packageVersion("ggplot2")

path <- "~/tools_dada2/MiSeq_SOP" 
fns <- list.files(path)
fns

# filter and trim
fastqs <- fns[grepl(".fastq$", fns)]
fastqs <- sort(fastqs)

fnFs <- fastqs[grepl("_R1", fastqs)]
fnRs <- fastqs[grepl("_R2", fastqs)]

sample.names <- sapply(strsplit(fnFs, "_"), "[", 1)

fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

# quality

plotQualityProfile(fnFs[[1]])

# filtering and triming

# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
# Filter
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]),
                    trimLeft=c(10, 10), truncLen=c(240,160), 
                    maxN=0, maxEE=2, truncQ=2, 
                    compress=TRUE, verbose=TRUE)
}

# dereplication

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

derepFs[[1]]

# sample infer

dadaFs <- dada(derepFs, err = NULL, selfConsist = TRUE)
dadaRs <- dada(derepRs, err = NULL, selfConsist = TRUE)

dadaFs[[1]]
plotErrors(dadaFs[[1]], nominalQ = TRUE)

# merge

mergers <- mergePairs(dadaFs, derepFs, dadaRs,derepRs, verbose = TRUE)
head(mergers[[1]])

# seq table
# similar to OTU table

seqtab <- makeSequenceTable(mergers[names(mergers) != "mock"])
dim(seqtab)
table(nchar(colnames(seqtab)))
# seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(230,236)])

# remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, verbose = TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim) / sum(seqtab)

# assign taxa

taxa <- assignTaxonomy(seqtab.nochim, refFasta = "rdp_train_set_14.fa.gz")
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
unname(head(taxa))
# New in 1.1: Accurate species-level assignment using exact matching with assignSpecies


# RData
save.image("dada2.RData")
