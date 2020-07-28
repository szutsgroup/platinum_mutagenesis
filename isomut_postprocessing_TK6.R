# Aim: find a Fischer score tresholds where at most 5 SNVs, at most 1 deletions and at most 1 insertions are left in any of the starting clones

library(tidyverse)
library(reshape2)
library(BSgenome.Hsapiens.NCBI.GRCh38)


treatments <- list(
  starting_clone = paste0("TK", c(4, 5)),
  mock = paste0("TK", c(11:13)),
  cisplatin = c("TK45", "TK48", "TK49"),
  carboplatin = c("TK46", "TK50", "TK51"),
  oxaliplatin = c("TK47", "TK52", "TK53"),
  other_samples = paste0("TK", c(1, 3, 29:37, 38:44))
)

# genotypes
genotypes <- list(
  WT = unname(unlist(treat_list))
)


#reading in isomut SNV output file
snv <- read.delim("isomut_output/all_SNVs_TK6.isomut")
snv$sample <- substr(snv$X.sample_name, 0, 5)

#cumulative plot: SNV counts vs. Fischer score
csnv <- matrix(0, ncol = length(table(snv$sample)), nrow = 301)
for (i in 1:301) { 
  csnv[i,] <- as.vector(tabulate(factor(filter(snv, score > (i-1)/10)$X.sample_name, levels = names(table(snv$X.sample_name))), nbins = length(table(snv$sample)))) 
}
csnv <- data.frame(csnv, stringsAsFactors = FALSE)
names(csnv) <- names(table(snv$sample))
csnv$treshold <- seq(0, 30, by = 0.1)
csnvm <- melt(csnv, id = "treshold")
csnvm$status <- ifelse(csnvm$variable %in% treat_list$starting_clone, "starting clone", "other")
ggplot(data = csnvm, aes(x = treshold, y = log(value, base = 10), group = variable, color = status)) + geom_line() + ylab("log10 (Mutation number) ") + xlab("Score treshold") + geom_hline(yintercept = log10(5), linetype = 2)

#calculate suitable Fischer treshold
nvScore <- csnv[ ,c(treat_list$starting_clone, "treshold")] %>% 
  apply(1, function(x) ifelse(sum(x[1:length(treat_list$starting_clone)] < 5) == length(treat_list$starting_clone), x[length(treat_list$starting_clone) + 1], NA)) %>% 
  min(na.rm = TRUE)

#filtered snv set
snv_filt <- filter(snv, score > snvScore, chr != "MT")
snv_filt$genotype <- NA
for (i in 1:length(geno_list)) {
  snv_filt[which(snv_filt$sample %in% geno_list[[i]]), "genotype"] <- names(geno_list)[i]
}
snv_filt$treatment <- NA
for (i in 1:length(treat_list)) {
  snv_filt[which(snv_filt$sample %in% treat_list[[i]]), "treatment"] <- names(treat_list)[i]
}



#reading in isomut indel output file
indel <- read.delim(file.path(outDir, "all_indels_TK6.isomut"))
indel$sample <- substr(indel$X.sample_name, 0, 5)

#Fischer scores separately for insertions and deletions
#first insertions
ins <- filter(indel, type == "INS")
cins <- matrix(0, ncol = samplenum, nrow = 201)
for (i in 1:201) {
  cins[i,] <- as.vector(tabulate(factor(filter(ins, score > (i-1)/10)$sample, levels = names(table(indel$sample))), nbins = samplenum)) 
}
cins <- data.frame(cins, stringsAsFactors = FALSE)
names(cins) <- names(table(indel$sample))
cins$treshold <- seq(0, 20, by = 0.1)
cinsm <- melt(cins, id = "treshold")
cinsm$status <- ifelse(cinsm$variable %in% setdiff(treat_list$starting_clone, problematic), "starting clone", "other")
ggplot(data = cinsm, aes(x = treshold, y = log(value, base = 10), group = variable, color = status)) + geom_line() + ylab("log10 (Insertion number)") + xlab("Score treshold") + geom_hline(yintercept = log10(1), linetype = 2)

#deletions
del <- filter(indel, type == "DEL")
cdel <- matrix(0, ncol = samplenum, nrow = 201)
for (i in 1:201) { 
  cdel[i,] <- as.vector(tabulate(factor(filter(del, score > (i-1)/10)$sample, levels = names(table(indel$sample))), nbins = samplenum)) 
}
cdel <- data.frame(cdel, stringsAsFactors = FALSE)
names(cdel) <- names(table(indel$sample))
cdel$treshold <- seq(0, 20, by = 0.1)
cdelm <- melt(cdel, id = "treshold")
cdelm$status <- ifelse(cdelm$variable %in% setdiff(treat_list$starting_clone, problematic), "starting clone", "other")
ggplot(data = cdelm, aes(x = treshold, y = log(value, base = 10), group = variable, color = status)) + geom_line() + ylab("log10 (Deletion number)") + xlab("Score treshold") + geom_hline(yintercept = log10(1), linetype = 2)

write.table(snv_filt, file = "all_SNVs_postprocessed_TK6.isomut", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(del_filt, file = "all_DELs_postprocessed_TK6.isomut", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ins_filt, file = "all_INSs_postprocessed_TK6.isomut", sep = "\t", row.names = FALSE, quote = FALSE)

# triplet spectra

snv_filt$chr <- gsub("MT", "M", snv_filt$chr)
snv_filt$chr <- factor(snv_filt$chr, levels = c(1:22, "X", "Y"))
snv_range <- GRanges(snv_filt$chr, ranges = IRanges(snv_filt$pos - 1, end = snv_filt$pos + 1))
snv_seq <- getSeq(Hsapiens, snv_range)
snv_filt$tripmut <- paste0(as.character(snv_seq), snv_filt$mut)
tmp <- DNAStringSet(snv_filt$mut)
snv_filt$tripmut96 <- ifelse((snv_filt$ref == "C" |snv_filt$ref == "T"), 
                             snv_filt$tripmut, 
                             paste0(unname(as.character(reverseComplement(snv_seq))), 
                                    unname(as.character(reverseComplement(tmp)))
                             )
)
snv_filt$sorting=paste0(substr(snv_filt$tripmut96, 2, 2), 
                        ">",
                        substr(snv_filt$tripmut96, 4, 4),
                        "::", 
                        substr(snv_filt$tripmut96, 1, 1), 
                        "_", 
                        substr(snv_filt$tripmut96, 3, 3))

snv_spectra <- as.data.frame.matrix(table(snv_filt$sample, snv_filt$sorting))
write.table(snv_spectra, "TK6_triplets.csv")
