# Aim: find a Fischer score tresholds where at most 5 SNVs, at most 1 deletions and at most 1 insertions are left in any of the starting clones

library(tidyverse)
library(reshape2)
library(BSgenome.Ggallus.UCSC.galGal4)


#treatments
treat_list <- list(
  starting_clone = c("PL01", "PL02", "PL03"),
  mock = paste0("PL0", 4:6),
  cisplatin_7.5uM = sprintf("PL%03d", 7:11),
  cisplatin_10uM = paste0("PL", 12:16),
  carboplatin_75uM = paste0("PL", 17:21),
  carboplatin_100uM = paste0("PL", 22:26),
  oxaliplatin_2.5uM = paste0("PL", 27:31),
  oxaliplatin_4.5uM = paste0("PL", 32:36)
)


# genotypes
geno_list <- list(
  WT1 = unname(unlist(treat_list))
)


#reading in isomut SNV output file
snv <- read.delim("isomut_output/all_SNVs_DT40.isomut")
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
indel <- read.delim(file.path(outDir, "all_indels_DT40.isomut"))
indel$sample <- substr(indel$X.sample_name, 0, 5)

#Fischer scores separately for insertions and deletions
#first insertions
ins <- filter(indel, type == "INS")
cins <- matrix(0, ncol = samplenum, nrow = 201)
for (i in 1:201) \{ \
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

write.table(snv_filt, file = "all_SNVs_postprocessed_DT40.isomut", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(del_filt, file = "all_DELs_postprocessed_DT40.isomut", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ins_filt, file = "all_INSs_postprocessed_DT40.isomut", sep = "\t", row.names = FALSE, quote = FALSE)

# triplet spectra

snv_filt$chr <- gsub("MT", "M", snv_filt$chr)
snv_filt$chr <- factor(snv_filt$chr, levels = c(1:28, "W", "Z", "M"))
snv_range <- GRanges(paste0("chr", snv_filt$chr), ranges = IRanges(snv_filt$pos - 1, end = snv_filt$pos + 1))
snv_seq <- getSeq(Ggallus, snv_range)
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


snv_filt$sorting <- factor(snv_filt$sorting) 

snv_spectra <- as.data.frame.matrix(table(snv_filt$sample, snv_filt$sorting))

write.table(snv_spectra, "DT40_triplets.csv")