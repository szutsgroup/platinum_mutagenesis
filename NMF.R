library(MutationalPatterns)
library(BSgenome.Ggallus.UCSC.galGal4)

# de novo NMF with DT40 and TK6 samples, aggregated by drugs and doses

filter(snv_filt, treatment != "starting_clone") -> snv2

vcf_list <- vector(mode = "list", length = length(unique(snv2$sample)))
names(vcf_list) <- unique(snv2$sample)

for (i in unique(snv2$sample)) {
  t <- filter(snv2, sample == i)
  
  vcf_list[[i]] <- with(t, GRanges(seqnames = paste0("chr", chr), ranges = IRanges(start = pos, end = pos), seqinfo = seqinfo(Ggallus)))
  names(vcf_list[[i]] ) <- paste0("chr", t$chr, ":", t$pos, "_", t$ref, "/", t$mut)
  mcols(vcf_list[[i]] )$paramRangeID <- factor(NA)
  mcols(vcf_list[[i]] )$REF <- DNAStringSet(t$ref)
  
  mcols(vcf_list[[i]] )$ALT <- DNAStringSetList(unname(sapply(t$mut, DNAStringSet)))
  mcols(vcf_list[[i]] )$QUAL <- as.numeric(NA)
  mcols(vcf_list[[i]])$FILTER <- "PASS"
}

mut_mat <- mut_matrix(vcf_list = vcf_list, ref_genome = "Ggallus")
mut_mat <- mut_mat[,c(treat_list$mock, treat_list$cisplatin_7.5uM, treat_list$cisplatin_10uM, treat_list$carboplatin_75uM, treat_list$carboplatin_100uM, treat_list$oxaliplatin_2.5uM, treat_list$oxaliplatin_4.5uM)]

tk6 = read.delim("TK6_triplets.csv", sep = " ")
colnames(tk6) <- rownames(mut_mat)
tk6 <- tk6[c("TK11", "TK12", "TK13", paste0("TK", 45:53)),]

mut_mat <- cbind.data.frame(mut_mat, t(tk6))

mut_mat <- mut_mat + 0.0001

nmf_res_2 <- extract_signatures(mut_mat, rank = 2, nrun = 100)
colnames(nmf_res_2$signatures) <- c("Signature A", "Signature B")
rownames(nmf_res_2$contribution) <- c("Signature A", "Signature B")
contribution <-  apply(nmf_res_2$contribution, 2, function(x) x/sum(x))
contribution <- sweep(contribution, 2, c(table(snv_filt$sample)[setdiff(unlist(treat_list), treat_list$starting_clone)], rowSums(tk6)), "*")

as.data.frame(contribution) %>% mutate(Sig = rownames(.)) %>% gather(Sample, Estimated, -Sig) -> con_molten
con_molten$Treatment <- NA
for (i in seq_along(treat_list)) {
  con_molten[which(con_molten$Sample %in% treat_list[[i]]), "Treatment"] <- names(treat_list)[i]
}
con_molten[which(con_molten$Sample %in% paste0("TK", 11:13)), "Treatment"] <- "TK6_mock"
con_molten[which(con_molten$Sample %in% paste0("TK", c(45, 48, 49))), "Treatment"] <- "TK6_cisplatin"
con_molten[which(con_molten$Sample %in% paste0("TK", c(46, 50, 51))), "Treatment"] <- "TK6_carboplatin"
con_molten[which(con_molten$Sample %in% paste0("TK", c(47, 52, 53))), "Treatment"] <- "TK6_oxaliplatin"


mutate(con_molten, Sig = factor(Sig, levels = c("Signature B", "Signature A")), Treatment = factor(Treatment, levels = c(names(treat_list)[-1], "TK6_mock", "TK6_cisplatin", "TK6_carboplatin", "TK6_oxaliplatin"))) %>% group_by(Treatment, Sig) %>% summarize(Mean = mean(Estimated), SEM = sd(Estimated)/sqrt(n())) %>% .["Mean"] %>% unlist() -> means
means[which(seq_along(means) %% 2 == 1)] <- means[which(seq_along(means) %% 2 == 1)] + means[which(seq_along(means) %% 2 == 0)]


mutate(con_molten, Sig = factor(Sig, levels = c("Signature B", "Signature A")), Treatment = factor(Treatment, levels = c(names(treat_list)[-1], "TK6_mock", "TK6_cisplatin", "TK6_carboplatin", "TK6_oxaliplatin"))) %>% group_by(Treatment, Sig) %>% summarize(Mean = mean(Estimated), SEM = sd(Estimated)/sqrt(n())) %>% ggplot(aes(x = Treatment, fill = Sig, y = Mean)) + geom_col() + geom_errorbar(aes(ymin = means - SEM, ymax = means + SEM), width = .2, size = 1, position = position_dodge()) + cowplot::theme_cowplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


# error

rmsd <- function(a, b) {
  ret <- sqrt(mean(unlist((a-b)^2)))
  return(ret)
}


error_df <- data.frame(RMSD = 0, rel_RMSD = 0, Treatment = c(names(treat_list)[-1], "TK6_mock", "TK6_cisplatin", "TK6_carboplatin", "TK6_oxaliplatin"))

TK6_mock <- paste0("TK", 11:13)
TK6_cisplatin <- paste0("TK", c(45, 48, 49))
TK6_carboplatin <- paste0("TK", c(46, 50, 51))
TK6_oxaliplatin <- paste0("TK", c(47, 52, 53))

for (i in 1:nrow(error_df)) {
  tt <- as.character(error_df[i, "Treatment"])
  if (substr(tt, 1, 2) != "TK") {  
    error_df[i, "RMSD"] <- rmsd(mut_mat[,treat_list[[tt]]], nmf_res_2$reconstructed[,treat_list[[tt]]])
    error_df[i, "rel_RMSD"] <- error_df[i, "RMSD"] / sum(mut_mat[,treat_list[[tt]]])
  } else {
    error_df[i, "RMSD"] <- rmsd(mut_mat[,get(tt)], nmf_res_2$reconstructed[,get(tt)])
    error_df[i, "rel_RMSD"] <- error_df[i, "RMSD"] / sum(mut_mat[,get(tt)])
  }
}

mutate(error_df, Treatment = factor(Treatment, levels = c(names(treat_list)[-1], "TK6_mock", "TK6_cisplatin", "TK6_carboplatin", "TK6_oxaliplatin"))) %>% ggplot(aes(x = Treatment, y = rel_RMSD)) + geom_col(fill = "steelblue") + cowplot::theme_cowplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + xlab("") + ylab("Relative reconstruction error")


# components

what <- apply(nmf_res_2$signatures, 2, function(x) x/sum(x))

spaces <- rep(c(2, rep(0.4, 15)), 6)

colors <- c(rgb(207/255, 19/255, 129/255), 
            rgb(111/255, 53/255, 132/255),
            rgb(90/255, 200/255, 209/255), 
            rgb(44/255, 53/255, 138/255), 
            rgb(20/255, 138/255, 76/255), 
            rgb(253/255, 209/255, 45/255))

trips <- rownames(what)

for (i in 1:ncol(what)) {
  pdf(paste0("NMF_", gsub(" ", "_", colnames(what)[i]), ".pdf"), width = 20, height = 7)
  opar <- par()
  par(mar = c(5.1, 7.1, 4.1, 2.1))
  yscale <- c(0, max(as.matrix(what)) * 1.2)
  barplot(as.numeric(what[,i]), 
          main = colnames(what[i]), 
          col = rep(colors, each = 16), 
          las = 2, 
          names.arg = trips, 
          cex.names = 1.3, 
          space = spaces, 
          cex.main = 2, 
          ylim = yscale, 
          border = NA, 
          cex.axis = 1.6, 
          family = "mono", 
          lwd = 8)
  for (j in 1:6) {
    rect(xleft = 2 + (j-1) * (2 + (16*1 + 15*0.4)), ybottom = (par("usr")[4] - par("usr")[3]) * 0.97, xright = j * (15*1.4 + 2 + 1), ytop = par("usr")[4], col = colors[j], border = NA)
    text(label = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")[j], x = (j - 1) * (15*1.4 + 2) + j + 11.5, y = (par("usr")[4] - par("usr")[3]) * 0.93, cex = 2, font = 2, family = "mono")
  }
  dev.off()
}



