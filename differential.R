# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: AD project - pairwise comparison
# Rev: Dec 30, 2014

library(EBSeq)
library(VennDiagram)
library(gplots)
library(ggplot2)

rm(list = ls())

setwd("~/Dropbox/AD")
load(file = "./R/brain2014.rdt")
load(file = "./R/cluster1.rdt")

de.month <- list()
for (month in c("4m", "5m", "6m")) {
  cat(month, "\n")
  
  tpm <- dt[, grep(month, colnames(dt))]
  tpm.wt <- tpm[, grep("WT", colnames(tpm))]
  tpm.app <- tpm[, grep("APP", colnames(tpm))]
  fc <- rowMeans(tpm.app) - rowMeans(tpm.wt) 
  
  treat <- gsub("^.*(WT|APP).*", "\\1", colnames(tpm))
  
  tt.pval <- apply(tpm, 1, function(x) pairwise.t.test(x, treat, p.adj = "none")$p.value)
  tt.qval <- p.adjust(tt.pval, method = "fdr")

  idx <- which(tt.pval < 0.05)
  name <- paste("de", month, sep = "")
  de.month[[name]] <- cbind(FC = fc, PVAL = tt.pval, QVAL = tt.qval)[idx, ]
}  # lm() and pairwise.t.test() do same


write.table(de.month[[1]], file = "./DE/4m.txt", sep = "\t", row.names = T, col.names = F, quote = F)
write.table(de.month[[2]], file = "./DE/5m.txt", sep = "\t", row.names = T, col.names = F, quote = F)
write.table(de.month[[3]], file = "./DE/6m.txt", sep = "\t", row.names = T, col.names = F, quote = F)

venn.diagram(list(de_4m = rownames(de.month$de4m), de_5m = rownames(de.month$de5m), de_6m = rownames(de.month$de6m)),
  fill = c("dodgerblue3", "firebrick1", "chartreuse3"), col = "transparent", alpha = .5, 
  cat.col = c("dodgerblue3", "firebrick1", "chartreuse3"), cat.cex = 2, cat.fontface = "bold",
  cex = 2, fontface = "bold", lwd = 0, file = "~/Dropbox/AD/Graphs/de_venn.tiff")

graph.dt <- rbind(de_4m = table(de.month$de4m[, "FC"] > 0),
                  de_5m = table(de.month$de5m[, "FC"] > 0),
                  de_6m = table(de.month$de6m[, "FC"] > 0))

pdf(file = "~/Dropbox/AD/Graphs/de_number.pdf", width = 4)
col.manual <- c("dodgerblue3", "firebrick1")
bar <- barplot(t(graph.dt), col = col.manual, border = NA, xlim = c(0, 10), ylim = c(0, 200), axes = F, font = 2, beside = T)
abline(0, 0, lwd = 5, col = "grey20")
text(x = bar, y = t(graph.dt) + 5, labels = t(graph.dt), font = 2, cex = 0.7, col = "grey20")
dev.off()

geneList <- NULL
for (month in c("de4m", "de5m", "de6m")) {
  geneList <- unique(c(geneList, rownames(de.month[[month]])))
  write.table(rownames(de.month[[month]]), row.names = FALSE, col.names = FALSE, quote = FALSE, 
              file = paste(paste("~/Dropbox/AD/DAVID", month, sep = "/"), "txt", sep = "."))
}

dt <- dt[geneList, ]
dt.wt<- dt[, grep("WT", colnames(dt))]
dt.app <- dt[, grep("APP", colnames(dt))]

hc1 <- hcluster(t(dt), method = "pearson", link = "average")
hc2 <- hcluster(t(dt.wt), method = "pearson", link = "average")
hc3 <- hcluster(t(dt.app), method = "pearson", link = "average")

pdf("./Graphs/phylo_sample.pdf", fonts = "Helvetica", width = 12)
par(mar = c(1, 4, 1, 2), mfrow = c(3, 1))
plot(as.phylo(hc1), edge.width = 2, font = 2, cex = 1, label.offset = 1e-2, direction = "downward")
plot(as.phylo(hc2), edge.width = 2, font = 2, cex = 1, label.offset = 1e-2, direction = "downward")
plot(as.phylo(hc3), edge.width = 2, font = 2, cex = 1, label.offset = 1e-2, direction = "downward")
dev.off()

# --- two month samples ---
gene.2m <- gene[["2m"]]
data2m.wt <- tpm2[gene.2m, grep("WT2m", colnames(tpm2))]
data2m.app <- tpm2[gene.2m, grep("APP2m", colnames(tpm2))]
mean.wt <- 2^(rowMeans(data2m.wt))
mean.app <- 2^(rowMeans(data2m.app))
log2.fc <- mean.app / mean.wt

howell <- data.frame(row.names = gene.2m, mean.wt, mean.app, log2.fc)
write.table(howell, file = "~/Dropbox/AD/2M/2m.tpm.txt", sep = "\t", row.names = T, col.names = T, quote = F)

# ---
intersect(gene.1, gene.2)
setdiff(gene.1, gene.2)  # only in gene.1
setdiff(gene.2, gene.1)  # only in gene.2

pdf(file = "~/Dropbox/AD/Figures/text2.pdf", width = 12)
textplot(matrix(c(intersect(gene.1, gene.2), rep(NA, 3)), nrow = 6), show.rownames = F, show.colnames = F)
dev.off()

# ---
tpm1 <- brain.tpm2
id1 <- which(apply(tpm1, 1, function(x) {min(x) < 0}))  # negative TPM from correction
tpm1[id1, ] <- tpm1[id1, ] + abs(apply(tpm1[id1, ], 1, min))
tpm2 <- log2(tpm1 + 1)  # log2 transformation

gene = list()
for (month in c("2m", "4m", "5m", "6m")) {  # EBSeq
  cat(month, "\n")
  tpm3 <- tpm1[, grep(month, colnames(tpm1))]
  group <- gsub("m.*", "m", colnames(tpm3))
  treat <- gsub("[2456]m", "", group)
  sizes = MedianNorm(tpm3)
  conditions <- as.factor(treat)
  ebout = EBTest(Data = tpm3, Conditions = conditions, sizeFactors = sizes, maxround = 5)
  ebpp = GetPPMat(ebout)
  ebfc = PostFC(ebout)
  gene[[month]] = which(ebpp[, "PPDE"] >= .75)
}  # EBSeq returns less
