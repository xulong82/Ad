# Copyright: Xulong Wang (xulong.wang@jax.org)
# PAIRWISE COMPARISON

library(gplots)
library(ggplot2)
library(VennDiagram)

rm(list = ls())
setwd("~/Dropbox/GitHub/Ad/data")
load("bc2014.rdt")

dt <- dt.bc
wt2m <- rowMeans(dt[, uid == "2m_WT"])
app2m <- rowMeans(dt[, uid == "2m_APP"])
wt456m <- dt[, uid %in% c("4m_WT", "5m_WT", "6m_WT")]
app456m <- dt[, uid %in% c("4m_APP", "5m_APP", "6m_APP")]
wt456m <- wt456m - wt2m
app456m <- app456m - app2m

dt <- cbind(wt456m, app456m)
hc1 <- hcluster(t(dt), method = "pearson", link = "average")
dev.off()
plot(as.phylo(hc1), type = "unrooted", cex = .5, font = 2, lab4ut = "axial")

dt <- dt.bc
pool <- c("2m", "4m", "5m", "6m")

dt <- cbind(wt456m, app456m)
pool <- c("4m", "5m", "6m")

dt <- dt[apply(dt, 1, sd) > 0.1, ]
de.month <- list()
for (month in pool) {
  cat(month, "\n")
  dt1 <- dt[, grep(month, colnames(dt))]
  dt1.wt <- dt1[, grep("WT", colnames(dt1))]
  dt1.app <- dt1[, grep("APP", colnames(dt1))]
  fc <- rowMeans(dt1.app) - rowMeans(dt1.wt) 
  treat <- gsub("^.*(WT|APP).*", "\\1", colnames(dt1))
  tt.pval <- apply(dt1, 1, function(x) pairwise.t.test(x, treat, p.adj = "none")$p.value)
  tt.qval <- p.adjust(tt.pval, method = "fdr")
  idx <- which(tt.qval < 0.05)
  name <- paste("DE", month, sep = "")
  de.month[[name]] <- cbind(FC = fc, PVAL = tt.pval, QVAL = tt.qval)[idx, ]
}  # lm() and pairwise.t.test() do same

lapply(de.month, write, "../TXT/de456.txt", append = T, ncolumns = 1e3)
lapply(de.month, write, "../TXT/de2456.txt", append = T, ncolumns = 1e3)

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

dev.off()
par(mar = c(1, 4, 1, 2), mfrow = c(3, 1))
plot(as.phylo(hc1), edge.width = 2, font = 2, cex = 1, label.offset = 1e-2, direction = "downward")
plot(as.phylo(hc2), edge.width = 2, font = 2, cex = 1, label.offset = 1e-2, direction = "downward")
plot(as.phylo(hc3), edge.width = 2, font = 2, cex = 1, label.offset = 1e-2, direction = "downward")

intersect(gene.1, gene.2)
setdiff(gene.1, gene.2)  # only in gene.1
setdiff(gene.2, gene.1)  # only in gene.2

textplot(matrix(c(intersect(gene.1, gene.2), rep(NA, 3)), nrow = 6), show.rownames = F, show.colnames = F)

# # --- EBSeq
# tpm1 <- brain.tpm2
# id1 <- which(apply(tpm1, 1, function(x) {min(x) < 0}))  # negative TPM from correction
# tpm1[id1, ] <- tpm1[id1, ] + abs(apply(tpm1[id1, ], 1, min))
# tpm2 <- log2(tpm1 + 1)  # log2 transformation
# 
# gene = list()
# for (month in c("2m", "4m", "5m", "6m")) {  # EBSeq
#   cat(month, "\n")
#   tpm3 <- tpm1[, grep(month, colnames(tpm1))]
#   group <- gsub("m.*", "m", colnames(tpm3))
#   treat <- gsub("[2456]m", "", group)
#   sizes = MedianNorm(tpm3)
#   conditions <- as.factor(treat)
#   ebout = EBTest(Data = tpm3, Conditions = conditions, sizeFactors = sizes, maxround = 5)
#   ebpp = GetPPMat(ebout)
#   ebfc = PostFC(ebout)
#   gene[[month]] = which(ebpp[, "PPDE"] >= .75)
# }  # EBSeq returns less
