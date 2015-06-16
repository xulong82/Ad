library(ape)
library(amap)
library(dplyr)
library(xlsx)
library(Biobase)
library(biomaRt)
library(ggplot2)
library(ggvis)
library(VennDiagram)
library(pheatmap)
library(contrast)
library(igraph)

rm(list = ls())
setwd("~/Dropbox/GitHub/Load/")

myList <- list()
load("../X/summary.rdt")
summary <- function(x) table %>% filter(query %in% x)
source("../X/kegg.R")  # KEGG

load("data/ge.rdt")
ge0 <- sapply(ge, function(x) x$TPM) %>% as.data.frame
rownames(ge0) <- ge[[1]]$gene_id
ge0 <- ge0[apply(ge0, 1, function(x) max(x) > 10 & sum(x > 0) > 2), ]
colnames(ge0) <- gsub("Bin-1", "Bin1", colnames(ge0))
grp1 <- c("ApoE4", "Apoe", "Bin1", "Cd2ap", "Clu") 
grp2 <- c("B6", "ApoE4", "Apoe", "Bin1", "Cd2ap", "Clu") 

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
martId <- getBM(c("ensembl_gene_id", "external_gene_name"), "ensembl_gene_id", rownames(ge0), mart)
ge0 <- apply(ge0[martId[, "ensembl_gene_id"], ], 2, function(x) tapply(x, martId[, 2], sum)) %>% as.data.frame

spInf <- data.frame(sample = colnames(ge0), grp = factor(gsub("-.*", "", colnames(ge0)), levels = grp2))
graph1 <- function(x) spInf %>% mutate(value = c(as.matrix(ge0[x, ]))) %>%
  ggplot(aes(x = grp, y = value, label = sample)) + geom_boxplot(aes(fill = grp), width = 0.5) +
  theme_bw() + xlab("") + ylab("") + ggtitle(x) + geom_text(size = 2)
graph1("Bin1")

ge0$"Bin1-6962" <- NULL  # Remove 6962
ge <- as.matrix(ge0)  # downstream analysis
geneId0 <- rownames(ge)[apply(ge, 1, function(x) max(x) > 50)]
grp3 <- colnames(ge)
grp4 <- factor(gsub("-.*", "", grp3), levels = grp2)

pdf(file = "Figure/hist.pdf")
par(mfrow = c(3, 4))
lapply(1:35, function(x) hist(log2(ge[, x] + 1), main = sample[x], xlab = "", ylab = ""))
dev.off()

# GLM fit
fit <- apply(log2(ge + 1), 1, function (x) lm(x ~ grp4))
fit.pv <- sapply(fit, function(x) summary(x)$coefficients[-1, "Pr(>|t|)"]) %>% t
fit.et <- sapply(fit, function(x) summary(x)$coefficients[-1, "Estimate"]) %>% t
colnames(fit.et) <- colnames(fit.pv) <- grp1

# MT vs B6 
vsB6 <- list()
vsB6$up <- up <- sapply(colnames(fit.pv), function(x) which(fit.pv[, x] < 0.05 & fit.et[, x] > 0.1) %>% names)
vsB6$down <- down <- sapply(colnames(fit.pv), function(x) which(fit.pv[, x] < 0.05 & fit.et[, x] < -0.1) %>% names)
vsB6$all <- all <- sapply(colnames(fit.pv), function(x) which(fit.pv[, x] < 0.05 & abs(fit.et[, x]) > 0.1) %>% names)
vsB6$kegg$up <- lapply(up, kegg)
vsB6$kegg$down <- lapply(down, kegg)

myList$vsB6 <- vsB6  # Save

(up_term <- lapply(up_kegg, function(x) x$Term[x$Pvalue < 0.01]))
y = sort(table(unlist(up_term)), decreasing = T)
(y = y[y > 1])
(down_term <- lapply(down_kegg, function(x) x$Term[x$Pvalue < 0.01]))
y = sort(table(unlist(down_term)), decreasing = T)
(y = y[y > 1])

x <- rbind(up = sapply(up, length), down = sapply(down, length))
df <- data.frame(number = c(x), gene = rep(colnames(x), each = 2), group = rep(rownames(x), 5))
geneId1 <- unlist(all) %>% unique

profile1 <- lapply(1:5, function(x) apply(fit.pv, 1, function(y) sum(y < 0.05) == x)  %>% which %>% names)
geneId2 <- unlist(profile1)  # Interactive Venn Diagram???

pdf("Figure/number.pdf", height = 3, width = 6)
ggplot(df, aes(x = gene, y = number, fill = group)) + geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values= c("firebrick1", "dodgerblue3")) + theme_bw() + xlab("") + ylab("") 
dev.off()

# EXCEL
vsB6 <- sapply(grp1, function(x) data.frame(symbol = all[[x]], ge[all[[x]], grep(paste0("^(", x, "|B6)"), sample)]))
vsB6 <- sapply(grp1, function(x) { y = vsB6[[x]]; y$symbol = as.character(y$symbol);
  mutate(y, B6Avg = rowMeans(y[grep("B6", colnames(y))]), MtAvg = rowMeans(y[grep(x, colnames(y))]), 
         Pvalue = fit.pv[y$symbol, x], Log2FC = fit.et[y$symbol, x])})
lapply(grp1, function (x) write.xlsx(vsB6[[x]], file = "data/vsB6.xlsx", sheetName = x, append = T))
lapply(grp1, function (x) write.xlsx(up_kegg[[x]], file = "data/up_kegg.xlsx", sheetName = x, append = T))
lapply(grp1, function (x) write.xlsx(down_kegg[[x]], file = "data/down_kegg.xlsx", sheetName = x, append = T))

# ApoE4 vs Apoe 
fit.Apoe <- sapply(fit, function(x) contrast(x, list(grp4 = "ApoE4"), list(grp4 = "Apoe")))[c("Contrast", "Pvalue"), ] %>% t
Apoe <- list(up = fit.Apoe[apply(fit.Apoe, 1, function(x) x[1] > 0.1 & x[2] < 0.05), ],
             down = fit.Apoe[apply(fit.Apoe, 1, function(x) x[1] < -0.1 & x[2] < 0.05), ])
Apoe <- lapply(Apoe, function(x) cbind(ge[rownames(x), grep("^Apo", colnames(ge))],
  ApoE4Avg = rowMeans(ge[rownames(x), grep("ApoE4", sample)]), ApoeAvg = rowMeans(ge[rownames(x), grep("Apoe", sample)]), x))
Apoe <- lapply(Apoe, function(x) {colnames(x) = gsub("Contrast", "Log2FC", colnames(x)); x})
(Apoe$kegg$up <- kegg(rownames(Apoe$up)))
(Apoe$kegg$down <- kegg(rownames(Apoe$down)))

myList$Apoe <- Apoe  # Save

write.xlsx(Apoe$up, file = "data/Apoe.xlsx", sheetName = "UP", append = T)
write.xlsx(Apoe$down, file = "data/Apoe.xlsx", sheetName = "DOWN", append = T)
write.xlsx(Apoe$kegg$up, file = "data/Apoe.xlsx", sheetName = "UP_KEGG", append = T)
write.xlsx(Apoe$kegg$down, file = "data/Apoe.xlsx", sheetName = "DOWN_KEGG", append = T)

# Key regulators
lapply(grp1, function(x) write.xlsx(all[[x]], file = "Regulator/InRegulator.xlsx", sheetName = x, append = T, col.names = F))

iRegulon <- function(file) {  # PARSE THE IREGULON OUTPUT
  ireg <- read.delim(paste("Regulator", file, sep = "/"), comment.char = ";", stringsAsFactors = F)
  factor <- lapply(1:nrow(ireg), function(x) unlist(strsplit(ireg$Transcription.factor[x], split = ",")))
  target <- lapply(1:nrow(ireg), function(x) unlist(strsplit(ireg$Target.genes[x], split = ",")))
  factor <- lapply(factor, function(x) x[x %in% geneId0])  # only mild expression, post-hoc filter required
  idx <- sapply(factor, length) > 0; factor <- factor[idx]; target <- target[idx]
  edges <- lapply(1:length(factor), function(x) expand.grid(factor[[x]], target[[x]], stringsAsFactors = F))
  edges <- do.call(rbind, edges); edges <- edges[! duplicated(edges), ]  %>% return
}

file <- list.files(path = "./Regulator", pattern = "*.tsv")
regulators <- lapply(file, iRegulon)
names(regulators) <- gsub(".tsv", "", file)

myList$regulators <- regulators

lapply(grp1, function(x) write.xlsx(regulators[[x]], file = "Regulator/OutRegulator.xlsx", sheetName = x, append = T, row.names = F))
factors <- sapply(grp1, function(x) unique(regulators[[x]][, "Var1"]))  # Interactive Venn Diagram???
venn.diagram(factors, imagetype = "png", file = "Regulator/factorsVenn.png")

iGraph <- function (edges) {
  igraph.dt <- graph.data.frame(edges)
  igraph.dt$layout <- layout.sphere
  V(igraph.dt)$color = rep("chartreuse3", length(V(igraph.dt)$name))
  V(igraph.dt)$color[V(igraph.dt)$name %in% edges$Var1] <- "gold"
  V(igraph.dt)$size = (degree(igraph.dt) - 1) / (max(degree(igraph.dt)) - 1) * 10 + 2
  V(igraph.dt)$label.cex = (degree(igraph.dt) - 1) / (max(degree(igraph.dt)) - 1) * 2 + 0.8
  V(igraph.dt)$label.color = rep("dodgerblue3", length(V(igraph.dt)$name))
  V(igraph.dt)$label.color[V(igraph.dt)$color == "gold"] = "firebrick1"
  plot.igraph(igraph.dt, vertex.frame.color = "white", edge.arrow.size = 0.3)
}

pdf(file = "Regulator/igraph.pdf")
lapply(regulators, function(x) iGraph(x))
dev.off()

# Two-way HC
mycol <- rep("black", ncol(ge))
mycol[grep("Apoe" ,colnames(ge))] <- "firebrick1" 
mycol[grep("ApoE4" ,colnames(ge))] <- "gold1" 
mycol[grep("Bin1" ,colnames(ge))] <- "dodgerblue3" 
mycol[grep("Clu" ,colnames(ge))] <- "chartreuse3" 
mycol[grep("Cd2ap" ,colnames(ge))] <- "darkorchid2" 
hc1 <- hcluster(t(ge[geneId1, ]), method = "pearson", link = "average") %>% as.phylo
pdf("Figure/hc1.pdf", height = 5)
par(mfrow = c(1, 1)); plot(hc1, edge.width=2, font=2, cex=0.7, label.offset=1e-3, tip.color = mycol, direction="downward")
dev.off()
