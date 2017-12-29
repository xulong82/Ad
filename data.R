library(dplyr)
library(xlsx)
library(biomaRt)

# ----------------------------------------------------------

setwd("/data/xwang/Load")  # GE @ CADILLAC 
name <- list.files(path = "./RSEM/", pattern = "*.genes.results")
ge <- lapply(name, function(x) {
  cat(x, "\n"); filepath <- file.path("./RSEM/", x)
  read.delim(filepath, stringsAsFactors = F)
}); names(ge) <- gsub("-GES.*", "", name)
save(ge, file = "~/Dropbox/GitHub/Load/data/ge.rdt")

# ----------------------------------------------------------

rm(list = ls())
setwd("~/Dropbox/GitHub/Load/")
load("data/ge.rdt")

# ----------------------------------------------------------

ge_tpm <- sapply(ge, function(x) x$TPM) %>% as.data.frame
rownames(ge_tpm) <- ge[[1]]$gene_id

bin1 <- ge_tpm[grepl("^B", colnames(ge_tpm))]

bin1 <- bin1[apply(bin1, 1, function(x) max(x) > 10 & sum(x > 0) > 2), ]
colnames(bin1) <- gsub("Bin-1", "Bin1", colnames(bin1))
group <- factor(gsub("^(B6|Bin1).*", "\\1", colnames(bin1)), levels = c("B6", "Bin1"))

fit <- apply(log2(bin1 + 1), 1, function (x) lm(x ~ group))

fit.r2 <- sapply(fit, function (x) summary(x)$r.squared)
fit.fs <- sapply(fit, function (x) summary(x)$fstatistic)
fit.pv <- apply(fit.fs, 2, function (x) pf(x[1], x[2], x[3], lower.tail = F)) 
fit.qv <- p.adjust(fit.pv, method = "fdr")

sig_fit <- fit[fit.pv < 0.05 & fit.r2 > 0.3]
sig_glm <- cbind(bin1, fit.pv)[names(sig_fit), ]

sig_glm <- sig_glm %>% mutate(log2fc = log2(Bin1_avg) - log2(B6_avg))

ensId <- names(sig_fit)
mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
biomart <- getBM(c("ensembl_gene_id", "external_gene_name"), "ensembl_gene_id", ensId, mart)

sig_glm$symbol <- biomart$external_gene_name[match(ensId, biomart$ensembl_gene_id)]
sig_glm$B6_avg <- rowMeans(sig_glm[grepl("B6-", colnames(sig_glm))])
sig_glm$Bin1_avg <- rowMeans(sig_glm[grepl("Bin1-", colnames(sig_glm))])

write.xlsx(sig_glm, file = "Bin1.xlsx", sheetName = "sig_glm2", append = T)

source("../X/function.R")

bin1_up <- sig_glm$symbol[sig_glm$Bin1_avg > sig_glm$B6_avg]
bin1_dn <- sig_glm$symbol[sig_glm$Bin1_avg < sig_glm$B6_avg]

bin1_up_gk <- myGK(bin1_up[! is.na(bin1_up)])
bin1_dn_gk <- myGK(bin1_dn[! is.na(bin1_dn)])

shiny <- list()
shiny$sig <- sig_glm
shiny$sig_up$gk <-bin1_up_gk
shiny$sig_dn$gk <-bin1_dn_gk

load("../X/summary.rdt")
summary <- table %>% filter(query %in% sig_glm$symbol)

shiny$summary <- summary

save(shiny, file = "shiny/shiny.rdt")

write.xlsx(bin1_up_gk$KEGG, file = "Bin1.xlsx", sheetName = "sig_up_KEGG", append = T)
write.xlsx(bin1_dn_gk$KEGG, file = "Bin1.xlsx", sheetName = "sig_dn_KEGG", append = T)

write.xlsx(bin1_up_gk$GO$BP, file = "Bin1.xlsx", sheetName = "sig_up_GO_BP", append = T)

# ----------------------------------------------------------

library(igraph)

# PARSE THE IREGULON OUTPUT

# *** OPTIONAL: N
file <- "Regulator/up_iregulon.tsv"
ensId <- rownames(ge_tpm)[rowMeans(ge_tpm[grepl("Bin1-", colnames(sig_glm))]) > 1e2]

file <- "Regulator/dn_iregulon.tsv"
ensId <- rownames(ge_tpm)[rowMeans(ge_tpm[grepl("B6-", colnames(sig_glm))]) > 1e2]

mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
biomart <- getBM(c("ensembl_gene_id", "external_gene_name"), "ensembl_gene_id", ensId, mart)
univ <- biomart$external_gene_name

ireg <- read.delim(file, comment.char = ";", stringsAsFactors = F)
factor <- lapply(1:nrow(ireg), function(x) unlist(strsplit(ireg$Transcription.factor[x], split = ",")))
target <- lapply(1:nrow(ireg), function(x) unlist(strsplit(ireg$Target.genes[x], split = ",")))

factor <- lapply(factor, function(x) x[x %in% univ])
idx <- sapply(factor, length) > 0

factor <- factor[idx]
target <- target[idx]

edges <- lapply(1:length(factor), function(x) expand.grid(factor[[x]], target[[x]], stringsAsFactors = F))
edges <- do.call(rbind, edges)
edges <- edges[! duplicated(edges), ]

# VISUALIZATION: IGRAPH
igraph.dt <- graph.data.frame(edges)

igraph.dt$layout <- layout.sphere
igraph.dt$layout <- layout.circle
igraph.dt$layout <- layout.fruchterman.reingold 

V(igraph.dt)$color = rep("chartreuse3", length(V(igraph.dt)$name))
V(igraph.dt)$color[V(igraph.dt)$name %in% unlist(factor)] <- "gold"
V(igraph.dt)$size = (degree(igraph.dt) - 1) / (max(degree(igraph.dt)) - 1) * 10 + 2
V(igraph.dt)$label.cex = (degree(igraph.dt) - 1) / (max(degree(igraph.dt)) - 1) * 2 + 0.5
V(igraph.dt)$label.color = rep("dodgerblue3", length(V(igraph.dt)$name))
V(igraph.dt)$label.color[V(igraph.dt)$color == "gold"] = "firebrick1"

plot.igraph(igraph.dt, vertex.frame.color = "white", edge.arrow.size = 0.3)

igraphList <- list()
igraphList$NN <- igraph.dt

# ----------------------------------------------------------

library(DESeq2)

ge_count <- sapply(ge, function(x) x$expected_count) %>% as.data.frame
rownames(ge_count) <- ge[[1]]$gene_id

bin1 <- ge_tpm[grepl("^B", colnames(ge_tpm))]
bin1 <- bin1[apply(bin1, 1, function(x) max(x) > 10 & sum(x > 0) > 2), ] %>% round %>% as.matrix

colnames(bin1) <- gsub("Bin-1", "Bin1", colnames(bin1))
group <- factor(gsub("^(B6|Bin1).*", "\\1", colnames(bin1)), levels = c("B6", "Bin1"))

colData <- data.frame(row.names = colnames(bin1), condition = group)
dds <- DESeqDataSetFromMatrix(countData = bin1, colData = colData, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
res
summary(res)
plotMA(res, main = "DESeq2", ylim = c(-2,2))

resSig <- subset(res, padj < 0.1) %>% as.data.frame
resSig <- subset(res, pvalue < 0.05) %>% as.data.frame

hist(res@listData$padj) 
table(res@listData$padj < 0.1) 
table(res@listData$pval < 0.05) 
