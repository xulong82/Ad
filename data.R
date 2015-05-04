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

# ----------------------------------------------------------

trim <- read.delim("data/trim.txt", header = F, stringsAsFactors = F)
rsem <- read.delim("data/rsem_c3h.txt", header = F, stringsAsFactors = F) 

qc <- full_join(trim, rsem, by = "V1")[, -3] 
colnames(qc) <- c("sample", "trim", "count", "bowtie")

group <- c("WIN", "MIN", "WNONP", "MNONP", "WPLM", "MPLM") 
qc$sample <- gsub("_.*", "", qc$sample)
qc$group <- gsub("[123]", "", qc$sample) %>% factor(levels = group)
qc$genotype <- gsub("^(M|W).*", "\\1", qc$sample) %>% factor(levels = c("W", "M"))

qc$count <- qc$count * 1e-6
qc$trim <- as.numeric(gsub("%", "", qc$trim)) * 1e-2
qc$bowtie<- as.numeric(gsub("%", "", qc$bowtie)) * 1e-2

qc$aligned <- qc$count * qc$bowtie
qc$group1 <- as.numeric(qc$group)

ggvis_boxplots <- function(x) {x %>% add_axis("x", values = NULL, title = "Group") %>%
  layer_boxplots(fill=~group, width = 0.3) %>% layer_text(text:=~sample)}

qc %>% ggvis(~group1, ~trim) %>% scale_numeric("y", domain=c(0.5, 1)) %>% ggvis_boxplots()
qc %>% ggvis(~group1, ~aligned) %>% scale_numeric("y", domain=c(15, 35)) %>% ggvis_boxplots()

qc$f_count <- qc$count / qc$count[1]
qc$f_aligned <- qc$aligned / qc$aligned[1]

qc$f_spike1 <- colSums(spike) / colSums(spike)[1]
qc$f_spike2 <- c(W1NONP = 1, apply(spike[, -1], 2, function(x) summary(lm(x ~ -1 + spike[, 1]))$coefficients[1, 1]))

txInf1 <- txInf %>% filter(! grepl("ERCC", gene_id))
raw <- apply(genes, 2, function(x) tapply(x, txInf1$symbol, sum))
raw <- raw[apply(raw, 1, function (x) length(x[x > 20]) > 2), ]

total <- t(apply(raw, 1, "/", qc$f_count))
aligned <- t(apply(raw, 1, "/", qc$f_aligned))
spike2 <- t(apply(raw, 1, "/", qc$f_spike2))

geList <- list()
geList$raw = raw; geList$total = total; geList$aligned = aligned; geList$spike2 = spike2; 
save(geList, file = "Shiny/geList.rdt")

std_err <- apply(spike[, -1], 2, function(x) summary(lm(x ~ -1 + spike[, 1]))$coefficients[1, 2])
pval <- apply(spike[, -1], 2, function(x) summary(lm(x ~ -1 + spike[, 1]))$coefficients[1, 4])

qc %>% ggvis(~group1, ~f_spike1) %>% scale_numeric("y", domain=c(.5, 2.0)) %>% ggvis_boxplots()

ensId <- txInf %>% filter(symbol == "Actb") %>% select(gene_id)
qc$actb <- colSums(genes[ensId$gene_id, ])

qc <- mutate(qc, actb_f_count = actb / f_count, actb_f_aligned = actb / f_aligned,
             actb_f_spike1 = actb / f_spike1, actb_f_spike2 = actb / f_spike2)

qc %>% ggvis(~group1, ~actb_f_spike1) %>% ggvis_boxplots()

ensId <- txInf %>% filter(symbol == "Hspa2") %>% select(gene_id)
qc$hspa2 = c(t(genes[ensId$gene_id, ]))

qc <- mutate(qc, hspa2_f_count = hspa2 / f_count, hspa2_f_aligned = hspa2 / f_aligned,
             hspa2_f_spike1 = hspa2 / f_spike1, hspa2_f_spike2 = hspa2 / f_spike2)

qc %>% ggvis(~group1, ~hspa2) %>% ggvis_boxplots()

ensId <- txInf %>% filter(symbol == "Actb") %>% select(gene_id)
qc$actb2 <- colSums(genes2[ensId$gene_id, ])

qc <- mutate(qc, actb2_f_spike3 = actb2 / f_spike3)

qc %>% ggvis(~group1, ~actb2) %>% ggvis_boxplots()
qc %>% ggvis(~group1, ~actb2_f_spike3) %>% ggvis_boxplots()

summary(sort(apply(genes, 2, sum)))
summary(sort(apply(genes2, 2, sum)))

hc1 <- hcluster(t(genes), method = "pearson", link = "average") %>% as.phylo
par(mfrow = c(1, 1)); plot(hc1, edge.width=2, font=2, cex=0.7, label.offset=1e-3, direction="downward")

