# Copyright: Xulong Wang (xulong.wang@jax.org)
# Identify CO-EXPRESSION NETWORKS

library(ape)
library(amap)
library(ggplot2)
library(ggvis)
library(WGCNA)
library(dynamicTreeCut)

rm(list = ls())
setwd("~/Dropbox/GitHub/Ad")
load("./data/brain_bc2014.rdt")
load("./data/retina_bc2014.rdt")
source("function.R")

dt <- dt.bc
aov.pval <- NULL  
for (i in 1:nrow(dt)) {
  if (i %% 1e3 == 0) cat(i, "\n")
  dat1 <- data.frame(y = as.matrix(dt)[i, ], age, group)
  aov <- aov(y ~ age * group, data = dat1)
  aov.pval <- c(aov.pval, min(summary(aov)[[1]][["Pr(>F)"]], na.rm = T))
}
dt.cond <- NULL
dt <- dt[aov.pval < 0.05, ]
conditions <- conditions[! grepl("4m", conditions)]  # retina only
for (idx in conditions) dt.cond <- cbind(dt.cond, rowMeans(dt[, uid == idx]))
colnames(dt.cond) <- conditions
dt.cond <- dt.cond[apply(dt.cond, 1, sd) > 0.1, ]
dt.cond <- dt.cond[apply(dt.cond, 1, sd) > 0.2, ]  # retina only
dt <- as.data.frame(dt.cond)
n.gene <- nrow(dt)
geneId <- rownames(dt)
str(dt.cond)

sft <- NULL  # --- Scale free topology
similarity <- cor(t(dt), method = "pearson")
diag(similarity) <- 0
for (beta in seq(1, 31, 2)) {
  adjacency <- abs(similarity)^beta
  k1 <- rowSums(adjacency)
  discretized.k = cut(k1, 20)  # breaks:20
  dk = tapply(k1, discretized.k, mean)
  p.dk = tapply(k1, discretized.k, length) / length(k1)
  dk = dk[! is.na(dk)]
  p.dk = p.dk[! is.na(p.dk)]
  fit <- summary(lm(log10(p.dk) ~ log10(dk)))
  sr2 <- -sign(fit$coefficients["log10(dk)", "Estimate"]) * fit$r.squared
  sft <- rbind(sft, c(beta, sr2))
}
plot(sft, main = "Scale free topology", xlab = "beta", ylab = "signed R2")
sft.beta <- 27
adjacency <- abs(similarity)^sft.beta
tom <- matrix(0, n.gene, n.gene)
rownames(tom) <- colnames(tom) <- geneId
k1 <- rowSums(adjacency)
for (i in 2:n.gene) {
  if (i %% 1e2 == 0) cat(i, "\n")
  for (j in 1:(i-1)) {
    kij <- min(k1[c(i, j)])
    lij <- sum(adjacency[i, ] * adjacency[j, ])
    wij <- (lij + adjacency[i, j]) / (kij + 1 - adjacency[i, j])
    tom[i, j] <- wij
  }
}
diag(tom) <- 1
tom[upper.tri(tom)] <- t(tom)[upper.tri(tom)]
diss.tom <- 1 - tom

# --- NETWORK ---
dendro <- hclust(as.dist(diss.tom), method = "average")
branch <- cutreeDynamic(dendro = dendro, distM = diss.tom, cutHeight = 0.99, minClusterSize = 30)
branch <- cutreeDynamic(dendro = dendro, distM = diss.tom, cutHeight = 0.999, minClusterSize = 30)  # retina only
names(branch) <- geneId
table(branch)
branch["App"]

mod <- length(table(branch)) - 1
geneId.network <- gk.network <- list()
for (idx in 1:mod) {
  gene1 <- names(branch[branch == idx])
  geneId.network[[idx]] <- gene1
  gk.network[[idx]] <- myGK(gene1)
}

eigene <- NULL
branchId <- names(table(branch))[-1]
for (idx in branchId) {
  dt1 <- dt[branch == idx, ]
  dt1 <- t(apply(dt1, 1, scale))
  svd1 <- svd(dt1)
  eigene <- rbind(eigene, svd1$v[, 1])
}
dimnames(eigene) <- list(branchId, conditions)

module <- paste("Module", branchId, sep = "")
group <- factor(conditions, levels = conditions)
gdt <- data.frame(value = c(eigene), module = rep(module, 6), group = rep(group, each = mod))
gdt$geno <- factor(gsub("^.*(WT|APP).*", "\\1", gdt$group), levels = c("WT", "APP"))
ggplot(gdt, aes(x = group, y = value, fill = geno)) + 
  geom_bar(stat = "identity") + facet_grid(. ~ module) +
  scale_fill_manual(values = c("red", "blue")) +
  theme_bw() + xlab("") + ylab("") + 
  theme(axis.text.x = element_text(angle = 90))

network <- list()
network$gdt <- gdt
network$symbol <- geneId.network
network$gk <- gk.network
save(network, file = "data/network_brain.rdt")
save(network, file = "data/network_retina.rdt")

x.dt = as.matrix(dt[branch == 6, ])
matplot(t(scale(x.dt)), type = "l")
x.dt <- t(apply(x.dt, 1, scale))
gdt <- data.frame(value = c(x.dt), gene = rep(rownames(x.dt), 8), group = rep(group, each = nrow(x.dt)))
gdt %>% ggvis(~group, ~value, stroke = ~gene) %>% layer_lines() %>% hide_legend("stroke") %>%
  add_tooltip(function (x) paste("Gene: ", x$gene), "hover")

tom.mod <- tom[branch == 6, branch == 6]
visant <- exportNetworkToVisANT(tom.mod, file = "./data/visant.txt", weighted = T, threshold = 0)
cyt <- exportNetworkToCytoscape(tom.mod, edgeFile = "./data/cytodge.txt", nodeFile = "./data/cytonode.txt", 
                               weighted = TRUE, threshold = 0.02)
summary.mod <- table[match(geneId[branch == 6], table$query), ]
