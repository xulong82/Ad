# Copyright: Xulong Wang (xulong.wang@jax.org)
# Identify CO-EXPRESSION NETWORKS

library(ape)
library(amap)
library(ggplot2)
library(dynamicTreeCut)

rm(list = ls())
setwd("~/Dropbox/AD")
load("R/bc2014.rdt")
dt <- dt.bc

# --- DATA --- 
age <- factor(gsub("^.*(2m|4m|5m|6m).*", "\\1", colnames(dt)), levels = c("2m", "4m", "5m", "6m"))
group <- factor(gsub("^.*(WT|APP).*", "\\1", colnames(dt)), levels = c("WT", "APP"))
uid <- paste(age, group, sep = "_")
conditions <- c("2m_WT", "2m_APP", "4m_WT", "4m_APP", "5m_WT", "5m_APP", "6m_WT", "6m_APP")

dt <- dt[apply(dt, 1, sd) > 0.1, ]
aov.pval <- NULL  # --- ANOVA ---
for (i in 1:nrow(dt)) {
  if (i %% 1e3 == 0) cat(i, "\n")
  dat1 <- data.frame(y = as.matrix(dt)[i, ], age, group)
  aov <- aov(y ~ age * group, data = dat1)
  aov.pval <- c(aov.pval, min(summary(aov)[[1]][["Pr(>F)"]], na.rm = T))
}
dt <- dt[aov.pval < 0.05, ]

byConditions <- NULL
for (idx in conditions) byConditions <- cbind(byConditions, rowMeans(dt[, uid == idx]))
colnames(byConditions) <- conditions

dt <- as.data.frame(byConditions)
n.gene <- nrow(dt)
geneId <- rownames(dt)

hc1 <- hcluster(t(dt), method = "pearson", link = "average")
plot(as.phylo(hc1), type = "unrooted", cex = .5, font = 2, lab4ut = "axial")

# --- TOPOLOGY OVERLAP MATRIX --- 
sft <- NULL
breaks <- 20
similarity <- cor(t(dt), method = "pearson")
diag(similarity) <- 0
for (beta in seq(1, 31, 2)) {
  adjacency <- abs(similarity)^beta
  k1 <- rowSums(adjacency)
  discretized.k = cut(k1, breaks)
  dk = tapply(k1, discretized.k, mean)
  p.dk = tapply(k1, discretized.k, length) / length(k1)
  sum(p.dk, na.rm = T)
  
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
branch <- cutreeDynamic(dendro = dendro, distM = diss.tom, 
                        cutHeight = 0.97, minClusterSize = 30, deepSplit = 2,
                        pamRespectsDendro = F)

names(branch) <- geneId
table(branch)

eigene <- NULL
branchId <- names(table(branch))[-1]
for (idx in branchId) {
  dt1 <- dt[branch == idx, ]
  dt1 <- t(apply(dt1, 1, scale))
  colnames(dt1) <- conditions
  svd1 <- svd(dt1)
  eigene <- rbind(eigene, svd1$v[, 1])
}
colnames(eigene) <- conditions
rownames(eigene) <- branchId

mod.sim <- cor(t(eigene))
diag(mod.sim) <- 0
mod.sim[upper.tri(mod.sim)] <- 0
which(mod.sim > 0.9, arr.ind = T)

barplot(eigene[11, ], las = 2)
pdf("~/Desktop/eigengene.pdf")
par(mfrow = c(1, 1))
for (i in 1:length(branchId)) barplot(eigen.gene[i, ], las = 2)
dev.off()

# --- STATIC CUTREE ---
branch <- cutree(tree1, h = 0.99)
true.branch <- table(branch) >= min.size
branch[!true.branch[branch]] <- 0;
branchId <- sort(unique(branch))
map <- data.frame(row.names = branchId, idx = 1:length(branchId) - 1)
branch <- map[as.character(branch), "idx"]
# --- WGCNA ---
library(WGCNA)
net <- blockwiseModules(t(dt), power = 27, minModuleSize = 30, deepSplit = 2, 
         reassignThreshold = 0, mergeCutHeight = 0.15, 
         numericLabels = T, pamRespectsDendro = F, verbose = 3)
