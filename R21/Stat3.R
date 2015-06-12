library(dplyr)
library(contrast)
library(multcomp)
library(biomaRt)
library(ggplot2)
library(ggvis)
library(KEGGREST)
library(xlsx)
library(png)

rm(list = ls())
setwd("~/Dropbox/GitHub/Ad")

load("./data/brain2014.rdt"); glm.dt <- brain.tpm
# load("./data/retina2014.rdt"); glm.dt <- retina.tpm
cutoff <- quantile(c(as.matrix(glm.dt)), 0.25)  # TPM level
glm.dt <- glm.dt[apply(glm.dt, 1, function(x) max(x) > cutoff & sum(x > 0) > round(ncol(glm.dt) / 10)), ]
sample <- colnames(glm.dt)
age <- factor(gsub("^.*(2m|4m|5m|6m).*", "\\1", sample), levels = c("2m", "4m", "5m", "6m"))
geno <- factor(gsub("^.*(WT|APP).*", "\\1", sample), levels = c("WT", "APP"))
batch <- factor(gsub("^.*(2014|mouse).*", "\\1", sample), levels = c("2014", "mouse"))
uid <- paste(age, geno, sep = "_")
grp <- intersect(c("2m_WT", "2m_APP", "4m_WT", "4m_APP", "5m_WT", "5m_APP", "6m_WT", "6m_APP"), uid)
nid <- factor(uid, levels = grp)

ggvis1 <- function(x)  # Visualize the Stat expression
  data.frame(sample, uid, age, geno, batch, nid) %>% mutate(value = c(as.matrix(glm.dt[x, ]))) %>%
  ggvis(~as.numeric(nid), ~value) %>% layer_boxplots(fill=~nid, width = 0.5)

ggvis1("Stat1")
ggvis1("Stat3")

# Model the data with GLM 
glm.fit <- apply(glm.dt, 1, function (x) summary(lm(x ~ age + geno + batch + age * geno)))
fit.batch <- sapply(glm.fit, function (x) x$coefficients["batchmouse", "Estimate"])
glm.dt_bc <- glm.dt - as.matrix(fit.batch) %*% (as.numeric(batch) - 1)
glm.fit <- apply(glm.dt_bc, 1, function (x) lm(x ~ age + geno + age * geno))
fit.r2 <- sapply(glm.fit, function (x) summary(x)$r.squared)
fit.fs <- sapply(glm.fit, function (x) summary(x)$fstatistic)
fit.pval <- apply(fit.fs, 2, function (x) pf(x[1], x[2], x[3], lower.tail = F)) 
fit.qval <- p.adjust(fit.pv, method = "fdr")

glm.fit <- glm.fit[fit.pval < 0.05 & fit.r2 > 0.3]
fit.eff <- do.call(rbind, lapply(glm.fit, function (x) summary(x)$coefficients[, "Estimate"]))
fit.pval <- do.call(rbind, lapply(glm.fit, function (x) summary(x)$coefficients[, "Pr(>|t|)"]))

lapply(glm.fit[grep("Stat", names(glm.fit))], summary)

# Identify signal: ANOVA
aov.fit <- lapply(glm.fit, anova)
aov.fit[grep("Stat", names(aov.fit))]
aovGene <- lapply(c("age", "geno", "age:geno"), function(x1) {
  y1 = sapply(aov.fit, function(x2) x2[x1, "Pr(>F)"]); y1[y1 < 0.05]
}); names(aovGene) <- c("age", "app", "age:app")

# Identify signal: GLM estiamtes
binary <- apply(fit.pval, 2, function (x) x < 0.05) & apply(fit.eff, 2, function (x) abs(x) > 0.1)
binary <- binary[apply(binary[, -1], 1, any), -1]
signal <- apply(binary, 1, function (x) paste(x, collapse = "-"))
signalTable <- sort(table(signal), decreasing = T)
(signalTable <- signalTable[signalTable > 30])
signalGene <- sapply(names(signalTable), function(x) names(signal)[signal == x])

op <- par(mar = c(5, 20, 4, 10))
bar <- barplot(signalTable, xlim = c(0, 750), axes = F, border = NA, horiz = T, las = 1, space = 0.75, cex.name = 0.5)
abline(v = 0, lwd = 1, col = "black")
text(y = bar, x = signalTable + 30, labels = signalTable)

binaryApp <- binary[, grep("APP", colnames(binary))]
binaryApp <- binaryApp[apply(binaryApp, 1, any), ]
signalApp <- apply(binaryApp, 1, function (x) paste(x, collapse = "-"))
signalTableApp <- sort(table(signalApp), decreasing = T)
(signalTableApp <- signalTableApp[signalTableApp > 30])
signalGeneApp <- sapply(names(signalTableApp), function(x) names(signalApp)[signalApp == x])

source("../X/function.R")
AppGK1 <- myGK(rownames(binaryApp))
AppGK2 <- lapply(signalGeneApp, myGK)

# STAT3-interacting proteins
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

interHg <- read.delim("R21/interaction.txt", stringsAsFactors = F, header = F)
interMs <- getLDS("hgnc_symbol", "hgnc_symbol", interHg$V1, human, attributesL = "external_gene_name", martL = mouse)
(interMs <- interMs[, 2])

intersect(interMs, rownames(binary))
intersect(interMs, rownames(binaryApp))
lapply(signalGeneApp, function(x) intersect(interMs, x))

keggFind("pathway", "STAT")  # KEGG: Jak-STAT signaling pathway
mmu04630 <- keggGet("mmu04630")[[1]]$GENE
mmu04630 <- gsub(";.*", "", matrix(mmu04630, byrow = T, ncol = 2)[, 2])

keggFind("pathway", "Alzheimer")  # KEGG: Alzheimer's disease
mmu05010 <- keggGet("mmu05010")[[1]]$GENE
mmu05010 <- gsub(";.*", "", matrix(mmu05010, byrow = T, ncol = 2)[, 2])

write(keggGet("hsa05010", "kgml"), file = "hsa05010.xml")  # Cytoscape
writePNG(keggGet("hsa05010", "image"), "hsa05010.png")

bgMs <- gsub("[,|;].*", "", keggList("mmu"))
myhyper <- function(g1, g2) {  # Hypergeometric
  if(length(intersect(g1, g2)) == 0) return(1)
  1 - phyper(length(intersect(g1, g2)) - 1, length(g2), length(setdiff(bgMs, g2)), length(g1))
}  # Pr(count >= length(intersect(g1, g2)))

myhyper(rownames(binaryApp), mmu05010)
myhyper(mmu05010, rownames(binaryApp))

lapply(signalGeneApp, function(x) myhyper(x, mmu05010))
lapply(signalGeneApp, function(x) myhyper(x, mmu04630))

# STAT3 targetome (iRegulon)
targets = read.delim("R21/STAT3_meta.csv", stringsAsFactors = F, sep = ",") %>% filter(Strength > 10)
targets <- getLDS("hgnc_symbol", "hgnc_symbol", targets$Target.Gene, human, attributesL = "external_gene_name", martL = mouse)
lapply(signalGeneApp, function(x) intersect(targets[, 2], x))

targets <- targets %>% filter(Associated.Gene.Name %in% rownames(binaryApp))
igraph.dt <- graph.data.frame(cbind("Stat3", targets$Associated.Gene.Name))  # IGRAPH
igraph.dt$layout <- layout.fruchterman.reingold
V(igraph.dt)$color = rep("chartreuse3", length(V(igraph.dt)$name))
V(igraph.dt)$color[V(igraph.dt)$name %in% mmu04630] <- "black"
V(igraph.dt)$color[V(igraph.dt)$name %in% "Stat3"] <- "gold"
V(igraph.dt)$label.color = rep("dodgerblue3", length(V(igraph.dt)$name))
V(igraph.dt)$label.color[V(igraph.dt)$color == "gold"] = "firebrick1"
plot.igraph(igraph.dt, vertex.frame.color = "white", edge.arrow.size = 0.1, vertex.size = 5, vertex.label.cex = 0.7)

# Comparisons in details
age_spe <- lapply(c("2m", "4m", "5m", "6m"), function(x1) {
  sapply(glm.fit, function(x2) 
    contrast(x2, list(age = x1, group = "APP"), list(age = x1, group = "WT")))[c("Contrast", "Pvalue"), ]
}); names(age_spe) <- c("2m", "4m", "5m", "6m")

lapply(age_spe, function(x) x[, grep("Stat", colnames(x))])

combn <- combn(c("2m", "4m", "5m", "6m"), 2)
geno_spe_wt <- lapply(1:ncol(combn), function(x1) {
  sapply(glm.fit, function(x2)
    contrast(x2, list(age = combn[2, x1], group = "WT"), list(age = combn[1, x1], group = "WT")))[c("Contrast", "Pvalue"), ]
}); names(geno_spe_wt) <- apply(combn, 2, function(x) paste0(x[2], "_vs_", x[1])) 

lapply(geno_spe_wt, function(x) x[, grep("Stat", colnames(x))])

sapply(1:ncol(combn), function(x1) {
  contrast(glm.fit[["Stat3"]], list(age = combn[2, x1], group = "WT"), list(age = combn[1, x1], group = "WT")) 
})[c("Contrast", "Pvalue"), ]

C <- matrix(c(0,0,0,0,0,0,-1,1), 1) # 6m APP vs 5m APP
C <- matrix(c(0,0,0,0,0,-1,0,1), 1) # 6m APP vs 4m APP
C <- matrix(c(0,0,0,0,-1,0,0,1), 1) # 6m APP vs 2m APP
summary(glht(glm.fit[["Stat3"]], linfct = C))  

# # Reactome & Common Pathway
# library(reactome.db)
# xx <- as.list(reactomePATHID2NAME)
# xx[which(unlist(lapply(xx, function(x) grepl("STAT", x))))]
# library(paxtoolsr)
# searchResults <- searchPc(q = "STAT3", type = "pathway")
# searchResultsDf <- ldply(xmlToList(searchResults), data.frame)
# 
# # PCA expression matrix on selected gene groups
# pca.dt <- lapply(signalGene, function(x) sapply(grp, function(y) rowMeans(glm.dt[x, ][, uid == y])))
# pca.dt <- lapply(aovGene, function(x) sapply(grp, function(y) rowMeans(glm.dt[names(x), ][, uid == y])))
# pca <- prcomp(pca.dt[[3]] - rowMeans(pca.dt[[3]]))
# barplot(pca$sdev)
# barplot(pca$rotation[, 1])
# barplot(pca$rotation[, 2])
# plot(pca$rotation[, 1], pca$rotation[, 2])
#   