# Copyright: Xulong Wang (xulong.wang@jax.org)

library(ape)
library(amap)
library(ggdendro)
library(MASS)
library(xtable)
library(ggplot2)
library(lattice)
library(pheatmap)
library(grid)
library(mygene)

rm(list = ls())
setwd("~/Dropbox/AD")
load("./R/brain2014.rdt")

dt <- brain.tpm
# dt <- dt[apply(dt, 1, function(x) max(x) - min(x) > 1), ]
cutoff <- quantile(c(as.matrix(dt)), 0.25)  # TPM level
dt <- dt[apply(dt, 1, function(x) max(x) > cutoff), ]
dt <- dt[apply(dt, 1, function(x) sum(x > 0) > round(ncol(dt) / 10)), ]

group <- factor(gsub("^.*(WT|APP).*", "\\1", colnames(dt)), levels = c("WT", "APP"))
age <- factor(gsub("^.*(2m|4m|5m|6m).*", "\\1", colnames(dt)), levels = c("2m", "4m", "5m", "6m"))
batch <- factor(gsub("^.*(2014|mouse).*", "\\1", colnames(dt)), levels = c("2014", "mouse"))
uid <- paste(age, group, sep = "_")
uid.level <- c("2m_WT", "2m_APP", "4m_WT", "4m_APP", "5m_WT", "5m_APP", "6m_WT", "6m_APP")

# --- model: regression ---
myfit <- list()
rg.dt <- as.matrix(dt)
for (i in 1:nrow(rg.dt)) {
  if (i %% 1e3 == 0) cat(i, "\n")
  y <- rg.dt[i, ]
  fit0 <- lm(y ~ age + group + batch + age * group)
  myfit[[i]] <- summary(fit0)
}

dt.bc <- dt
for (i in 1:nrow(dt.bc)) 
  dt.bc[i, ] <- dt[i, ] - myfit[[i]]$coefficients["batchmouse", "Estimate"] * (as.numeric(batch) - 1)
  
save(dt.bc, group, age, uid, uid.level, file = "./Shiny/dt_batch_correction.rdt")

names(myfit) <- rownames(rg.dt)
r2 <- sapply(myfit, function (x) x$r.squared)
fval <- sapply(myfit, function (x) x$fstatistic)
pval <- apply(fval, 2, function (x) pf(x[1], x[2], x[3], lower.tail = F))
qval <- p.adjust(pval, method = "fdr")

geneId <- names(myfit)[qval < 0.05 & r2 > 0.5]
fit <- myfit[geneId]

fit.pval <- NULL
fit.est <- NULL
for (id in geneId) {
  fit.pval <- rbind(fit.pval, fit[[id]]$coefficients[, "Pr(>|t|)"])
  fit.est <- rbind(fit.est, fit[[id]]$coefficients[, "Estimate"])
} 
rownames(fit.pval) <- rownames(fit.est) <- geneId

sel <- apply(fit.pval, 2, function (x) x < 0.05) & apply(fit.est, 2, function (x) abs(x) > 0.2)
geneId.app <- geneId[as.logical(rowSums(sel[, grep("APP", colnames(fit.pval))]))]
geneId.age <- geneId[as.logical(rowSums(sel[, grep("age", colnames(fit.pval))]))]
write.table(geneId.app, file = "./Gene/geneId_app.txt", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(geneId.age, file = "./Gene/geneId_age.txt", sep = "\t", row.names = F, col.names = F, quote = F)

profile.app <- sel[geneId.app, grep("APP", colnames(sel))]
profile.app.str <- apply(profile.app, 1, function (x) paste(x, collapse = "-"))
profile.app.table <- sort(table(profile.app.str))
profile.app.Id <- names(profile.app.table)

geneId.app.profile <- list()
for (idx in profile.app.Id)
  geneId.app.profile[[idx]] <- geneId.app[profile.app.str == idx]

lapply(geneId.app.profile, write, "./Gene/glm_app_profile.txt", append = TRUE, ncolumns = 1e3)
  
profile.age <- sel[geneId.age, grep("age", colnames(sel))]
profile.age.str <- apply(profile.age, 1, function (x) paste(x, collapse = "-"))
profile.age.table <- sort(table(profile.age.str))
profile.age.Id <- names(profile.age.table)
  
geneId.age.profile <- list()
for (idx in profile.age.Id)
  geneId.age.profile[[idx]] <- geneId.age[profile.age.str == idx]
  
lapply(geneId.age.profile, write, "./Gene/glm_age_profile.txt", append = TRUE, ncolumns = 1e3)
  
mus2hg <- read.delim("~/Dropbox/X/hg2mus.map", header = F, stringsAsFactors = F)
annotation <- queryMany(geneId.app, scopes="symbol", species="mouse", 
  return.as = "DataFrame", fields = c("name", "summary", "go", "kegg"))
annotation <- queryMany(geneId.age, scopes="symbol", species="mouse", 
  return.as = "DataFrame", fields = c("name", "summary", "go", "kegg"))
table <- annotation[c("query", "name", "summary")]
table <- table[!duplicated(table$query), ]

table$hg <- mus2hg$V1[match(geneId.app, mus2hg$V4)]
table$hg <- mus2hg$V1[match(geneId.age, mus2hg$V4)]
annotation.hg <- queryMany(table$hg, scopes="symbol", species="human", 
  return.as = "DataFrame", fields = c("name", "summary", "go", "kegg"))
table.hg <- annotation.hg[c("query", "name", "summary")]
table.hg <- table.hg[!duplicated(table.hg$query), ]
table$summary_hg <- table.hg$summary[match(table$hg, table.hg$query)]
table$summary <- paste("Mus", table$summary, sep = ":")
table$summary_hg <- paste("Hg", table$summary_hg, sep = ":")
table$summary <- paste(table$summary, table$summary_hg)
table <- as.data.frame(table[c("query", "name", "hg", "summary")])

write.table(table, file = "./Gene/table_app.txt", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(table, file = "./Gene/table_age.txt", sep = "\t", row.names = F, col.names = T, quote = F)

write(print(xtable(table)), file = "./Gene/summary_app.tex", append = T)
write(print(xtable(table)), file = "./Gene/summary_age.tex", append = T)

# --- cluster ---
dt.hc <- dt.bc[geneId, ]
dt.hc <- dt.bc[unique(c(geneId.app, geneId.age)), ]

susId <- c("Gm13394", "Lamr1-ps1", "Gnat1", "Rho")
dt.sus <- dt.bc[susId, ]
# dt.hc <- dt.hc[! rownames(dt.hc) %in% susId, ]

hc1 <- hcluster(t(dt.hc), method = "pearson", link = "average")
hc2 <- hcluster(dt.hc, method = "pearson", link = "average")

plot(as.phylo(hc1), edge.width = 2, font = 2, cex = 0.7, label.offset = 1e-4, direction = "downward")

conditions <- unique(uid)
byConditions <- NULL
for (idx in conditions) byConditions <- cbind(byConditions, rowMeans(dt.bc[, uid == idx]))
colnames(byConditions) <- conditions

tef <- byConditions["Tef", ]
c1qa <- byConditions["C1qa", ]

PlotGroups(data.abiotic["STMCL34", ], edesign = edesign.abiotic)
see.genes(sigs$sig.genes$ColdvsControl, main = "ColdvsControl", show.fit = T, dis = design$dis, 
          cluster.method = "kmeans", cluster.data = 1, k = 9)

# --- ANOVA ---
aov.pval <- NULL
for (i in 1:nrow(dt)) {
  if (i %% 1e3 == 0) cat(i, "\n")
  y <- dt[i, ]
  aov <- aov(y ~ age * group * batch)
  aov.pval <- c(aov.pval, min(summary(aov)[[1]][["Pr(>F)"]], na.rm = T))
}
dt.aov <- dt[aov.pval < 0.05, ]

list <- read.delim("./list.txt", header = F, stringsAsFactors = F)
list$V1[list$V1 %in% geneId.de]
