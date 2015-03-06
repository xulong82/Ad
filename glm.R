# Copyright: Xulong Wang (xulong.wang@jax.org)
# Generalized Linear Model (GLM)

library(ape)
library(amap)
library(ggdendro)
library(VennDiagram)
library(MASS)
library(xtable)
library(ggplot2)
library(ggvis)
library(pheatmap)
library(grid)
library(xlsx)

rm(list = ls())
setwd("~/Dropbox/GitHub/Ad")
load("./data/brain2014.rdt")
load("./data/retina2014.rdt")
source("function.R")

dt <- brain.tpm
dt <- retina.tpm

cutoff <- quantile(c(as.matrix(dt)), 0.25)  # TPM level
dt <- dt[apply(dt, 1, function(x) max(x) > cutoff), ]
dt <- dt[apply(dt, 1, function(x) sum(x > 0) > round(ncol(dt) / 10)), ]

age <- factor(gsub("^.*(2m|4m|5m|6m).*", "\\1", colnames(dt)), levels = c("2m", "4m", "5m", "6m"))
group <- factor(gsub("^.*(WT|APP).*", "\\1", colnames(dt)), levels = c("WT", "APP"))
batch <- factor(gsub("^.*(2014|mouse).*", "\\1", colnames(dt)), levels = c("2014", "mouse"))
uid <- paste(age, group, sep = "_")
conditions <- c("2m_WT", "2m_APP", "4m_WT", "4m_APP", "5m_WT", "5m_APP", "6m_WT", "6m_APP")

# --- GLM regression ---
glm.dt <- as.matrix(dt)
glm.fit <- lapply(1:nrow(glm.dt), function (x) summary(lm(glm.dt[x, ] ~ age + group + batch + age * group)))
names(glm.fit) <- rownames(glm.dt)
dt.bc <- dt
for (i in 1:nrow(dt.bc)) 
  dt.bc[i, ] <- dt[i, ] - glm.fit[[i]]$coefficients["batchmouse", "Estimate"] * (as.numeric(batch) - 1)
save(dt.bc, age, group, uid, conditions, file = "./data/brain_bc2014.rdt")
save(dt.bc, age, group, uid, conditions, file = "./data/retina_bc2014.rdt")

r2 <- sapply(glm.fit, function (x) x$r.squared)
fval <- sapply(glm.fit, function (x) x$fstatistic)
pval <- apply(fval, 2, function (x) pf(x[1], x[2], x[3], lower.tail = F))
qval <- p.adjust(pval, method = "fdr")

fit.qr <- glm.fit[qval < 0.05 & r2 > 0.5]
fit.qr <- glm.fit[qval < 0.05 & r2 > 0.4]

fit.est <- lapply(fit.qr, function (x) x$coefficients[, "Estimate"])
fit.est <- do.call(rbind, fit.est)
fit.pval <- lapply(fit.qr, function (x) x$coefficients[, "Pr(>|t|)"])
fit.pval <- do.call(rbind, fit.pval)
  
# --- GLM cohorts ---
geneId <- rownames(fit.pval)
logit <- apply(fit.pval, 2, function (x) x < 0.05) & apply(fit.est, 2, function (x) abs(x) > 0.2)
logit <- apply(fit.pval, 2, function (x) x < 0.05) & apply(fit.est, 2, function (x) abs(x) > 0.1)
geneId.app <- geneId[as.logical(rowSums(logit[, grep("APP", colnames(fit.pval))]))]
geneId.age <- geneId[as.logical(rowSums(logit[, grep("age", colnames(fit.pval))]))]

gk.app <- myGK(geneId.app)
gk.age <- myGK(geneId.age)

profile <- logit[geneId.app, grep("APP", colnames(logit))]
profile <- logit[geneId.age, grep("age", colnames(logit))]

profile.str <- apply(profile, 1, function (x) paste(x, collapse = "-"))
profile.table <- sort(table(profile.str))
profile.Id <- names(profile.table)

geneId.profile <- list()
for (idx in profile.Id)
  geneId.profile[[idx]] <- geneId.app[profile.str == idx]
  geneId.profile[[idx]] <- geneId.age[profile.str == idx]

glm <- list()
save(glm, file = "data/glm_brain.rdt")
save(glm, file = "data/glm_retina.rdt")
  
glm$less$app$symbol <- geneId.app
glm$less$app$go_bp <- gk.app$BP
glm$less$app$go_mf <- gk.app$MF
glm$less$app$go_cc <- gk.app$CC
glm$less$app$go_kegg <- gk.app$KEGG
glm$less$app$symbol_by_pattern <- geneId.profile

glm$less$age$symbol <- geneId.age
glm$less$age$go_bp <- gk.age$BP
glm$less$age$go_mf <- gk.age$MF
glm$less$age$go_cc <- gk.age$CC
glm$less$age$go_kegg <- gk.age$KEGG
glm$less$age$symbol_by_pattern <- geneId.profile

glm$more$app$symbol <- geneId.app
glm$more$app$go_bp <- gk.app$BP
glm$more$app$go_mf <- gk.app$MF
glm$more$app$go_cc <- gk.app$CC
glm$more$app$go_kegg <- gk.app$KEGG
glm$more$app$symbol_by_pattern <- geneId.profile

glm$more$age$symbol <- geneId.age
glm$more$age$go_bp <- gk.age$BP
glm$more$age$go_mf <- gk.age$MF
glm$more$age$go_cc <- gk.age$CC
glm$more$age$go_kegg <- gk.age$KEGG
glm$more$age$symbol_by_pattern <- geneId.profile
  
grid.newpage()                                                                                                                                       
draw.pairwise.venn(length(geneId.age), length(geneId.app), length(intersect(geneId.age, geneId.app)), 
                   category = c("Age", "APP"), fill = c("light blue", "pink"))
                                                                                                                                                        
tile.dt <- NULL
for (i in 1:length(profile.Id))
  tile.dt <- rbind(tile.dt, as.logical(unlist(strsplit(profile.Id[i], "-"))))
tile.dt <- data.frame(value = c(tile.dt), 
  profile = factor(rep(profile.Id, ncol(tile.dt)), levels = profile.Id),
  group = factor(rep(colnames(profile), each = nrow(tile.dt)), levels = colnames(profile)))

ggplot(tile.dt, aes(x = group, y = profile, fill = value)) + geom_tile(colour = "white") +
  theme_bw() + xlab("") + ylab("") + coord_flip() +
  scale_x_discrete(labels = c("2m:APP", "4m:APP", "5m:APP", "6m:APP")) +
# scale_x_discrete(labels = c("4m:WT", "5m:WT", "6m:WT", "4m:APP", "5m:APP", "6m:APP")) +
  scale_y_discrete(labels = profile.table) +
  scale_fill_manual(values = c("grey80", "firebrick1")) 

# --- HC on GLM cohorts ---
dt.hc <- dt.bc[unique(c(geneId.app, geneId.age)), ]

susId <- c("Gm13394", "Lamr1-ps1", "Gnat1", "Rho")
dt.sus <- dt.bc[susId, ]

hc1 <- hcluster(t(dt.hc), method = "pearson", link = "average")
hc2 <- hcluster(dt.hc, method = "pearson", link = "average")

mycol <- rep("grey50", ncol(dt.hc))
mycol[grep("APP", colnames(dt.hc))] <- "firebrick1"
plot(as.phylo(hc1), edge.width = 2, font = 2, cex = 0.7, label.offset = 1e-4, tip.color = mycol, direction = "downward")

group1<- c("APP5m1558.2014", "APP5m2751.2014", "APP5m1633.2014","APP5m1636.2014","APP6m1684.2014","APP5m1648.2014", "APP5m1623.2014")
group2<- c("mouse_3346_6m_APP", "mouse_3440_6m_APP", "mouse_3351_6m_APP", "mouse_6585_6m_APP", "APP5m1647.2014", "APP5m1738.2014")

dt.g1 <- dt.bc[, group1]
dt.g2 <- dt.bc[, group2]
colnames(dt.g1) <- paste("group1", colnames(dt.g1), sep = "_")
colnames(dt.g2) <- paste("group2", colnames(dt.g2), sep = "_")
dt.g12 <- cbind(dt.g1, dt.g2)
fc <- rowMeans(dt.g1) - rowMeans(dt.g2) 
treat <- gsub("^.*(group1|group2).*", "\\1", colnames(dt.g12))
tt.pval <- apply(dt.g12, 1, function(x) pairwise.t.test(x, treat, p.adj = "none")$p.value)
tt.qval <- p.adjust(tt.pval, method = "fdr")
idx <- which(tt.qval < 0.05 & abs(fc) > 0.2)
geneId <- rownames(dt.g12)[idx]

gk.5m <- myGK(geneId)

# --- PCA GLM estimate
svd.dt <- fit.est[, -c(1, 6)]
feature <- factor(colnames(svd.dt), levels = colnames(svd.dt))
geneId <- rownames(svd.dt)
PC <- paste("PC", 1:ncol(svd.dt), sep = "")
svd <- svd(svd.dt)
barplot(svd$d)
rownames(svd$v) <- feature
colnames(svd$v) <- PC

gdt <- data.frame(value = c(svd$v), PC = rep(PC, each = 7), feature = rep(feature, 7))
gdt$geno <- rep("WT", nrow(gdt))
gdt$geno[grep("APP", gdt$feature)] <- "APP"
ggplot(gdt, aes(x = feature, y = value, fill = geno)) + 
  geom_bar(stat = "identity") + facet_grid(. ~ PC) +
  scale_fill_manual(values = c("red", "blue")) +
  theme_bw() + xlab("") + ylab("") + 
  theme(axis.text.x = element_text(angle = 90))

barplot(svd$u[, 1]) 
cut <- quantile(abs(svd$u[, 1]), 0.95)
abline(h = -cut, col = "red"); abline(h = cut, col = "red")

dt.cond <- NULL
for (idx in conditions) dt.cond <- cbind(dt.cond, rowMeans(dt.bc[, uid == idx]))
colnames(dt.cond) <- conditions

x = geneId[abs(svd$u[, 1]) > cut]
x.dt = dt.cond[x, ]
gdt <- data.frame(value = c(x.dt), gene = rep(rownames(x.dt), 8), group = rep(conditions, each = nrow(x.dt)))
gdt %>% ggvis(~group, ~value, stroke = ~gene) %>% layer_lines() %>%
  add_tooltip(function (x) paste("Gene: ", x$gene), "hover")

mycol <- c("grey70", "firebrick1", "chartreuse3", "dodgerblue3", "gold1", "darkorchid2")
clusts = cutree(hc2, 6)
plot(as.phylo(hc2), type = "unrooted", tip.color = mycol[clusts], cex = 0.5, font = 2, lab4ut = "axial")
