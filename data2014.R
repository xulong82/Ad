# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: AD project - clustering

library(ape)
library(amap)
library(mygene)
library(ggplot2)
library(ggdendro)
library(pheatmap)

rm(list = ls())
setwd("~/Dropbox/GitHub/Ad")
load(file = "./data/complete_tpm.rdt")
col.manual <- c("grey70", "firebrick1", "chartreuse3", "dodgerblue3", "gold1", "darkorchid2")

# --- 2014 brain samples
brain.tpm <- brain.tpm[, c(grep("2014", colnames(brain.tpm)), grep("mouse", colnames(brain.tpm)))]
brain.tpm <- brain.tpm[apply(brain.tpm, 1, function (x) max(x) > 5), ]
brain.tpm <- log2(brain.tpm + 1)

retina.tpm <- retina.tpm[apply(retina.tpm, 1, function (x) max(x) > 5), ]
retina.tpm <- log2(retina.tpm + 1)

mycol <- rep("blue", ncol(brain.tpm))
mycol[grep("mouse", colnames(brain.tpm))] <- "red"
hc1 <- hcluster(t(brain.tpm), method = "pearson", link = "average")
par(mar = c(0, 0, 0, 0))
plot(as.phylo(hc1), type = "unrooted", tip.color = mycol, lab4ut = "axial")
save(hc1, mycol, file = "./markdown/brain_hc.rdt")

mycol <- rep("blue", ncol(retina.tpm))
mycol[grep("mouse", colnames(retina.tpm))] <- "red"
hc1 <- hcluster(t(retina.tpm), method = "pearson", link = "average")
par(mar = c(0, 0, 0, 0))
plot(as.phylo(hc1), type = "unrooted", cex = 1, tip.color = mycol, lab4ut = "axial")

# --- mouse 1559
brain.tpm1 <- brain.tpm[apply(brain.tpm, 1, function(x) max(x) - min(x) > 1), ]
par(mfrow = c(2, 1), mar = c(2, 4, 1, 2))
plot(brain.tpm1$WT6m1484.2014 - brain.tpm1$APP5m1558.2014, ylim = c(-6, 6), xlab = "", ylab = "")
title("WT6m1484 - APP5m1558", line = -2, col.main = "firebrick1") +
abline(0, 0, lwd = 1, col = "firebrick1")
plot(brain.tpm1$APP5m1559.2014 - brain.tpm1$APP5m2633.2014, ylim = c(-6, 6), xlab = "", ylab = "") 
title("APP5m1559 - APP5m2633", line = -2, col.main = "firebrick1" ) +
abline(0, 0, lwd = 1, col = "firebrick1")
save(brain.tpm1, file = "./markdown/mouse1559.rdt")

brain.tpm <- brain.tpm[, -grep("APP5m1559.2014", colnames(brain.tpm))]
dev.off()
brain.hist <- hist(c(as.matrix(brain.tpm)))
save(brain.hist, file = "./markdown/brain2014.rdt")

save(brain.tpm, file = "./data/brain2014.rdt")
save(retina.tpm, file = "./data/retina2014.rdt")

# --- Annotation
geneId <- union(rownames(brain.tpm), rownames(retina.tpm))
mus2hg <- read.delim("~/Dropbox/X/hg2mus.map", header = F, stringsAsFactors = F)
table <- queryMany(geneId, scopes="symbol", species="mouse", fields = c("name", "summary"))
table <- table[!duplicated(table$query), ]
table$hg <- mus2hg$V1[match(geneId, mus2hg$V4)]
query.hg <- queryMany(table$hg, scopes="symbol", species="human", fields = c("name", "summary"))
query.hg <- query.hg[!duplicated(query.hg$query), ]
summary_hg <- query.hg$summary[match(table$hg, query.hg$query)]
table$summary <- paste("Mus:", table$summary, "Hg:", summary_hg)
table <- table[c("query", "name", "hg", "summary")]
rownames(table) <- NULL
table <- as.data.frame(table)
save(table, file = "./data/gene_summary.rdt")

# --------------------------------------------------------------------
# rm(list = ls())
# load("./R/brain2014.rdt")
# dt <- brain.tpm
# 
# # --- batch correction
# batch <- rep("2014", ncol(dt))
# batch[grep("mouse", colnames(dt))] <- "2014-new"
# 
# dt.mean <- apply(dt, 1, mean)
# dt.center <- dt - dt.mean
# dt.res <- t(apply(dt.center, 1, function (x) {lm(x ~ as.factor(batch) - 1)$res}))
# dt <- dt.res + dt.mean
# 
# save(dt, file = "./R/batch2014.rdt")
# 
# mycol <- rep(col.manual[2], ncol(dt))
# mycol[grep("4m", colnames(dt))] <- col.manual[3]
# mycol[grep("5m", colnames(dt))] <- col.manual[4]
# mycol[grep("6m", colnames(dt))] <- col.manual[5]
# 
# dt <- dt[apply(dt, 1, function(x) max(x) - min(x) > 1), ] 
# hc1 <- hcluster(t(dt), method = "pearson", link = "average")  # clust samples
# pdf("./Graphs/batch2014dendro.pdf")
# plot(as.phylo(hc1), type = "unrooted", tip.color = mycol, cex = .5, font = 2, lab4ut = "axial")
# dev.off()
# 
# load("./R/batch2014.rdt")
# # --- ANOVA --- ANOVA's return confusing
# treat <- gsub("^.*(WT|APP).*", "\\1", colnames(dt))
# month <- gsub("^.*(2m|4m|5m|6m).*", "\\1", colnames(dt))
# aov.pval <- NULL
# for (i in 1:nrow(dt)) {
#   if (i %% 1e3 == 0) cat(i, nrow(dt), "\n")
#   aov1 = aov(as.matrix(dt)[i, ] ~ month * treat)
#   aov.pval <- c(aov.pval, min(summary(aov1)[[1]][["Pr(>F)"]], na.rm = T))
# }
# dt.aov <- dt[aov.pval < 0.05, ]  # Select genes
# 
# dt.fc <- dt[apply(dt, 1, function(x) max(x) - min(x) > 1), ]  # FC works better
# 
# dt <- dt.aov
# dt <- dt.fc
# 
# # ---
# dt <- dt[, grep("6m", colnames(dt))]
# hc1 <- hcluster(t(dt), method = "pearson", link = "average")  # clust samples
# plot(as.phylo(hc1), type = "unrooted", cex = .5, font = 2, lab4ut = "axial")
# # ---
# hc1 <- hcluster(t(complete456m), method = "pearson", link = "average")
# plot(as.phylo(hc1), type = "unrooted", cex = .5, font = 2, lab4ut = "axial")
# 
# # ---
# rm(list = ls())
# load(file = "~/Dropbox/AD/R/complete_tpm.rdt")
# load(file = "~/Dropbox/AD/R/cluster1.rdt")
# 
# DEG <- unique(c(rownames(DE.month$DE4m), rownames(DE.month$DE5m), rownames(DE.month$DE6m)))
# dt <- dt[DEG, ]
# 
# hc1 <- hcluster(t(dt), method = "pearson", link = "average")  # clust samples
# plot(as.phylo(hc1), type = "unrooted", cex = .5, font = 2, lab4ut = "axial")
# 
# dt <- dt[, grep("APP", colnames(dt))]
# col.manual <- c("firebrick1", "chartreuse3", "dodgerblue3", "gold1", "darkorchid2")
# hc1 <- hcluster(t(dt), method = "pearson", link = "average")
# clusts = cutree(hc1, 3)
# pdf(file = "~/Dropbox/AD/Graphs/cluster2_phylo.pdf")
# plot(as.phylo(hc1), type = "unrooted", tip.color = col.manual[clusts], cex = .5, font = 2, lab4ut = "axial")
# dev.off()
# 
# #-------------
# tpm1 <- brain.tpm2
# id1 <- which(apply(tpm1, 1, function(x) {min(x) < 0}))  # negative TPM from correction
# tpm1[id1, ] <- tpm1[id1, ] + abs(apply(tpm1[id1, ], 1, min))
# # id2 <- which(rownames(tpm1) %in% c("App", "Prnp"))  # App & Prnp
# # tpm1 <- tpm1[-id2, ]
# tpm2 <- log2(tpm1 + 1)  # log2 transformation
# tpm3 <- tpm2[apply(tpm2, 1, function (x) {min(x) > 3}), ]
# tpm3 <- tpm3[apply(tpm3, 1, function (x) {max(x) - min(x) > 1}), ]
# tpm4 <- tpm3 - apply(tpm3, 1, mean)  # center
# # tpm3.sd <- apply(tpm3, 1, sd)
# # tpm4 <- tpm3[tpm3.sd > quantile(tpm3.sd, 0.75), ]  # 25% of high SD 
# 
# tpm5 <- tpm4[, grep("WT", colnames(tpm4))]  # WT
# tpm6 <- tpm4[, grep("APP", colnames(tpm4))]  # APP
# tpm7 <- tpm4[, grep("5m", colnames(tpm4))]  # 5 month
# 
# hc1 <- hcluster(t(tpm7), method = "pearson", link = "average")  # clust samples
# hc2 <- hcluster(tpm7, method = "pearson", link = "average")  # clust genes
# 
# plot(as.phylo(hc1), type = "unrooted", lab4ut = "axial")
# 
# tile1.dat1 <- tpm3[hc2$order, hc1$order]
# # tile1.dat1 <- tile1.dat1[1:1000, ]
# tile1.dat2 <- t(apply(tile1.dat1, 1, rank))
# tile1.dat3 <- data.frame(value = c(tile1.dat2), 
#                          gene = rep(rownames(tile1.dat2), time = ncol(tile1.dat2)),
#                          sample = rep(colnames(tile1.dat2), each = nrow(tile1.dat2)))
# tile1.dat3$sample <- factor(tile1.dat3$sample, levels = colnames(tile1.dat1))
# tile1.dat3$gene <- factor(tile1.dat3$gene, levels = rownames(tile1.dat1))
# 
# pdf("~/Dropbox/AD/Figures/tile1.pdf", width = 15, height = 5)
# ggplot(tile1.dat3, aes(x = gene, y = sample, fill = value)) + 
#   geom_tile() + scale_fill_gradient(low="green", high="red") +
#   theme_bw() +
#   xlab("") + ylab("") +
#   theme(panel.border = element_rect(size = 1, color = "black")) +
#   theme(axis.text.x = element_blank(),
#         axis.text.y = element_text(size = 5),
#         axis.ticks = element_blank()) +
#   theme(legend.position = "none")
# dev.off()

# #--- Archieve
# colnames(tpm) <- gsub("m.*", "m", colnames(tpm))
# time <- c("2m", "4m", "5m", "6m")
# 
# tpm.wt <- tpm[, grep("WT", colnames(tpm))]
# tpm.app <- tpm[, grep("APP", colnames(tpm))]
# 
# tpm.wt1 <- matrix(nrow = nrow(tpm.wt), ncol = length(time), dimnames = list(rownames(tpm.wt), time))
# tpm.app1 <- matrix(nrow = nrow(tpm.app), ncol = length(time), dimnames = list(rownames(tpm.app), time))
# for (i in time) {  # app.mean, wt.mean
#   tpm.wt1[, i] <- rowMeans(tpm.wt[, grep(i, colnames(tpm.wt))])
#   tpm.app1[, i] <- rowMeans(tpm.app[, grep(i, colnames(tpm.app))])
# }
# tpm.app2 <- tpm.app1 - tpm.wt1  # app.mean - wt.mean
# 
# tpm.app3 <- tpm.app
# for (i in 1:ncol(tpm.app)) {  # app - wt.mean
#   idx <- gsub("APP", "", colnames(tpm.app)[i])
#   tpm.app3[, i] <- tpm.app[, i] - tpm.wt1[, grep(idx, colnames(tpm.wt1))] 
# }
# 
# #--- Kmeans
# app2.svd <- svd(tpm.app2[idx, ])  # SVD
# plot(app2.svd$d, type = "b")  # dimension reduction
# 
# nc.max <- 30
# withinss <- cbind(nc = c(2:nc.max), withinss = rep(NA, nc.max - 1))
# for (nc in 2:nc.max) {  # choose K based on within group ss
#   cat(nc, "/", nc.max, "\n")
# # app2.kmeans <- kmeans(tpm.app2[idx, ], centers = nc, iter.max = 50, nstart = 10)  # Kmeans
#   app2.kmeans <- kmeans(app2.svd$u[, 1:3], centers = nc, iter.max = 50, nstart = 10)  # Kmeans
#   withinss[nc-1, 2] <- sum(app2.kmeans$withinss)
# }
# plot(withinss[, 1:2], type = "b")
# 
# nc <- 20  # number of clusters
# app2.kmeans <- kmeans(app2.svd$u[, 1:3], centers = nc, iter.max = 50, nstart = 10)  # Kmeans
# app2.kmeans <- kmeans(tpm.app2[idx, ], centers = nc, iter.max = 50, nstart = 10)  # Kmeans
# table(app2.kmeans$cluster)
# table(as.vector(table(app2.kmeans$cluster)) > 20)
# 
# pdf("~/Dropbox/AD/Figures/kmean1.pdf", width = 10)
# colors <- palette(rainbow(nc))
# par(mfrow = c(4, 5), mar = c(5, 4, 4, 2))
# for (i in 1:nc) {
#   name <- paste("Cluster", i, sep = " ")
#   length <- table(app2.kmeans$cluster)[i]
#   temp <- hclust.dat2[app2.kmeans$cluster == i, ]
# # temp <- tpm.app2[idx, ][app2.kmeans$cluster == i, ]
#   
#   if (length == 1) {
#     plot(temp, type = "b", col = colors[i], ylim = c(-10, 10),
#          main = name, xlab = "", ylab = "") 
#   } else {
#     plot(colMeans(temp), type = "b", ylim = c(-10, 10),
#          main = name, xlab = "", ylab = "") 
#     apply(temp, 1, function(x) {lines(x, type = "l", col = colors[i])})
#   }
# }
# dev.off()
# 
# i = 2
# temp <- tpm.app2[idx, ][app2.kmeans$cluster == i, ]
# par(mfrow = c(1, 1))
# plot(colMeans(temp), type = "b") 
# apply(temp, 1, function(x) {lines(x, type = "l", col = colors[i])})
# 

#-------
# for (month in c("4m", "5m", "6m")) {
#   cat(month, "\n")
#   
#   tpm1 <- brain.tpm1[, grep(month, colnames(brain.tpm1))]  # select samples by month
#   year <- rep("2014", ncol(tpm1))
#   year[grep("2013", colnames(tpm1))] <- "2013"
#   
#   tpm1.mean <- apply(tpm1, 1, mean)
#   tpm1.center <- tpm1 - tpm1.mean
#   tpm1.svd <- svd(tpm1.center)  # SVD
#   
#   tpm1.res <- t(apply(tpm1.center, 1, function (x) {lm(x ~ as.factor(year) - 1)$res}))  # regress on year
#   # tpm1.res <- t(apply(tpm1.center, 1, function (x) {lm(x ~ tpm1.svd$v[, 1])$res}))  # regress on PC1
#   tpm2 <- tpm1.res + tpm1.mean
#   
#   assign(paste("tpm.rgd", month, sep = "_"), tpm2)
#   assign(paste("tpm.res", month, sep = "_"), tpm1.res)
#   assign(paste("tpm.svd", month, sep = "_"), tpm1.svd)
# }
# 
# #------- plot
# type.4m <- rep("WT", ncol(tpm.res_4m))
# type.4m[grep("APP", colnames(tpm.res_4m))] <- "APP"
# type.5m <- rep("WT", ncol(tpm.res_5m))
# type.5m[grep("APP", colnames(tpm.res_5m))] <- "APP"
# type.6m <- rep("WT", ncol(tpm.res_6m))
# type.6m[grep("APP", colnames(tpm.res_6m))] <- "APP"
# 
# year.4m <- rep("2013", ncol(tpm.res_4m))
# year.4m[grep("2014", colnames(tpm.res_4m))] <- "2014"
# year.5m <- rep("2013", ncol(tpm.res_5m))
# year.5m[grep("2014", colnames(tpm.res_5m))] <- "2014"
# year.6m <- rep("2013", ncol(tpm.res_6m))
# year.6m[grep("2014", colnames(tpm.res_6m))] <- "2014"
# 
# pca1.dat <- data.frame(value = c(tpm.svd_4m$v[, 1], tpm.svd_5m$v[, 1], tpm.svd_6m$v[, 1]),
#                        month = c(rep("4m", ncol(tpm.res_4m)), rep("5m", ncol(tpm.res_5m)), rep("6m", ncol(tpm.res_6m))), 
#                        type =  c(type.4m, type.5m, type.6m),
#                        year =  c(year.4m, year.5m, year.6m),
#                        sample = c(colnames(tpm.res_4m), colnames(tpm.res_5m), colnames(tpm.res_6m)))
# pdf("~/Dropbox/AD/Figures/pca1.pdf", height = 15, width = 25)
# ggplot(pca1.dat, aes(x = sample, y = value, fill = type)) + 
#   geom_bar(stat = "identity") + facet_grid(month ~ year) +
#   theme_bw() +
#   geom_hline(yintercept = 0) +
#   scale_colour_brewer(palette="Set1") +
#   xlab("") + ylab("PCA1") +
#   theme(panel.border = element_rect(size = 2, color = "black")) +
#   theme(axis.text = element_text(size = 0, angle = -90, face = "bold"),
#         axis.title = element_text(size = 20, face = "bold"),
#         strip.text = element_text(size = 20, face = "bold")) +
#   theme(legend.position = "top", legend.direction = "horizontal", 
#         legend.text = element_text(size = 20, face = "bold"),
#         legend.title = element_blank(), legend.key = element_blank())
# dev.off()
# 
# tpm.2m <- as.matrix(brain.tpm1[, grep("2m", colnames(brain.tpm1))])  # No 2m sample in 2013
# tpm.ctr_2m <- tpm.2m - apply(tpm.2m, 1, mean)
# 
# tpm.res <- cbind(tpm.res_4m, tpm.res_5m, tpm.res_6m)
# tpm.rgd <- cbind(tpm.rgd_4m, tpm.rgd_5m, tpm.rgd_6m)
# 
# brain.tpm2 <- cbind(tpm.2m, tpm.rgd)
# save(brain.tpm2, file = "~/Dropbox/AD/R/brain.tpm2.rdt")
# 
# pdf("~/Dropbox/AD/Figures/heatmap2m.pdf", height =5, width = 5)
# pheatmap(cor(tpm.ctr_2m, method = "pearson"), display_number = T, treeheight_row = 0)
# dev.off()
# 
# pdf("~/Dropbox/AD/Figures/heatmap4m.pdf", height = 10, width = 10)
# pheatmap(cor(tpm.res_4m, method = "pearson"), display_number = T, treeheight_row = 0)
# dev.off()
# 
# pdf("~/Dropbox/AD/Figures/heatmap5m.pdf", height = 15, width = 15)
# pheatmap(cor(tpm.res_5m, method = "pearson"), display_number = T, treeheight_row = 0)
# dev.off()
# 
# pdf("~/Dropbox/AD/Figures/heatmap6m.pdf", height = 8, width = 8)
# pheatmap(cor(tpm.res_6m, method = "pearson"), display_number = T, treeheight_row = 0)
# dev.off()
# 
# pdf("~/Dropbox/AD/Figures/heatmap6m1.pdf", height = 8, width = 8)
# pheatmap(cor(tpm[, grep("6m", colnames(tpm))], method = "pearson"), display_number = T, treeheight_row = 0)
# dev.off()
# 
# pdf("~/Dropbox/AD/Figures/heatmap456m.pdf", height = 15, width = 15)
# pheatmap(cor(tpm.rgd, method = "spearman"), display_number = T, 
#          treeheight_row = 0, treeheight_col = 200, fontsize = 10, fontsize_number = 5)
# dev.off()
# 
# pca3.dat <- tpm2.svd$v[, c(1, 4)]
# rownames(pca3.dat) <- colnames(tpm2)
# colnames(pca3.dat) <- c("PC1", "PC4")
# pca3.dat <- as.data.frame(pca3.dat)
# pca3.dat$group <- gsub("m.*", "m", colnames(tpm2))
# mycolors <- c("darkorchid4", "darkgreen", "blue2", "gold", "chocolate1", "brown2")
# pdf("~/Dropbox/AD/Figures/pca3.pdf", fonts = "Helvetica", width = 10)
# ggplot(pca3.dat, aes(x = PC1, y = PC4)) + 
#   geom_point(aes(color = group), size = 5) +
#   theme_bw() + 
#   scale_colour_brewer(palette="Set1") +
#   # scale_colour_manual(values = mycolors) +
#   theme(panel.border = element_rect(size = 1, color = "black")) +
#   theme(axis.text = element_text(size = 10, face = "bold"),
#         axis.title = element_text(size = 10, face = "bold")) +
#   theme(legend.position = "top", legend.direction = "horizontal", 
#         legend.text = element_text(size = 10, face = "bold"),
#         legend.title = element_blank(), legend.key = element_blank()) 
# dev.off()
