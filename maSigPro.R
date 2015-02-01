# --- maSigPro model

library(maSigPro)
library(xtable)

rm(list = ls())

setwd("~/Dropbox/AD")
load("./R/batch2014.rdt")
load("./R/cluster1.rdt")

# dt <- dt[apply(dt, 1, function(x) max(x) - min(x) > 1), ]

treat <- gsub("^.*(WT|APP).*", "\\1", colnames(dt))
month <- gsub("^.*(2m|4m|5m|6m).*", "\\1", colnames(dt))
uid <- paste(treat, month, sep = "")

edesign <- cbind(Time = as.numeric(gsub("m", "", month)), Replicate = as.numeric(as.factor(uid)), 
                 Control = as.numeric(treat == "WT"), MT = as.numeric(treat == "APP"))
rownames(edesign) <- colnames(dt)
design <- make.design.matrix(edesign, degree = 3)

fit <- p.vector(dt, design, Q = 0.05, MT.adjust = "BH")
tstep <- T.fit(fit, step.method = "backward", alfa = 0.05)
tstep$coefficients

sigs <- get.siggenes(tstep, rsq = 0.5, vars = "groups")
sigs <- get.siggenes(tstep, rsq = 0.5, vars = "all")
sigs <- get.siggenes(tstep, rsq = 0.5, vars = "each")

see.genes(sigs$sig.genes$MTvsControl, main = "MTvsControl", show.fit = T, dis = design$dis, 
          cluster.method = "kmeans", cluster.data = 1, k = 3)

print(xtable(sigs$summary[, -1]))

idx <- unique(c(as.matrix(sigs$summary)))
idx <- idx[idx != " "]
dt <- dt[idx, ]

hc1 <- hcluster(t(dt), method = "pearson", link = "average")
hc2 <- hcluster(dt, method = "pearson", link = "average")

hc1 <- hcluster(t(dt), method = "correlation", link = "centroid")
hc2 <- hcluster(dt, method = "correlation", link = "centroid")

plot(as.phylo(hc1), type = "unrooted", lab4ut = "axial")

pdf("./Graphs/phylo_sample.pdf", fonts = "Helvetica", height = 4, width = 10)
plot(as.phylo(hc1), edge.width = 2, font = 2, cex = 0.7, label.offset = 1e-2, direction = "downward")
dev.off()

pdf("./Graphs/phylo_gene.pdf", fonts = "Helvetica", width = 3, height = 9)
col.manual <- c("grey70", "firebrick1", "chartreuse3", "dodgerblue3", "gold1", "darkorchid2")
clusts = cutree(hc2, 6)
plot(as.phylo(hc2), tip.color = col.manual[clusts], cex = .5, 
     edge.width = 2, adj = 0.5, font = 2, label.offset = 1e-2, direction = "rightward")
dev.off()

tile.dt1 <- dt[hc2$order, hc1$order]
tile.dt2 <- t(apply(tile.dt1, 1, scale))
colnames(tile.dt2) <- colnames(tile.dt1)
tile.dt3 <- data.frame(value = c(tile.dt2), 
                       gene = rep(rownames(tile.dt2), time = ncol(tile.dt2)),
                       sample = rep(colnames(tile.dt2), each = nrow(tile.dt2)))
tile.dt3$sample <- factor(tile.dt3$sample, levels = colnames(tile.dt1))
tile.dt3$gene <- factor(tile.dt3$gene, levels = rownames(tile.dt1))

pdf("./Graphs/heatmap.pdf", width = 10, height = 6)
ggplot(tile.dt3, aes(x = sample, y = gene, fill = value, alpha = abs(value))) + 
  geom_tile() + guides(alpha = F) +
  scale_fill_gradient(low="blue", high="red") +
  theme_bw() + xlab("") + ylab("") + 
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 8), axis.ticks = element_blank()) +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.margin = unit(c(0, 0, -6, 0), "cm"), 
        legend.text = element_text(size = 10, face = "bold"), legend.title = element_blank())
dev.off()


suma2Venn(sigs$summary[, 1:2])
PlotGroups(data.abiotic["STMCL34", ], edesign = edesign.abiotic)
see.genes(sigs$sig.genes$ColdvsControl, main = "ColdvsControl", show.fit = T, dis = design$dis, 
          cluster.method = "kmeans", cluster.data = 1, k = 9)

# ---
myfit <- lm(y ~ ., as.data.frame(dt))
myfit <- stepAIC(myfit, direction = "backward")
