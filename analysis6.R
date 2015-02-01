# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: AD project - Genes & Pathways
# Rev: June 16, 2014

library(ggplot2)
library(gplots)

rm(list = ls())
load("~/Dropbox/AD/R/brain.tpm1.rdt")  # brain TPM without batch correction
load("~/Dropbox/AD/R/brain.tpm2.rdt")  # brain TPM with batch correction
#-------------
tpm1 <- brain.tpm2
id1 <- which(apply(tpm1, 1, function(x) {min(x) < 0}))  # negative TPM from correction
tpm1[id1, ] <- tpm1[id1, ] + abs(apply(tpm1[id1, ], 1, min))
tpm2 <- tpm1

#--- Choose genes by PCA
tpm2.ctr <- tpm2 - apply(tpm2, 1, mean)
tpm2.svd <- svd(tpm2.ctr)
barplot(tpm2.svd$v[, 1])
barplot(tpm2.svd$v[, 2])

type <- gsub("[2456]m.*", "", colnames(tpm2))
cors <- apply(tpm2.svd$v, 2, function (x) {cor(x, as.numeric(as.factor(type)))})

bar3.dat <- data.frame(pc = paste("PC", 1:length(cors), sep = ""), value = abs(cors))
bar3.dat$pc <- factor(bar3.dat$pc, levels = paste("PC", 1:length(cors), sep = ""))

pdf("~/Dropbox/AD/Figures/bar3.pdf", width = 12)
ggplot(bar3.dat, aes(x = pc, y = value)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = .2) +
  theme_bw() + 
  xlab("") + ylab("Pearson's Correlation (abs)") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text.x = element_text(size = 10, angle = -90, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 10, face = "bold"))
dev.off()

id1.pc <- c(3, 4)
id2.pc <- paste("PC", id1.pc, sep = "")
pca4.dat <- data.frame(value = c(tpm2.svd$v[, id1.pc]), 
                       pc = rep(id2.pc, each = ncol(tpm2)),
                       sample = rep(colnames(tpm2), times = length(id1.pc)),
                       type = rep(type, times = length(id1.pc))) 
pca4.dat$order <- reorder(colnames(tpm2), 1:ncol(tpm2))
pdf("~/Dropbox/AD/Figures/pca4.pdf", height = 15, width = 25)
ggplot(pca4.dat, aes(x = order, y = value, fill = type)) + 
  geom_bar(stat = "identity") + facet_grid(pc ~ .) +
  theme_bw() +
  geom_hline(yintercept = 0) +
  xlab("") + ylab("") +
  theme(panel.border = element_rect(size = 2, color = "black")) +
  theme(axis.text.x = element_text(size = 15, angle = -90, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 15, face = "bold"),
        strip.text = element_text(size = 15, face = "bold")) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 15, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank())
dev.off()

tpm2.u <- tpm2.svd$u[, id1.pc]
dimnames(tpm2.u) <- list(rownames(tpm2), id2.pc)
ids <- apply(tpm2.u, 2, function (x) {abs(x) > sd(x)})
gns <- rownames(tpm2.u)[apply(ids, 1, function (x) {sum(as.numeric(x)) > 0})]

#-------------------------------
tpm2 <- log2(tpm1 + 1)  # log2 transformation
write.table(tpm2, file = "~/Dropbox/AD/cluster/data1.txt", sep = "\t", quote = F)
write.table(tpm2[gns, ], file = "~/Dropbox/AD/cluster/data3.txt", sep = "\t", quote = F)

tpm3 <- matrix(nrow = nrow(tpm2), ncol = 4,  # means per month 
               dimnames = list(rownames(tpm2), c("2m", "4m", "5m", "6m")))
for (month in c("2m", "4m", "5m", "6m")) {
  cat(month, "\n")
  dat1 = tpm2[, grep(month, colnames(tpm2))]
  tpm3[, month] = rowMeans(dat1)
}

tpm4 <- tpm2[, grep("5m", colnames(tpm2))]  # 5 month
tpm4 <- tpm4[apply(tpm4, 1, max) > log2(ncol(tpm4) + 1), ]
tpm4.sd <- apply(tpm4, 1, sd)
tpm4 <- tpm4[tpm4.sd > quantile(tpm4.sd, 0.75), ]  # 25% of high SD 

cs1 <- as.character(read.delim("~/Dropbox/AD/List/gene.cs1.txt", header = F)$V1)
cs2 <- as.character(read.delim("~/Dropbox/AD/List/gene.cs2.txt", header = F)$V1)
cs3 <- as.character(read.delim("~/Dropbox/AD/List/gene.cs3.txt", header = F)$V1)
cs4 <- as.character(read.delim("~/Dropbox/AD/List/gene.cs4.txt", header = F)$V1)
cs5 <- as.character(read.delim("~/Dropbox/AD/List/gene.cs5.txt", header = F)$V1)

g2mup <- as.character(read.delim("~/Dropbox/AD/List/gene.05.2m-up.txt", header = F)$V1)
g2mdn <- as.character(read.delim("~/Dropbox/AD/List/gene.05.2m-dn.txt", header = F)$V1)
g4mup <- as.character(read.delim("~/Dropbox/AD/List/gene.05.4m-up.txt", header = F)$V1)
g4mdn <- as.character(read.delim("~/Dropbox/AD/List/gene.05.4m-dn.txt", header = F)$V1)
g5mup <- as.character(read.delim("~/Dropbox/AD/List/gene.05.5m-up.txt", header = F)$V1)
g5mdn <- as.character(read.delim("~/Dropbox/AD/List/gene.05.5m-dn.txt", header = F)$V1)
g6mup <- as.character(read.delim("~/Dropbox/AD/List/gene.05.6m-up.txt", header = F)$V1)
g6mdn <- as.character(read.delim("~/Dropbox/AD/List/gene.05.6m-dn.txt", header = F)$V1)

ids <- unique(c(g2mup, g2mdn, g4mup, g4mdn, g6mup, g6mdn))
ids <- unique(c(g2mup, g2mdn, g4mup, g4mdn, g5mup, g5mdn, g6mup, g6mdn))
dat <- tpm2[ids, ]
dat <- dat[apply(dat, 1, function (x) {max(x) > 3}), ]
dat <- dat[apply(dat, 1, function (x) {max(x) - min(x) > 1}), ]
write.table(dat, file = "~/Dropbox/AD/cluster/data2.txt", sep = "\t", quote = F)
write.table(dat[, grep("APP", colnames(dat))], file = "~/Dropbox/AD/cluster/data4.txt", sep = "\t", quote = F)

hc1 <- hcluster(t(dat), method = "pearson", link = "average")  # clust samples
hc1 <- hcluster(t(dat), method = "euclidean", link = "average")  # clust samples
hc2 <- hcluster(dat, method = "pearson", link = "average")  # clust genes
plot(as.phylo(hc1), type = "unrooted", lab4ut = "axial")
plot(as.phylo(hc1))

css <- c(cs1, cs2, cs3, cs4, cs5)
updns <- unique(c(g2mup, g2mdn, g4mup, g4mdn, g5mup, g5mdn, g6mup, g6mdn))
updns245 <- unique(c(g2mup, g2mdn, g4mup, g4mdn, g5mup, g5mdn))
updns6 <- unique(c(g6mup, g6mdn))
setdiff(updns6, updns245)
setdiff(g6mup, updns245)
setdiff(g6mdn, updns245)

bar3.dat <- tpm2[setdiff(updns6, updns245), grep("6m", colnames(tpm2))]
group <- gsub("[2456]m.*", "", colnames(bar3.dat))

bar3.dat <- data.frame(tpm = c(bar3.dat), 
                       gene = rep(rownames(bar3.dat), time = ncol(bar3.dat)), 
                       group = rep(group, each = nrow(bar3.dat)))
pdf("~/Dropbox/AD/Figures/box3.pdf", width = 12)
ggplot(bar3.dat, aes(x = gene, y = tpm, fill = group)) + 
  geom_boxplot() +
  theme_bw() +
  scale_colour_brewer(palette="Set1") +
  xlab("") + ylab("Log2(TPM)") +
  theme(panel.border = element_rect(size = 2, color = "black")) +
  theme(axis.text = element_text(size = 10, angle = -90, face = "bold"),
        axis.title = element_text(size = 12, face = "bold")) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank()) 
dev.off()

homolog <- read.delim("~/Dropbox/X/homolog.hm.map", header = F)
homolog <- homolog[!duplicated(homolog$V1), ]
homolog <- data.frame(row.names = homolog$V1, mouse = homolog$V4)

tf.rike <- as.character(read.delim("~/Dropbox/X/TFdb.Riken.txt", header = F)$V1)
ad.tops <- as.character(read.delim("~/Dropbox/AD/tophit.csv", header = F)$V1)
ad.gwas <- as.character(read.delim("~/Dropbox/AD/gwas.csv", header = F)$V1)
ad.gwas <- as.character(homolog[ad.gwas, "mouse"])
ad.gwas <- ad.gwas[!is.na(ad.gwas)]

intersect(ad.gwas, ad.tops)
union(ad.gwas, ad.tops)
intersect(union(ad.gwas, ad.tops), tf.rike)

intersect(css, updns)

intersect(css, tf.rike)
intersect(updns, tf.rike)

box2.dat <- tpm2[intersect(updns, union(ad.gwas, ad.tops)), ]
group <- gsub("m.*", "m", colnames(box2.dat))
month <- gsub("WT", "", gsub("APP", "", group))
treat <- gsub("[2456]m", "", group)
box2.dat <- data.frame(tpm = c(box2.dat), 
                       gene = rep(rownames(box2.dat), time = ncol(box2.dat)), 
                       sample = rep(group, each = nrow(box2.dat)),
                       month = rep(month, each = nrow(box2.dat)))
box2.ord <- data.frame(row.names = unique(group)[c(1:2, 4, 3, 6, 5, 8, 7)], order = 1:8)
box2.dat$order <- reorder(box2.dat$sample, box2.ord[as.character(box2.dat$sample), ])
pdf("~/Dropbox/AD/Figures/box2.pdf", width = 12)
ggplot(box2.dat, aes(x = order, y = tpm, fill = sample)) + 
  geom_boxplot() + facet_grid(. ~ gene) +
  theme_bw() +
  scale_colour_brewer(palette="Set1") +
  xlab("") + ylab("Log2(TPM)") +
  theme(panel.border = element_rect(size = 2, color = "black")) +
  theme(axis.text = element_text(size = 10, angle = -90, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 20, face = "bold")) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank()) 
dev.off()

ens.ucsc <- read.delim("~/Dropbox/X/myRefGene.tx")
ens.ucsc <- ens.ucsc[!duplicated(ens.ucsc$name2), ]
ens.ucsc <- data.frame(row.names  = ens.ucsc$name2, 
                       chrom      = ens.ucsc$chrom,
                       chromStart = ens.ucsc$txStart,
                       chromEnd   = ens.ucsc$txEnd)
ens.ucsc$chrom <- as.character(ens.ucsc$chrom)

clust1.dat <- ens.ucsc[updns, ]
clust1.dat <- clust1.dat[!is.na(clust1.dat$chrom), ]
clust1.ord <- data.frame(row.names = c(paste("chr", 1:19, sep = ""), "chrX", "chrY"), order = 1:21)
clust1.dat$order <- reorder(clust1.dat$chrom, clust1.ord[clust1.dat$chrom, ])
pdf(file = "~/Dropbox/AD/Figures/hist3.pdf", width = 12, height = 8)
ggplot(clust1.dat, aes(x = order)) +
  geom_bar() +
  theme_bw() +
  xlab("") + ylab("Count") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text = element_text(size = 10, face = "bold"),
        strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 10, face = "bold"))
dev.off()

cluster.5m.list <- as.character(read.delim("~/Dropbox/AD/Cluster/data4_list.txt", header = F)$V1)
pdf("~/Dropbox/AD/Figures/text3.pdf")
textplot(matrix(cluster.5m.list, ncol = 5), show.rownames = F, show.colnames = F, cex = .5)
dev.off()

id.5m1 <- c("226987.2013", "226989.2013", "226996.2013", "1647.2014", "1738.2014", "2735.2014")
id.5m2 <- c("226954.2013", "1559.2014", "1558.2014", "1636.2014", "1623.2014", "1648.2014", "1633.2014", "2751.2014")
app.5m <- tpm2[, grep("APP5m", colnames(tpm2))]
colnames(app.5m) <- gsub("APP5m", "", colnames(app.5m))
app.5m <- app.5m[, c(id.5m1, id.5m2)]
treat <- c(rep("CS1", length(id.5m1)), rep("CS2", length(id.5m2)))

tt.pval <- apply(app.5m, 1, function(x) {pairwise.t.test(x, treat, p.adj = "none")$p.value})
tt.qval <- p.adjust(tt.pval, method = "fdr")
id.gene <- names(tt.qval[which(tt.qval < .05)])
data.c1 <- app.5m[id.gene, id.5m1]
data.c2 <- app.5m[id.gene, id.5m2]
mean.c1 <- 2^(rowMeans(data.c1))
mean.c2 <- 2^(rowMeans(data.c2))
log2.fc <- mean.c2 / mean.c1

howell <- data.frame(row.names = id.gene, mean.c1, mean.c2, log2.fc, qval = tt.qval[id.gene])

write.table(howell, file = "~/Dropbox/AD/app5m.tpm.txt", sep = "\t", row.names = T, col.names = T, quote = F)
write.table(id.gene, file = "~/Dropbox/AD/app5m.txt", sep = "\t", row.names = F, col.names = F, quote = F)
#--- This goes to DAVID ---
