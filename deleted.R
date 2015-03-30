myfitAic <- list()
dt.aic <- dt[which(qval < 0.05), ]
for (i in 1:nrow(dt.aic)) {
  if (i %% 1e2 == 0) cat(i, "\n")
  y <- dt.aic[i, ]
  fitLr <- lm(y ~ age + group + batch + age * group)
  fitAic <- stepAIC(fitLr, direction = "backward", trace = 0)
  myfitAic[[i]] <- summary(fitAic)
}

r2 <- sapply(myfitAic, function (x) x$r.squared)
fval <- sapply(myfitAic, function (x) x$fstatistic)
pval <- apply(fval, 2, function (x) pf(x[1], x[2], x[3], lower.tail = F))
qval <- p.adjust(pval, method = "fdr")

cond <- c("age4m", "age5m", "age6m", "groupAPP", "batchmouse", "age4m:groupAPP", "age5m:groupAPP", "age6m:groupAPP")
grouping <- matrix(FALSE, nrow = length(geneId), ncol = length(cond), dimnames = list(geneId, cond))
for (i in 1:length(geneId)) {
  if (i %% 1e3 == 0) cat(i, "\n")
  fit0 <- myfitAic[[i]]
  fit0.r2 <- fit0$r.squared
  fit0.coef <- fit0$coefficients
  cond0 <- rownames(fit0.coef)[fit0.coef[, "Pr(>|t|)"] < 0.01]
  # grouping[i, ] <- cond %in% cond0
  grouping[i, ] <- cond %in% cond0 & fit0.r2 > 0.5
}

grouping.geneId <- list()
for (idx in colnames(grouping)) grouping.geneId[[idx]] <- c(idx, rownames(grouping)[grouping[, idx]])

number <- sapply(grouping.geneId, length)
col.manual <- c("grey70", "firebrick1", "chartreuse3", "dodgerblue3", "gold1", "darkorchid2")
col.manual <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf")
op <- par(mar = c(15, 4, 4, 2))
bar <- barplot(number, las = 3, col = col.manual, border = NA, ylim = c(0, 300), axes = F)
abline(0, 0, lwd = 5, col = "black")
text(x = bar, y = number + 5, labels = number, font = 2, col = "grey20")

lapply(grouping.geneId, write, "./DE/glm.txt", append = TRUE, ncolumns = 1000)

geneId.de = geneId[as.logical(apply(grouping[, -5], 1, max))]
write.table(geneId.de, file = "./DE/geneId.txt", sep = "\t", row.names = F, col.names = F, quote = F)

# --- optional: 2014 samples only
brain.tpm <- brain.tpm[, c(grep("2014", colnames(brain.tpm)), grep("mouse", colnames(brain.tpm)))]

# --- PCA
brain.tpm_ts <- log2(brain.tpm + 1)
brain.tpm_ts <- brain.tpm_ts - apply(brain.tpm_ts, 1, mean)
pca_brain <- prcomp(brain.tpm_ts)

barplot(pca_brain$rotation[, 1])

genotype <- rep("WT", ncol(brain.tpm))
genotype[grep("APP", colnames(brain.tpm))] <- "APP"
graph.dt <- data.frame(value = c(pca_brain$rotation[, 1:3]), 
                       pc = rep(paste("PC", 1:3, sep = ""), each = ncol(brain.tpm)), 
                       sample = rep(colnames(brain.tpm), 3),
                       type = rep(genotype, 3))
graph.dt$order <- reorder(colnames(brain.tpm), 1:ncol(brain.tpm))

# col.manual <- c("grey70", "firebrick1", "chartreuse3", "dodgerblue3", "gold1", "darkorchid2")
ggplot(graph.dt, aes(x = order, y = value, fill = type)) + 
  geom_bar(stat = "identity", width = .75) + facet_grid(pc ~ .) +
  scale_fill_manual(values = c("firebrick1", "chartreuse3")) +
  theme_bw() +
  geom_hline(yintercept = 0) +
  xlab("") + ylab("") +
  theme(panel.border = element_rect(size = 1, color = "grey30")) +
  theme(axis.text.x = element_text(size = 6, angle = -90, face = "bold", color = "grey30"),
        axis.text.y = element_text(size = 6, face = "bold"),
        axis.title = element_text(size = 8, face = "bold"),
        strip.text = element_text(size = 8, face = "bold")) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 8, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank())

#--- Choose genes by PCA
tpm2.ctr <- tpm2 - apply(tpm2, 1, mean)
tpm2.svd <- svd(tpm2.ctr)
barplot(tpm2.svd$v[, 1])
barplot(tpm2.svd$v[, 2])

type <- gsub("[2456]m.*", "", colnames(tpm2))
cors <- apply(tpm2.svd$v, 2, function (x) {cor(x, as.numeric(as.factor(type)))})

bar3.dat <- data.frame(pc = paste("PC", 1:length(cors), sep = ""), value = abs(cors))
bar3.dat$pc <- factor(bar3.dat$pc, levels = paste("PC", 1:length(cors), sep = ""))

ggplot(bar3.dat, aes(x = pc, y = value)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = .2) +
  theme_bw() + 
  xlab("") + ylab("Pearson's Correlation (abs)") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text.x = element_text(size = 10, angle = -90, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 10, face = "bold"))

tpm2.u <- tpm2.svd$u[, id1.pc]
dimnames(tpm2.u) <- list(rownames(tpm2), id2.pc)
ids <- apply(tpm2.u, 2, function (x) {abs(x) > sd(x)})
gns <- rownames(tpm2.u)[apply(ids, 1, function (x) {sum(as.numeric(x)) > 0})]

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
ggplot(clust1.dat, aes(x = order)) +
  geom_bar() +
  theme_bw() +
  xlab("") + ylab("Count") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text = element_text(size = 10, face = "bold"),
        strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 10, face = "bold"))
