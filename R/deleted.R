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

# Mixed Age and App effects
mixed_eff <- lapply(1:ncol(combn), function(x1) {
  y1 = sapply(glm.fit, function(x2) {
    y2 = contrast(x2, list(age = combn[2, x1], group = "APP"), list (age = combn[1, x1], group = "APP")) 
    c(y2$Contrast, y2$Pvalue)})
  rownames(y1) = c("contrast", "pval"); t(y1)
}); names(mixed_eff) <- apply(combn, 2, function(x) paste0(x[2], x[1])) 

fit.coef <- lapply(glm.fit, function(x) summary(x)$coefficients)
fit.r2 <- sapply(glm.fit, function (x) summary(x)$r.squared)
fit.qval <- sapply(glm.fit, function (x) summary(x)$fstatistic) %>% 
  apply(2, function (x) pf(x[1], x[2], x[3], lower.tail = F)) %>% p.adjust(method = "fdr")

fit.qr <- glm.fit[fit.qval < 0.05 & fit.r2 > 0.5]
fit.qr <- glm.fit[fit.qval < 0.05 & fit.r2 > 0.4]

fit.est <- lapply(fit.qr, function (x) x$coefficients[, "Estimate"])
fit.est <- do.call(rbind, fit.est)
fit.pval <- lapply(fit.qr, function (x) x$coefficients[, "Pr(>|t|)"])
fit.pval <- do.call(rbind, fit.pval)
  
# --- GLM cohorts ---
geneId <- names(fit.qr)
logit <- apply(fit.pval, 2, function (x) x < 0.05) & apply(fit.est, 2, function (x) abs(x) > 0.5)
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
  
# --- GRAPH
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
# scale_x_discrete(labels = c("2m:APP", "4m:APP", "5m:APP", "6m:APP")) +
  scale_x_discrete(labels = c("2m:APP", "5m:APP", "6m:APP")) +
# scale_x_discrete(labels = c("4m:WT", "5m:WT", "6m:WT", "4m:APP", "5m:APP", "6m:APP")) +
# scale_x_discrete(labels = c("5m:WT", "6m:WT", "5m:APP", "6m:APP")) +
  scale_y_discrete(labels = profile.table) +
  scale_fill_manual(values = c("grey80", "firebrick1")) 

# --- BUILD THE RESULT OBJECT
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

theId <- theGK <- list()
theId[[1]] <- glm$more$app$symbol_by_pattern
theId[[2]] <- glm$more$age$symbol_by_pattern
theId[[3]] <- glm$less$app$symbol_by_pattern
theId[[4]] <- glm$less$age$symbol_by_pattern

for (i in 1:length(theId)) theGK[[i]] <- lapply(theId[[i]], myGK)

glm$more$app$gk_by_pattern <- theGK[[1]]
glm$more$age$gk_by_pattern <- theGK[[2]]
glm$less$app$gk_by_pattern <- theGK[[3]]
glm$less$age$gk_by_pattern <- theGK[[4]]

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
require(HiveR)
require(plyr)
require(colorspace)
require(classInt)

d = ggplot2::diamonds
d = d[,c(1:4,7)]
head(d); dim(d)

# separate carat-size data into equal interval groups
brks = classIntervals(d$carat, n=11, style="quantile")$brks[1:11] # also try 'equal' style
d$carat = findInterval(d$carat, brks)

## NODES DATA

nodegroups = list()
for(i in 1:4){
  vals = as.numeric(unique(d[[i]]))
  nodegroup = data.frame(id = 1:length(vals), lab = unique(d[[i]]), vals = vals,
                         radius = 100 * vals/max(vals), axis = i)
  sizes = table(d[[i]])
  nodegroup$size = as.numeric(sizes[ match(nodegroup$lab, names(sizes)) ])
  nodegroup$size = 2 * nodegroup$size / max(nodegroup$size)
  if(i>1) nodegroup$id = nodegroup$id + max(nodegroups[[i-1]]$id)
  nodegroups[[ names(d)[i] ]] = nodegroup
}
nodegroups

nodes = rbind(nodegroups[[1]], nodegroups[[2]], nodegroups[[3]], nodegroups[[4]])
nodes$lab = as.character(nodes$lab)
nodes$axis = as.integer(nodes$axis)
nodes$radius = as.numeric(nodes$radius)
nodes$color = "#ffffff"
head(nodes)

## EDGES DATA

# first update edge data with new node IDs
head(d)
for(i in 1:4) {
  header = paste0(names(nodegroups)[i], 'id')
  d[[header]] = nodegroups[[i]]$id[ match(as.numeric(d[[i]]), nodegroups[[i]]$vals) ]
}
head(d)

# edges between the 4 axes in terms of node IDs
for(i in 6:8){
  edgegroup = data.frame(id1 = d[[i]], id2 = d[[i+1]], price = d[[5]])
  if(i==6) all_edges = edgegroup else all_edges = rbind(all_edges, edgegroup)
}
head(all_edges); dim(all_edges)

# summarise edge data
edges = aggregate(all_edges$price, by=list(all_edges$id1, all_edges$id2), FUN='mean')
names(edges) = c('id1','id2','price')
edges = edges[with(edges, order(id1,id2)),] # reorder

# set edge weights (stroke thickness)
weights = count(all_edges, vars = c('id1', 'id2')) # summary data
weights = weights[with(weights, order(id1,id2)),] # reorder to match egdes
all(weights$id1 == edges$id1, weights$id2 == edges$id2) # check all IDs match up
edges$weight = weights$freq * 0.004
edges$weight = pmax(edges$weight, 0.2) # set min edge weight to still visible
range(weights$freq)
range(edges$weight)

# normalise prices for each group of edges (to utilise full colour range)
p = edges$price
edges$colorvals = 0

for(i in nodegroups[1:3]){
  sel = edges$id1 %in% range(i$id)[1] : range(i$id)[2]
  edges$colorvals[sel] = (p[sel] - min(p[sel])) / (max(p[sel]) - min(p[sel]))
}

edges$color = paste0(hex(HSV(edges$colorvals * 300, 1, 1)), '60') # set alpha
edges = edges[order(edges$weight, decreasing=T),] # draw thin edges last

head(edges)

hpd = list()
hpd$nodes = nodes
hpd$edges = edges
hpd$type = "2D"
hpd$desc = "Diamonds"
hpd$axis.cols = rep('#00000000', 4) # make invisible
hpd$axLabs = c("carats","cut","colour","clarity")
class(hpd) = "HivePlotData"

# Check data correctly formatted
chkHPD(hpd, confirm = TRUE)

plotHive(hpd, axLabs = hpd$axLabs, ch = 0.1)
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
