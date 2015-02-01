# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: AD project: B6xC3H
# Rev: July 7, 2014

library(ggplot2)

rm(list = ls())
#--- 2014 brain ---
name1 <- list.files(path = "~/Dropbox/AD/B6xC3H/emase2",
                    pattern = "*.emase.groups")
name2 <- paste("M", gsub("_GES14_.*", "", name1), sep = "")

ens2sym <- read.delim("~/Dropbox/X/ensembl2symbol.map", header = F, row.names = 1)
ens2sym$V2 <- as.character(ens2sym$V2)
for (i in 1:length(name1)) {
  filepath <- file.path("~/Dropbox/AD/B6xC3H/emase2", name1[i])
  cat(i, "/", length(name1), name1[i], "\n")
  temp <- read.delim(filepath)
  temp$locus <- as.character(temp$locus)
  temp$symbol <- ens2sym[temp$locus, ]
  temp <- temp[!duplicated(temp$symbol), ]
  # temp <- aggregate(. ~ symbol, temp, sum)
  temp <- temp[!is.na(temp$symbol), ]
  temp <- data.frame(row.names = temp$symbol, B6 = temp$B6, C3H = temp$C3H)
  assign(name2[i], temp)
}

map1 <- read.delim("~/Dropbox/AD/map1.csv", sep = " ", header = F)[, c(1, 6, 5)]
map1$sym <- paste(paste(map1$V6, map1$V5, sep = ""), map1$V1, sep = "")
map1 <- data.frame(row.names = paste("M", map1$V1, sep = ""), sym = map1$sym)
map1$sym <- as.character(map1$sym)
dt.emase <- list()
for (i in 1:length(name1)) {
  dt.emase[[i]] <- get(name2[i])
}
names(dt.emase) <- map1[name2, ]
save(dt.emase, file = "~/Dropbox/AD/R/emase2.rdt")

dt.emase.wt<- dt.emase[grep("WT", names(dt.emase))]
dt.emase.app <- dt.emase[grep("APP", names(dt.emase))]

dt.emase.wt.c3h <- lapply(dt.emase.wt, function (x) {y1 = x[x$C3H > x$B6, ]
                                                     y2 = y1[rowSums(y1) > 10, ]
                                                     y3 = y2[(y2$C3H + 1) / (y2$B6 + 1) > 1.5, ]
                                                     return(y3)})
dt.emase.app.c3h <- lapply(dt.emase.app, function (x) {y1 = x[x$C3H > x$B6, ]
                                                       y2 = y1[rowSums(y1) > 10, ]
                                                       y3 = y2[(y2$C3H + 1) / (y2$B6 + 1) > 1.5, ]
                                                       return(y3)})

chr.ucsc <- read.delim("~/Dropbox/X/chromInfo.txt", header = F)  # chromosomes length
chrlen <- cumsum(as.numeric(chr.ucsc$V2)) * 1e-6
names(chrlen) <- 1:21
chrmid <- diff(c(0, chrlen)) * 0.5 + c(0, chrlen[-length(chrlen)])
ens.ucsc <- read.delim("~/Dropbox/X/myRefGene.tx")
ens.ucsc <- ens.ucsc[!duplicated(ens.ucsc$name2), ]
ens.ucsc <- data.frame(row.names  = ens.ucsc$name2, 
                       chrom      = ens.ucsc$chrom,
                       chromStart = ens.ucsc$txStart,
                       chromEnd   = ens.ucsc$txEnd)
ens.ucsc$chrom <- as.character(ens.ucsc$chrom)
ens.ucsc$chrom <- gsub("chrX", "chr20", ens.ucsc$chrom)
ens.ucsc$chrom <- gsub("chrY", "chr21", ens.ucsc$chrom)
ens.ucsc$chrom <- gsub("chr", "", ens.ucsc$chrom)
ens.ucsc$chrom <- as.numeric(ens.ucsc$chrom)
chrlen1 <- c(0, chrlen)
ens.ucsc.pos <- chrlen1[ens.ucsc$chrom] + rowMeans(ens.ucsc[, 2:3]) * 1e-6
ens.ucsc <- cbind(ens.ucsc, pos = ens.ucsc.pos)

dt.emase.c3h <- lapply(dt.emase, function (x) {y1 = x[x$C3H > x$B6, ]
                                               y2 = y1[rowSums(y1) > 10, ]
                                               y3 = y2[(y2$C3H + 1) / (y2$B6 + 1) > 1.5, ]
                                               return(y3)})
map1.emase <- data.frame(value = dt.emase.c3h[[1]]$C3H, 
                         gene = rownames(dt.emase.c3h[[1]]),
                         sample = rep(names(dt.emase.c3h)[1], nrow(dt.emase.c3h[[1]])))
for (i in 2:length(dt.emase.c3h)) {
  temp <- data.frame(value = dt.emase.c3h[[i]]$C3H, 
                     gene = rownames(dt.emase.c3h[[i]]),
                     sample = rep(names(dt.emase.c3h)[i], nrow(dt.emase.c3h[[i]])))
  map1.emase <- rbind(map1.emase, temp)
}
map1.emase$gene <- as.character(map1.emase$gene)
map1.emase$sample <- as.character(map1.emase$sample)
map1.emase$type <- gsub("[2456]m.*", "", map1.emase$sample)
map1.emase$pos <- ens.ucsc[map1.emase$gene, "pos"]
map1.emase <- map1.emase[!is.na(map1.emase$pos), ]
                         
map2.emase <- map1.emase[map1.emase$type == "WT", ]
map3.emase <- map1.emase[map1.emase$type == "APP", ]

pdf("~/Dropbox/AD/B6xC3H/map2.pdf", width = 10, height = 10)  # WT
ggplot(map1.emase, aes(x = pos, y = value, color = type)) +
  geom_line() +
  facet_grid(sample ~ .) +
  theme_bw() + 
  xlab("") + ylab("Depth") +
  scale_x_continuous(breaks = chrmid, labels = names(chrlen)) +
# scale_y_continuous(limits=c(0, 8)) +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text = element_text(size = 6, face = "bold"),
        strip.text = element_text(size = 0, face = "bold"),
        axis.title = element_text(size = 5, face = "bold")) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 5, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank())
dev.off()
