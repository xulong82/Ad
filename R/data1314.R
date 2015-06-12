# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: AD project - import the data from HPC
# Rev: June 2, 2014

rm(list = ls())

#--- 2014 brain ---
name1 <- list.files(path = "/data/xwang/AD/RSEM/howell_2014/brain",
                    pattern = "*.genes.results")
name2 <- paste("M", gsub("_GES14_.*", "", name1), sep = "")
for (i in 1:length(name1)) {
  filepath <- file.path("/data/xwang/AD/RSEM/howell_2014/brain", name1[i])
  cat(i, "/", length(name1), name1[i], "\n")
  assign(name2[i], read.delim(filepath)[, -2])
}

#--- 2014 retina ---
name3 <- list.files(path = "/data/xwang/AD/RSEM/howell_2014/retina",
                    pattern = "*.genes.results")
name4 <- paste("M", gsub("_GES14_.*", "", name3), sep = "")
for (i in 1:length(name3)) {
  filepath <- file.path("/data/xwang/AD/RSEM/howell_2014/retina", name3[i])
  cat(i, "/", length(name3), name3[i], "\n")
  assign(name4[i], read.delim(filepath)[, -2])
}

#--- 2013 brain ---
name5 <- list.files(path = "/data/xwang/AD/RSEM/howell_2013",
                    pattern = "*.genes.results")
name6 <- gsub(".genes.results", "", name5)
for (i in 1:length(name5)) {
  filepath <- file.path("/data/xwang/AD/RSEM/howell_2013", name5[i])
  cat(i, "/", length(name5), name5[i], "\n")
  assign(name6[i], read.delim(filepath)[, -2])
}

#-------------
map1 <- read.delim("~/Dropbox/AD/map1.csv", sep = " ", header = F)[, c(1, 6, 5)]
map1$sym <- paste(paste(map1$V6, map1$V5, sep = ""), map1$V1, sep = "")
map1 <- data.frame(row.names = paste("M", map1$V1, sep = ""), sym = map1$sym)

brain2014.tpm <- matrix(nrow = nrow(get(name2[1])), ncol = length(name2), 
                        dimnames = list(get(name2[1])$gene_id, map1[name2, ]))

for (i in 1:length(name2)) {
  brain2014.tpm[, i] <- get(name2[i])$TPM
}

map2 <- read.delim("~/Dropbox/AD/map2.csv", sep = " ", header = F)
map2$sym <- paste(paste(paste(map2$V2, map2$V1, sep = ""), map2$V3, sep = ""), map2$V4, sep = "")
map2 <- data.frame(row.names = paste(paste("M", map2$V3, sep = ""), map2$V4, sep = ""), sym = map2$sym)

retina2014.tpm <- matrix(nrow = nrow(get(name4[1])), ncol = length(name4), 
                         dimnames = list(get(name4[1])$gene_id, map2[name4, ]))

for (i in 1:length(name4)) {
  retina2014.tpm[, i] <- get(name4[i])$TPM
}

brain2013.tpm <- matrix(nrow = nrow(get(name6[1])), ncol = length(name6), 
                        dimnames = list(get(name6[1])$gene_id, name6))

for (i in 1:length(name6)) {
  brain2013.tpm[, i] <- get(name6[i])$TPM
}

# for (i in 1:length(name3)) {
#   name.tmp1 <- name3[i]  # mouse id
#   name.tmp2 <- name2[grep(name.tmp1, name2)]  # data objects
#   cat(name.tmp1, ":", name.tmp2, "\n")
#   if (length(name.tmp2) == 1) {
#     brain2014.tpm[, i] <- get(name.tmp2)$TPM
#     brain2014.cou[, i] <- get(name.tmp2)$expected_count
#   } else {
#     tpm.tmp1 <- get(name.tmp2[1])$TPM
#     tpm.tmp2 <- get(name.tmp2[1])$expected_count
#     for (j in 2:length(name.tmp2)) {
#       tpm.tmp1 <- tpm.tmp1 + get(name.tmp2[j])$TPM
#       tpm.tmp2 <- tpm.tmp2 + get(name.tmp2[j])$expected_count
#     }
#     brain2014.tpm[, i] <- tpm.tmp1 / length(name.tmp2)  # mean of TPM in replicates
#     brain2014.cou[, i] <- tpm.tmp2 / length(name.tmp2)  # mean of TPM in replicates
#   }
# }

#-------------
ens2sym <- read.delim("~/Dropbox/X/ensembl2symbol.map", header = F, row.names = 1)
brain2013.tpm <- as.data.frame(brain2013.tpm)
brain2014.tpm <- as.data.frame(brain2014.tpm)
retina2014.tpm <- as.data.frame(retina2014.tpm)
brain2013.tpm$symbol <- ens2sym[rownames(brain2013.tpm), ]
brain2014.tpm$symbol <- ens2sym[rownames(brain2014.tpm), ]
retina2014.tpm$symbol <- ens2sym[rownames(retina2014.tpm), ]
brain2013.tpm <- aggregate(. ~ symbol, brain2013.tpm, sum)
brain2014.tpm <- aggregate(. ~ symbol, brain2014.tpm, sum)
retina2014.tpm <- aggregate(. ~ symbol, retina2014.tpm, sum)

colnames(brain2013.tpm) <- gsub("Het", "APP", colnames(brain2013.tpm))
colnames(brain2013.tpm)[-1] <- paste(colnames(brain2013.tpm)[-1], "2013", sep = ".")
colnames(brain2014.tpm)[-1] <- paste(colnames(brain2014.tpm)[-1], "2014", sep = ".")
colnames(retina2014.tpm)[-1] <- paste(colnames(retina2014.tpm)[-1], "2014", sep = ".")

save(brain2013.tpm, brain2014.tpm, retina2014.tpm, file = "~/Dropbox/AD/R/data1314.rdt")
load("~/Dropbox/AD/R/data1314.rdt")  # tpm per sample

write.table(retina2014.tpm, file = "~/Dropbox/AD/file1.txt")

#-------------
heatmap(as.matrix(brain[, -1]))
