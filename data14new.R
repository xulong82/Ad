# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: AD project - import the data from HPC
# Rev: Dec 12, 2014

rm(list = ls())

#--- 2014 new brain ---
name1 <- list.files(path = "~/Dropbox/AD/RSEM/howell_new/brain", pattern = "*.genes.results")
name2 <- paste("mouse_", gsub("_GES14_.*", "", name1), sep = "")
for (i in 1:length(name1)) {
  filepath <- file.path("/data/xwang//AD/RSEM/howell_new/brain", name1[i])
  cat(i, "/", length(name1), name1[i], "\n")
  assign(name2[i], read.delim(filepath, stringsAsFactors = F)[, -2])
}

#--- 2014 new retina ---
name3 <- list.files(path = "/data/xwang//AD/RSEM/howell_new/retina", pattern = "*.genes.results")
name4 <- paste("mouse_", gsub("_GES14_.*", "", name3), sep = "")
for (i in 1:length(name3)) {
  filepath <- file.path("/data/xwang//AD/RSEM/howell_new/retina", name3[i])
  cat(i, "/", length(name3), name3[i], "\n")
  assign(name4[i], read.delim(filepath, stringsAsFactors = F)[, -2])
}

#-------------
map <- read.delim("~/Dropbox/AD/map3.csv", header = F, stringsAsFactors = F)
map$V1 <- gsub(" ", "_", map$V1)
map$V2 <- paste("mouse", map$V2, sep = "_")
map$symbol <- paste(map$V2, map$V1, sep = "_")

symbol_brain <- map[match(name2, map$V2), "symbol"]

name5 <- gsub("mouse_3551R", "mouse_3351R", name4)
symbol_retina <- map[match(gsub("R", "", name5), map$V2), "symbol"]

# --- 3551R_GES14_04146_GTCCGC.genes.results

brain14new.tpm <- matrix(nrow = nrow(get(name2[1])), ncol = length(name2), 
  dimnames = list(get(name2[1])$gene_id, symbol_brain))

for (i in 1:length(name2)) {
  brain14new.tpm[, i] <- get(name2[i])$TPM
}

retina14new.tpm <- matrix(nrow = nrow(get(name4[1])), ncol = length(name4), 
  dimnames = list(get(name4[1])$gene_id, symbol_retina))

for (i in 1:length(name4)) {
  retina14new.tpm[, i] <- get(name4[i])$TPM
}

#-------------
ens2sym <- read.delim("~/Dropbox/X/ensembl2symbol.map", header = F, row.names = 1, stringsAsFactors = F)

brain14new.tpm <- as.data.frame(brain14new.tpm)
brain14new.tpm$symbol <- ens2sym[rownames(brain14new.tpm), ]
brain14new.tpm <- aggregate(. ~ symbol, brain14new.tpm, sum)

retina14new.tpm <- as.data.frame(retina14new.tpm)
retina14new.tpm$symbol <- ens2sym[rownames(retina14new.tpm), ]
retina14new.tpm <- aggregate(. ~ symbol, retina14new.tpm, sum)

save(brain14new.tpm, retina14new.tpm, file = "~/Dropbox/AD/R/data14new.rdt")

#-------------
load("~/Dropbox/AD/R/data14new.rdt")  # tpm per sample
heatmap(cor(brain14new.tpm[, -1]))
