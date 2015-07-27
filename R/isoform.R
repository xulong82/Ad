rm(list = ls())

# 2014 new samples
fname <- list.files(path = "/data/xwang/AD/RSEM/howell_new/brain", pattern = "*.isoforms.results")
brain2014new <- lapply(fname, function(x) { cat(x, "\n")
  filepath = file.path("/data/xwang/AD/RSEM/howell_new/brain", x)
  read.delim(filepath, stringsAsFactors = F)[, -2]
}); names(brain2014new) <- paste0("mouse_", gsub("_GES14_.*", "", fname))

brain2014new_tpm <- sapply(brain2014new, function(x) x$TPM)
rownames(brain2014new_tpm) <- brain2014new[[1]]$transcript_id

map2014new <- read.delim("~/Dropbox/GitHub/Ad/Samples/map3.csv", header = F, stringsAsFactors = F)
map2014new$V1 <- gsub(" ", "_", map2014new$V1)
map2014new$V2 <- paste("mouse", map2014new$V2, sep = "_")
map2014new$symbol <- paste(map2014new$V2, map2014new$V1, sep = "_")
colnames(brain2014new_tpm) <- map2014new$symbol[match(colnames(brain2014new_tpm), map2014new$V2)]

# 2014 first 
fname <- list.files(path = "/data/xwang/AD/RSEM/howell_2014/brain", pattern = "*.isoforms.results")
brain2014 <- lapply(fname, function(x) { cat(x, "\n")
  filepath = file.path("/data/xwang/AD/RSEM/howell_2014/brain", x)
  read.delim(filepath, stringsAsFactors = F)[, -2]
}); names(brain2014) <- paste0("M", gsub("_GES14_.*", "", fname))

brain2014_tpm <- sapply(brain2014, function(x) x$TPM)
rownames(brain2014_tpm) <- brain2014[[1]]$transcript_id

map2014 <- read.delim("~/Dropbox/GitHub/Ad/Samples/map1.csv", header = F, stringsAsFactors = F, sep = " ")
map2014$symbol1 <- paste0("M", map2014$V1)
map2014$symbol2 <- paste0(paste0(map2014$V6, map2014$V5), map2014$V1)
colnames(brain2014_tpm) <- map2014$symbol2[match(colnames(brain2014_tpm), map2014$symbol1)]
colnames(brain2014_tpm) <- paste(colnames(brain2014_tpm), "2014", sep = ".")

all(rownames(brain2014new_tpm) == rownames(brain2014_tpm))
brain2014All <- cbind(brain2014_tpm, brain2014new_tpm)

setwd("~/Dropbox/GitHub/Ad")
save(brain2014All, file = "data/brain2014All_isoform.rdt")

