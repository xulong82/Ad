# --- ADDITIONAL GENE INFORMATION ---
ad.alz <- read.delim("./data/alzgene20120917.tsv", stringsAsFactors = F)$Gene
hg2mus <- read.delim("~/Dropbox/X/hg2mus.map", stringsAsFactors = F, header = F)
ad.alz <- hg2mus$V4[match(ad.alz, hg2mus$V1)]
ad.alz <- unique(ad.alz[!is.na(ad.alz)])
save(ad.alz, file = "data/alzgene.rdt")
ad.kegg <- unique(read.delim("./mmu05010.txt", stringsAsFactors = F, header = F)$V1)
TF <- read.delim("~/Dropbox/X/TFdb.Riken.txt", header = F, stringsAsFactors = F)$V1
