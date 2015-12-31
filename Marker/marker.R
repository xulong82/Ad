library(dplyr)
library(rstan)  # R/3.1.1
library(pheatmap)

rm(list = ls())
cadillac <- "/data/xwang/Load"
github <- "~/GitHub/Load"

setwd(cadillac)  # ------------------------------------
load("GLM/fit_sampling.rdt")  # GLM fit with Stan

grp <- c("B6", "ApoE4", "Apoe", "Bin1", "Cd2ap", "Clu") 

# Gene expression as GLM estimations's mode
mode <- sapply(fit, function(x) {
  as.data.frame(x) %>% select(contains("beta")) %>%
  apply(2, function(y) {  # mode
    dens <- density(y); dens$x[which.max(dens$y)]
  })  # each cell's mode estimate of a single gene
}) %>% t
colnames(mode) <- grp 

# 95 credible interval: low end
summary <- lapply(fit, function(x) summary(x)$summary)
ci95_lo <- sapply(summary, function(x) x[, "2.5%"]) %>% t %>% as.data.frame %>% select(contains("beta"))
ci95_hi <- sapply(summary, function(x) x[, "97.5%"]) %>% t %>% as.data.frame %>% select(contains("beta"))
colnames(ci95_lo) <- colnames(ci95_hi) <- grp

# Markers
profile0 <- sapply(grp, function(x) rownames(ci95_lo)[ci95_lo[, x] > 6])
profile1 <- lapply(profile0, function(x) x[x%in% names(which(table(do.call(c, profile0)) < 3))])
profile2 <- lapply(profile0, function(x) x[x%in% names(which(table(do.call(c, profile0)) < 6))])

lapply(profile1, write, "Marker/1.txt", append = TRUE, ncolumns = 1e3)
lapply(profile2, write, "Marker/2.txt", append = TRUE, ncolumns = 1e3)

setwd(github)  # ------------------------------------

# Enrichment
load("../Brain/shiny/shinyList.rdt")  # background
bg <- shinyList$bg
marker1 <- shinyList$marker1  # less than 47
marker2 <- shinyList$marker2  # less than 10

myhyper <- function(g1, g2) {  # Hypergeometric
  if(length(intersect(g1, g2)) == 0) return(1)
  1 - phyper(length(intersect(g1, g2)) - 1, length(g2), length(setdiff(bg, g2)), length(g1))
}  # Pr(count >= length(intersect(g1, g2)))

enrich1 <- sapply(profile1, function(y) sapply(marker2, function(x) myhyper(y, x)))
enrich2 <- sapply(profile2, function(y) sapply(marker1, function(x) myhyper(y, x)))

save(profile1, profile2, enrich, file = "Marker/marker.rdt")

pdf(file = "Figure/enrich1.pdf", width = 12)
pheatmap(t(enrich1))
dev.off()

pdf(file = "Figure/enrich2.pdf", width = 15)
pheatmap(t(enrich2), display_numbers = T)
dev.off()
