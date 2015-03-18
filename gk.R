rm(list = ls())
setwd("~/Dropbox/GitHub/Ad")
source("function.R")

# --- load either of the below
load("data/glm_brain.rdt")
load("data/glm_retina.rdt")

# -------------------------------------------------------------------------
theId <- theGK <- list()
theId[[1]] <- glm$more$app$symbol_by_pattern
theId[[2]] <- glm$more$age$symbol_by_pattern
theId[[3]] <- glm$less$app$symbol_by_pattern
theId[[4]] <- glm$less$age$symbol_by_pattern

for(i in 1:length(theId)) {
  mygk1 <- list()
  gene1 <- theId[[i]]
  for (j in 1:length(gene1)) {
    gene0 <- gene1[[j]]
    cat(i, j, length(gene1), "\n")
    if (length(gene0) > 1)
      mygk1[[j]] <- myGK(gene0)
    else
      mygk1[[j]] <- "Not tested, only one gene!"
  }
  theGK[[i]] <- mygk1
}

glm$more$app$gk_by_pattern <- theGK[[1]]
glm$more$age$gk_by_pattern <- theGK[[2]]
glm$less$app$gk_by_pattern <- theGK[[3]]
glm$less$age$gk_by_pattern <- theGK[[4]]
