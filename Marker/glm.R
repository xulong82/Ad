library(rstan)  # R/3.1.1
library(dplyr)
library(parallel)

rm(list = ls())
setwd("~/Dropbox/GitHub/Load")

glm <- stan_model(file = "./GLM/glm.stan", model_name='GLM')
myGLM <- function(x) {  # GLM model
  fit1 <- lapply(1:nrow(x), function(i) { glm.dt$y <- x[i, ] 
    sampling(glm, data = glm.dt, warmup = 3e2, iter = 6e2, chains = 3)
  }); names(fit1) <- rownames(x); return(fit1)
}

load("./GLM//ge.rdt")  # RNA-seq data
grp1 <- c("B6", "ApoE4", "Apoe", "Bin1", "Cd2ap", "Clu") 
grp2 <- gsub("-.*", "", colnames(ge))
x <- sapply(grp1, function(x) as.numeric(grp2 == x))
glm.dt <- list(N = nrow(x), K = ncol(x), x = x)

ge <- ge[apply(ge, 1, function(x) sum(x > 50) > 1), ]  # Choose genes
ge <- log2(ge + 1)

u.core <- round(nrow(ge) / 20)  # PC
idx1 <- (1:20 - 1) * u.core + 1; idx2 <- c(1:19 * u.core, nrow(ge))
geList <- lapply(1:20, function(x) ge[idx1[x]:idx2[x], ])

y <- mclapply(geList, myGLM, mc.cores = 20)
fit <- do.call(c, y)

save(fit, file = "/data/xwang/Load/GLM/fit_sampling.rdt")

