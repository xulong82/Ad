library(xlsx)
library(HiveR)
library(igraph)

rm(list = ls())
setwd("~/Dropbox/GitHub/Ad")

load("./data/glm_brain.rdt")

# ***
symbol <- glm$less$app$symbol_by_pattern[["FALSE-FALSE-FALSE-TRUE"]]
symbol <- symbol[sapply(fit.coef[symbol], function(x) x["age6m:groupAPP", "Estimate"] > 0)]
# ***
symbol <- glm$less$app$symbol_by_pattern[["TRUE-FALSE-FALSE-FALSE"]]
symbol <- symbol[sapply(fit.coef[symbol], function(x) x["groupAPP", "Estimate"] > 0)]
# ***
symbol <- glm$more$app$symbol_by_pattern[["FALSE-FALSE-FALSE-TRUE"]]
symbol <- symbol[sapply(fit.coef[symbol], function(x) x["age6m:groupAPP", "Estimate"] > 0)]
# ***
symbol <- glm$more$app$symbol_by_pattern[["TRUE-FALSE-FALSE-FALSE"]]
symbol <- symbol[sapply(fit.coef[symbol], function(x) x["groupAPP", "Estimate"] > 0)]
# ***
symbol <- glm$less$age$symbol_by_pattern[["FALSE-TRUE-TRUE-FALSE-FALSE-FALSE"]]
symbol <- symbol[sapply(fit.coef[symbol], function(x) x["age5m", "Estimate"] > 0)]
# ***
symbol <- glm$more$age$symbol_by_pattern[["FALSE-TRUE-TRUE-FALSE-FALSE-FALSE"]]
symbol <- symbol[sapply(fit.coef[symbol], function(x) x["age5m", "Estimate"] > 0)]
# ***
symbol <- glm$less$age$symbol_by_pattern[["FALSE-FALSE-TRUE-FALSE-FALSE-FALSE"]]
symbol <- symbol[sapply(fit.coef[symbol], function(x) x["age6m", "Estimate"] > 0)]
# ***
symbol <- glm$more$age$symbol_by_pattern[["FALSE-FALSE-TRUE-FALSE-FALSE-FALSE"]]
symbol <- symbol[sapply(fit.coef[symbol], function(x) x["age6m", "Estimate"] > 0)]
# ***

write.xlsx(symbol, file = "Pathway/symbols.xlsx", sheetName = "new", append = T)

# symbol <- glm$more$app$symbol_by_pattern; symbol <- symbol[sapply(symbol, length) > 10]
# symbol <- glm$more$age$symbol_by_pattern; symbol <- symbol[sapply(symbol, length) > 10]
# symbol <- lapply(1:length(symbol), function (x) c(names(symbol)[x], symbol[[x]]))
# ***

system("rm Pathway/1.txt"); lapply(symbol, write, "Pathway/1.txt", append = T, ncolumns = 1e3)

# --- PARSE THE IREGULON OUTPUT. MORE THAN 10 GENES REQUIRED FOR IREGULON
load("./data/fit_brain.rdt")

file <- "Pathway/regulon_brain_less_app_0001.tsv"
file <- "Pathway/regulon_brain_less_app_1000.tsv"
file <- "Pathway/regulon_brain_more_app_0001.tsv"
file <- "Pathway/regulon_brain_more_app_1000.tsv"
file <- "Pathway/regulon_brain_less_age_011000.tsv"
file <- "Pathway/regulon_brain_more_age_011000.tsv"
file <- "Pathway/regulon_brain_less_age_001000.tsv"
file <- "Pathway/regulon_brain_more_age_001000.tsv"

ireg <- read.delim(file, comment.char = ";", stringsAsFactors = F)

pars <- c("(Intercept)", "groupAPP")  # 1000
pars <- c("(Intercept)", "age6m:groupAPP")  # 0001
pars <- c("(Intercept)", "age5m")  # 011000
pars <- c("(Intercept)", "age6m")  # 001000

univ <- names(fit.coef)[sapply(fit.coef, function (x) sum(x[pars, 1]) > 5 & x[pars[2], 1] > 0)]

factor <- lapply(1:nrow(ireg), function(x) unlist(strsplit(ireg$Transcription.factor[x], split = ",")))
target <- lapply(1:nrow(ireg), function(x) unlist(strsplit(ireg$Target.genes[x], split = ",")))
factor <- lapply(factor, function(x) x[x %in% univ])
idx <- sapply(factor, length) > 0; factor <- factor[idx]; target <- target[idx]
edges <- lapply(1:length(factor), function(x) expand.grid(factor[[x]], target[[x]], stringsAsFactors = F))
edges <- do.call(rbind, edges); edges <- edges[! duplicated(edges), ]

edgesList <- list()
edgesList[[gsub("Pathway/", "", file)]] <- edges
save(edgesList, file = "Shiny/edges.rdt")

# VISUALIZATION: IGRAPH
igraph.dt <- graph.data.frame(edges)

igraph.dt$layout <- layout.sphere
igraph.dt$layout <- layout.circle
igraph.dt$layout <- layout.fruchterman.reingold 

V(igraph.dt)$color = rep("chartreuse3", length(V(igraph.dt)$name))
V(igraph.dt)$color[V(igraph.dt)$name %in% unlist(factor)] <- "gold"
V(igraph.dt)$size = (degree(igraph.dt) - 1) / (max(degree(igraph.dt)) - 1) * 10 + 2
V(igraph.dt)$label.cex = (degree(igraph.dt) - 1) / (max(degree(igraph.dt)) - 1) * 2 + 0.5
V(igraph.dt)$label.color = rep("dodgerblue3", length(V(igraph.dt)$name))
V(igraph.dt)$label.color[V(igraph.dt)$color == "gold"] = "firebrick1"

plot.igraph(igraph.dt, vertex.frame.color = "white", edge.arrow.size = 0.3)

# --- VISUALIZATION: ARC DIAGRAM, HIVE PLOT
# nodes <- data.frame(lab = unique(c(unlist(factor), unlist(target))))
# nodes$lab <- as.character(nodes$lab)
# nodes$axis <- rep(2, nrow(nodes))
# nodes$axis[nodes$lab %in% unlist(factor)] <- 1
# nodes$axis <- as.integer(nodes$axis)
# nodes$radius <- 0
# nodes$radius[nodes$lab %in% unlist(factor)] <- 1:sum(nodes$lab %in% unlist(factor))
# nodes$radius[nodes$lab %in% unlist(target)] <- 1:sum(nodes$lab %in% unlist(target))
# nodes$id <- 1:nrow(nodes)
# nodes$color <- "grey30"
# nodes$color[nodes$lab %in% unlist(factor)] <- "firebrick1"
# nodes$size <- rep(1, nrow(nodes))
# mytable <- table(c(as.matrix(edges[, 1:2])))
# nodes$size <- as.numeric(mytable[nodes$lab]) / max(mytable) * 2
# 
# edges <- NULL
# for (i in 1:length(factor)) 
#   edges <- rbind(edges, expand.grid(factor[[i]], target[[i]], stringsAsFactors = F))
# edges$id1 <- nodes$id[match(edges$Var1, nodes$lab)]
# edges$id2 <- nodes$id[match(edges$Var2, nodes$lab)]
# edges$color <- "grey60"
# edges$weight <- as.numeric(table(edges$Var1)[edges$Var1]) / 10
# edges$weight <- rep(1, nrow(edges))
# 
# hpd = list()
# hpd$nodes = nodes
# hpd$edges = edges
# hpd$type = "2D"
# hpd$desc = "Hive Plot"
# hpd$axis.cols = "white"
# hpd$axLabs = c("TF","Target")
# class(hpd) = "HivePlotData"
# chkHPD(hpd, confirm = TRUE)
# 
# mylab = nodes$lab[nodes$axis == 1]
# plotHive(hpd, axLabs = hpd$axLabs, axLab.pos = 2, ch = 0.1, bkgnd = "white",
#          axLab.gpar = gpar(col = "chartreuse3", fontsize = 12))
# grid.text(mylab, x = 2, y = 0.25 + (1:length(mylab)), gp = gpar(fontsize = 10))
