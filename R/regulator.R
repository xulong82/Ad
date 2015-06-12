library(xlsx)
library(HiveR)
library(igraph)

rm(list = ls())
setwd("~/Dropbox/GitHub/Ad")

load("./data/glmList.rdt")

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

write.xlsx(symbol, file = "symbols.xlsx", sheetName = "new", append = T)

# symbol <- glm$more$age$symbol_by_pattern; symbol <- symbol[sapply(symbol, length) > 10]
# symbol <- lapply(1:length(symbol), function (x) c(names(symbol)[x], symbol[[x]]))
# system("rm Regulator/1.txt"); lapply(symbol, write, "Regulator/1.txt", append = T, ncolumns = 1e3)

myEdges <- function(file) {  # PARSE THE IREGULON OUTPUT
  ireg <- read.delim(paste("Regulator", file, sep = "/"), comment.char = ";", stringsAsFactors = F)
  
  factor <- lapply(1:nrow(ireg), function(x) unlist(strsplit(ireg$Transcription.factor[x], split = ",")))
  target <- lapply(1:nrow(ireg), function(x) unlist(strsplit(ireg$Target.genes[x], split = ",")))
  factor <- lapply(factor, function(x) x[x %in% names(fit.coef)])  # only mild expression, post-hoc filter required
  idx <- sapply(factor, length) > 0; factor <- factor[idx]; target <- target[idx]
  edges <- lapply(1:length(factor), function(x) expand.grid(factor[[x]], target[[x]], stringsAsFactors = F))
  edges <- do.call(rbind, edges); edges <- edges[! duplicated(edges), ]
  
  return(edges)
}

file <- list.files(path = "./Regulator", pattern = "*.tsv")
edgesList <- lapply(file, myEdges)
names(edgesList) <- gsub("regulon_(.*).tsv", "\\1", file)
save(edgesList, file = "Shiny/edges.rdt")

# VISUALIZATION: IGRAPH
myIgraph <- function (edges) {
  igraph.dt <- graph.data.frame(edges)
  igraph.dt$layout <- layout.sphere
  
  V(igraph.dt)$color = rep("chartreuse3", length(V(igraph.dt)$name))
  V(igraph.dt)$color[V(igraph.dt)$name %in% edges$Var1] <- "gold"
  V(igraph.dt)$size = (degree(igraph.dt) - 1) / (max(degree(igraph.dt)) - 1) * 10 + 2
  V(igraph.dt)$label.cex = (degree(igraph.dt) - 1) / (max(degree(igraph.dt)) - 1) * 2 + 0.8
  V(igraph.dt)$label.color = rep("dodgerblue3", length(V(igraph.dt)$name))
  V(igraph.dt)$label.color[V(igraph.dt)$color == "gold"] = "firebrick1"
  
  plot.igraph(igraph.dt, vertex.frame.color = "white", edge.arrow.size = 0.3)
}

myIgraph(edgesList[[1]])

sapply(edgesList, function(x) unique(x$Var1[grep("Stat", x$Var1)]))  # Stat* appears frequently

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
