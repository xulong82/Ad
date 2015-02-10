library(VennDiagram)
library(gplots)
library(ggplot2)

# --- GENE INTERSECT
venn.diagram(list(APP = geneId.app, AGE = geneId.age),
             fill = c("dodgerblue3", "firebrick1"), col = "transparent", alpha = .5, 
             cat.col = c("dodgerblue3", "firebrick1"), cat.cex = 2,
             cex = 3, lwd = 1, file = "./Graphs/glm_venn.tiff")

# --- TILE:APP
tile.app <- NULL
for (i in 1:length(profile.app.Id))
  tile.app <- rbind(tile.app, as.logical(unlist(strsplit(profile.app.Id[i], "-"))))
  tile.app <- data.frame(value = c(tile.app), 
                        profile = factor(rep(profile.app.Id, ncol(tile.app)), levels = profile.app.Id),
                        group = factor(rep(colnames(profile.app), each = nrow(tile.app)), levels = colnames(profile.app)))
  
pdf(file = "./Graphs/tile_app.pdf", width = 2, height = 4)
  ggplot(tile.app, aes(x = group, y = profile, fill = value)) + geom_tile(colour = "white") +
    theme_bw() + xlab("") + ylab("") +
    theme(panel.border = element_rect(size = .5, color = "grey30")) +
  theme(axis.text.y = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold", angle = 90)) +
  scale_x_discrete(labels = c("2m:APP", "4m:APP", "5m:APP", "6m:APP")) +
  scale_y_discrete(labels = profile.app.table) +
  scale_fill_manual(values = c("grey80", "firebrick1")) +
  theme(legend.position = "none") 
dev.off()

# --- TILE:AGE
tile.age <- NULL
for (i in 1:length(profile.age.Id))
  tile.age <- rbind(tile.age, as.logical(unlist(strsplit(profile.age.Id[i], "-"))))
  tile.age <- data.frame(value = c(tile.age), 
                        profile = factor(rep(profile.age.Id, ncol(tile.age)), levels = profile.age.Id),
                        group = factor(rep(colnames(profile.age), each = nrow(tile.age)), levels = colnames(profile.age)))
  
pdf(file = "./Graphs/tile_age.pdf", width = 3, height = 8)
  ggplot(tile.age, aes(x = group, y = profile, fill = value)) + geom_tile(colour = "white") +
    theme_bw() + xlab("") + ylab("") +
    theme(panel.border = element_rect(size = .5, color = "grey30")) +
  theme(axis.text.y = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold", angle = 90)) +
  scale_x_discrete(labels = c("4m:WT", "5m:WT", "6m:WT", "4m:APP", "5m:APP", "6m:APP")) +
  scale_y_discrete(labels = profile.age.table) +
  scale_fill_manual(values = c("grey80", "firebrick1")) +
  theme(legend.position = "none") 
dev.off()

# --- PHYLO
mycol <- rep("grey50", ncol(dt.hc))
mycol[grep("APP", colnames(dt))] <- "firebrick1"

pdf("./Graphs/glm_phylo.pdf", fonts = "Helvetica", width = 10)
plot(as.phylo(hc1), edge.width = 2, font = 2, cex = 0.7, label.offset = 1e-4, tip.color = mycol, direction = "downward")
dev.off()

pdf("./Graphs/glm_phylo_gene.pdf", fonts = "Helvetica", width = 3, height = 9)
col.manual <- c("grey70", "firebrick1", "chartreuse3", "dodgerblue3", "gold1", "darkorchid2")
clusts = cutree(hc2, 6)
plot(as.phylo(hc2), tip.color = col.manual[clusts], cex = .2, 
     edge.width = 2, adj = 0.5, font = 2, label.offset = 1e-2, direction = "rightward")
dev.off()

tile.dt1 <- dt.hc[hc2$order, hc1$order]
tile.dt2 <- t(apply(tile.dt1, 1, scale))
tile.dt3 <- data.frame(value = c(tile.dt2), 
                       gene = factor(rep(rownames(tile.dt1), time = ncol(tile.dt1)), levels = rownames(tile.dt1)),
                       sample = factor(rep(colnames(tile.dt1), each = nrow(tile.dt1)), levels = colnames(tile.dt1)))

pdf("./Graphs/heatmap.pdf", width = 15, height = 10)
ggplot(tile.dt3, aes(x = sample, y = gene, fill = value, alpha = abs(value))) + 
  geom_tile() + guides(alpha = F) +
  scale_fill_gradient(low="blue", high="red") +
  theme_bw() + xlab("") + ylab("") + 
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 3), axis.ticks = element_blank()) +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.margin = unit(c(0, 0, -6, 0), "cm"), 
        legend.text = element_text(size = 10, face = "bold"), legend.title = element_blank())
dev.off()
