# -------------------------------------------------------------------------
myfitAic <- list()
dt.aic <- dt[which(qval < 0.05), ]
for (i in 1:nrow(dt.aic)) {
  if (i %% 1e2 == 0) cat(i, "\n")
  y <- dt.aic[i, ]
  fitLr <- lm(y ~ age + group + batch + age * group)
  fitAic <- stepAIC(fitLr, direction = "backward", trace = 0)
  myfitAic[[i]] <- summary(fitAic)
}

r2 <- sapply(myfitAic, function (x) x$r.squared)
fval <- sapply(myfitAic, function (x) x$fstatistic)
pval <- apply(fval, 2, function (x) pf(x[1], x[2], x[3], lower.tail = F))
qval <- p.adjust(pval, method = "fdr")

# --- grouping ---
cond <- c("age4m", "age5m", "age6m", "groupAPP", "batchmouse", "age4m:groupAPP", "age5m:groupAPP", "age6m:groupAPP")
grouping <- matrix(FALSE, nrow = length(geneId), ncol = length(cond), dimnames = list(geneId, cond))
for (i in 1:length(geneId)) {
  if (i %% 1e3 == 0) cat(i, "\n")
  fit0 <- myfitAic[[i]]
  fit0.r2 <- fit0$r.squared
  fit0.coef <- fit0$coefficients
  cond0 <- rownames(fit0.coef)[fit0.coef[, "Pr(>|t|)"] < 0.01]
  # grouping[i, ] <- cond %in% cond0
  grouping[i, ] <- cond %in% cond0 & fit0.r2 > 0.5
}

grouping.geneId <- list()
for (idx in colnames(grouping)) grouping.geneId[[idx]] <- c(idx, rownames(grouping)[grouping[, idx]])

number <- sapply(grouping.geneId, length)
col.manual <- c("grey70", "firebrick1", "chartreuse3", "dodgerblue3", "gold1", "darkorchid2")
col.manual <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf")
pdf(file = "./Graphs/glm_number.pdf", width = 6, height = 12)
op <- par(mar = c(15, 4, 4, 2))
bar <- barplot(number, las = 3, col = col.manual, border = NA, ylim = c(0, 300), axes = F)
abline(0, 0, lwd = 5, col = "black")
text(x = bar, y = number + 5, labels = number, font = 2, col = "grey20")
dev.off()

lapply(grouping.geneId, write, "./DE/glm.txt", append = TRUE, ncolumns = 1000)

geneId.de = geneId[as.logical(apply(grouping[, -5], 1, max))]
write.table(geneId.de, file = "./DE/geneId.txt", sep = "\t", row.names = F, col.names = F, quote = F)
