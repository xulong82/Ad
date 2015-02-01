rm(list = ls())

setwd("~/Dropbox/AD")
load(file = "./R/complete_tpm.rdt")

# --- C1qa
# brain.tpm <- brain.tpm[, -grep("APP5m1559.2014", colnames(brain.tpm))]

C1qa.brain <- brain.tpm["C1qa", ]
C1qa.retina <- retina.tpm["C1qa", ]

dt <- C1qa.brain
dt <- C1qa.retina

dt6m2013 <- dt[grep("2013", names(dt))]
dt6m2013 <- dt[grep("6m",names(dt6m2013))]
mygroup <- factor(gsub("^.*(WT|APP).*", "\\1", names(dt6m2013)), levels = c("WT", "APP"))
xx <- pairwise.t.test(c(as.matrix(dt6m2013)), mygroup)
xx <- lm(log2(c(as.matrix(dt6m2013)) + 1) ~ mygroup)

age <- factor(gsub("^.*(2m|4m|5m|6m).*", "\\1", names(dt)), levels = c("2m", "4m", "5m", "6m"))
group <- factor(gsub("^.*(WT|APP).*", "\\1", names(dt)), levels = c("WT", "APP"))
batch <- gsub("^.*(2013|2014|mouse).*", "\\1", names(dt))
batch <- gsub("mouse", "2014new", batch)
batch <- factor(batch, levels = c("2013", "2014", "2014new"))
uid <- paste(age, group, sep = "_")

graph.dt <- data.frame(value = c(as.matrix(dt)), age, group, batch)

pdf("./Graphs/C1qa_all.pdf", width = 10, height = 5)
ggplot(graph.dt, aes(x = age, y = value)) + 
  geom_boxplot(aes(shape = group, fill = group)) +
# facet_grid(. ~ batch) +
  theme_bw() + xlab("") + ylab("TPM") +
  theme(panel.border = element_rect(size = 1, color = "grey30")) +
  scale_fill_manual(values = c("chartreuse3", "firebrick1")) +
  theme(axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold", vjust = 1)) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank()) 
dev.off()

y <- log2(c(as.matrix(dt)) + 1)
fit0 <- lm(y ~ age)
fit0 <- lm(y ~ age + batch)
fit0 <- lm(y ~ group)
fit0 <- lm(y ~ group + batch)
fit0 <- lm(y ~ age + group)
fit0 <- lm(y ~ age + group + batch)

summary(fit0)

c1qa2m <- mean(as.matrix(c1qa[grep("2m_APP", uid)])) - mean(as.matrix(c1qa[grep("2m_WT", uid)]))
c1qa4m <- mean(as.matrix(c1qa[grep("4m_APP", uid)])) - mean(as.matrix(c1qa[grep("4m_WT", uid)]))
c1qa5m <- mean(as.matrix(c1qa[grep("5m_APP", uid)])) - mean(as.matrix(c1qa[grep("5m_WT", uid)]))
c1qa6m <- mean(as.matrix(c1qa[grep("6m_APP", uid)])) - mean(as.matrix(c1qa[grep("6m_WT", uid)]))
