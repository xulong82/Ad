library(GOstats)
library(KEGG.db)
library(org.Mm.eg.db)
library(Category)
library(dplyr)

# GO AND KEGG ENRICHMENT
myGK <- function (geneId) {
  entrezId <- mget(geneId, org.Mm.egSYMBOL2EG, ifnotfound = NA) %>% unlist %>% na.omit %>% as.character
  if(! length(entrezId) > 1) return(NA) # use nsFilter() for additional filtering if required
  
  uni_go <- get("org.Mm.egGO") %>% Lkeys
  uni_kg <- get("org.Mm.egPATH") %>% Lkeys
  
  GO <- lapply(c("BP", "MF", "CC"), function(category) {
    over <- new("GOHyperGParams", geneIds = entrezId, universeGeneIds = uni_go, annotation = "org.Mm.eg.db", 
                ontology = category, pvalueCutoff = 0.001, testDirection = "over") %>% hyperGTest 
    glist <- sapply(geneIdsByCategory(over), function(x) {y <- mget(x, envir=org.Mm.egSYMBOL); paste(y, collapse=";")})
    go <- summary(over); go$Symbols <- glist[as.character(go[, paste0("GO", category, "ID")])]; go
  }); names(GO) <- c("BP", "MF", "CC")
  
  over <- new("KEGGHyperGParams", geneIds=entrezId, universeGeneIds=uni_kg, annotation="org.Mm.eg.db", 
              categoryName="KEGG", pvalueCutoff = 0.01, testDirection="over") %>% hyperGTest
  glist <- sapply(geneIdsByCategory(over), function(x) {y <- mget(x, envir=org.Mm.egSYMBOL); paste(y, collapse=";")})
  kegg <- summary(over); kegg$Symbols <- glist[as.character(kegg$KEGGID)]
    
  gk <- list(); gk$GO <- GO; gk$KEGG <- kegg; return(gk)
}
