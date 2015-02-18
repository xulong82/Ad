library(shiny)
library(ggplot2)

load("./brain_bc2014.rdt")
load("./gene_summary.rdt")

shinyServer(function(input, output) {
  
  output$boxplot <- renderPlot({
    tissue <- input$tissue
    geneId <- input$gene
    
    if (tissue == 2) load("./retina_bc2014.rdt")
    
    if (geneId %in% rownames(dt.bc))
      expr <- t(dt.bc[geneId, ])
    else
      expr <- rep(0, length(uid))
    
    graph.dt <- data.frame(value = expr, 
      cond = factor(uid, levels = conditions),
      geno = factor(group, levels = c("WT", "APP")))
    rownames(graph.dt) <- NULL
    colnames(graph.dt) <- c("value", "cond", "geno")
    
    ggplot(graph.dt, aes(x = cond, y = value, fill = geno)) + geom_boxplot() +
    xlab("") + ylab("") + ggtitle(input$gene) +
    scale_fill_manual(values = c("white", "firebrick1")) +
    theme(plot.title = element_text(size=20, face="bold", vjust=2))
  })

  output$summary <- renderTable({
    geneId <- input$gene
    table[table$query == geneId, ]
  }, include.rownames = FALSE)
  
  output$downloadData <- downloadHandler(
    filename = "feb17talk.html", content = function (file) file.copy("./talk.html", file)
  )
})
