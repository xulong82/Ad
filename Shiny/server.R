library(shiny)
library(ggplot2)
library(mail)

load("./brain_bc2014.rdt")
load("./gene_summary.rdt")
load("./glm_brain.rdt")
load("./network_brain.rdt")
load("./five.rdt")

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
  
  output$glm <- renderTable({
    if (input$tissue2 == 2) load("./glm_retina.rdt")

    table1 <- glm[[input$number]][[input$type]]
    if (input$result == "symbol") {
      x <- table1$symbol
      table2 <- table[table$query %in% x, ]
      rownames(table2) <- 1:nrow(table2)
    } else if (input$result == "symbol_by_pattern") {
      x <- table1[["symbol_by_pattern"]] 
      x <- lapply(x, function (x) paste(x, collapse = ", "))
      table2 <- do.call(rbind, x)
      colnames(table2) <- "Gene Symbol"
    } else {
      table2 <- table1[[input$result]]
      table2$Symbols = gsub(";", "; ", table2$Symbols)
    } 
    table2
  }, include.rownames = TRUE)

  output$network_graph <- renderPlot({
    if (input$tissue3 == 2) load("./network_retina.rdt")
    gdt <- network$gdt
    ggplot(gdt, aes(x = group, y = value, fill = geno)) + 
      geom_bar(stat = "identity") + facet_grid(. ~ module) +
      scale_fill_manual(values = c("red", "blue")) +
      theme_bw() + xlab("") + ylab("") + 
      theme(axis.text.x = element_text(angle = 90))
  })

  output$five_graph <- renderPlot({
    mycol <- five$gdt$mycol
    myhc <- five$gdt
    plot(myhc, edge.width = 2, font = 2, cex = 0.7, label.offset = 1e-4, tip.color = mycol, direction = "downward")
  })

  output$five_table <- renderTable({
      if (input$result3 == "symbol") {
        x <- five$symbol
        table1 <- table[table$query %in% x, ]
        rownames(table1) <- 1:nrow(table1)
      } else {
        table1 <- five$gk$KEGG
        table1$Symbols = gsub(";", "; ", table1$Symbols)
      }
      table1
  })

  output$network_table <- renderTable({
      module <- as.numeric(input$select1)
      if (input$tissue3 == 2) {
        load("./network_retina.rdt")
        module <- as.numeric(input$select2)
      }
      if (input$result2 == "symbol") {
        x <- network$symbol[[module]]
        table1 <- table[table$query %in% x, ]
        rownames(table1) <- 1:nrow(table1)
      } else {
        table1 <- network$gk[[module]]$KEGG
        table1$Symbols = gsub(";", "; ", table1$Symbols)
      }
      table1
  })

  output$downloadData <- downloadHandler(
    filename = "feb17talk.html", content = function (file) file.copy("./talk.html", file)
  )

  output$message <- renderText({
    message <- input$message
    sendmail("xulong.wang@jax.org", "AD Notice", as.character(message))
  })

})

