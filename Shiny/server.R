library(shiny)
library(ggplot2)
library(mail)
library(ape)
library(amap)

load("./brain_bc2014.rdt")
load("./gene_summary.rdt")
load("./glm_brain.rdt")
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
  
  output$glm_graph <- renderPlot({
    geneId.profile <- glm[[input$select3]][[input$type]]$symbol_by_pattern
    profile.Id <- names(geneId.profile)
    if (input$type == "app") {
      mylab <- c("2m:APP", "4m:APP", "5m:APP", "6m:APP")
    } else {
      mylab <- c("4m:WT", "5m:WT", "6m:WT", "4m:APP", "5m:APP", "6m:APP")
    }

    if (input$tissue2 == 2) {
      load("./glm_retina.rdt")
      geneId.profile <- glm[[input$select4]][[input$type]]$symbol_by_pattern
      profile.Id <- names(geneId.profile)
      if (input$type == "app") {
        mylab <- c("2m:APP", "5m:APP", "6m:APP")
      } else {
        mylab <- c("5m:WT", "6m:WT", "5m:APP", "6m:APP")
      }
    } 

    tile.dt <- NULL
    for (i in 1:length(profile.Id))
      tile.dt <- rbind(tile.dt, as.logical(unlist(strsplit(profile.Id[i], "-"))))
    tile.dt <- data.frame(value = c(tile.dt), 
      profile = factor(rep(profile.Id, ncol(tile.dt)), levels = profile.Id),
      group = factor(rep(mylab, each = nrow(tile.dt)), levels = mylab))

    ggplot(tile.dt, aes(x = group, y = profile, fill = value)) + geom_tile(colour = "white") +
      theme_bw() + xlab("") + ylab("") + coord_flip() +
      scale_y_discrete(labels = sapply(geneId.profile, length)) +
      scale_fill_manual(values = c("grey80", "firebrick1")) 
  })

  output$glm_table <- renderTable({
    table1 <- glm[[input$select3]][[input$type]]
    if (input$tissue2 == 2) {
      load("./glm_retina.rdt")
      table1 <- glm[[input$select4]][[input$type]]
    } 

    if (input$result == "symbol_by_pattern") {
      x <- table1[["symbol_by_pattern"]] 
      x <- lapply(x, function (x) paste(x, collapse = ", "))
      table2 <- do.call(rbind, x)
      colnames(table2) <- "Gene Symbol"
      table2 <- apply(table2, 2, rev)
    } else {
      table2 <- table1[[input$result]]
      table2$Symbols = gsub(";", "; ", table2$Symbols)
    } 
    table2
  }, include.rownames = TRUE)

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
        table1 <- five$gk[[input$result3]]
        table1$Symbols = gsub(";", "; ", table1$Symbols)
      }
      table1
  })

  output$downloadData <- downloadHandler(
    filename = "feb17talk.html", content = function (file) file.copy("./talk.html", file)
  )

})
