tabPanel("Modules",
         fluidPage(
           p(strong("Modules:"), "Identify co-expression modules (CEM) with App and Age phenotype.", strong(a(href="network.html", "Click here for the protocol"))),
           hr(),
           sidebarLayout(
             sidebarPanel(
               radioButtons("tissue3", label = h3("Choose tissue"), choices = list("Brain" = 1, "Retina" = 2), selected = 1),
               selectInput("select1", label = h3("Select a brain module"), 
                           choices = list("Module 1" = 1, "Module 2" = 2, "Module 3" = 3, "Module 4" = 4, 
                                          "Module 5" = 5, "Module 6" = 6, "Module 7" = 7, "Module 8" = 8), selected = 6),
               selectInput("select2", label = h3("Select a retina module"), 
                           choices = list("Module 1" = 1, "Module 2" = 2, "Module 3" = 3, "Module 4" = 4), selected = 2),
               radioButtons("result2", label = h3("Choose a table to display"), 
                            choices = list("Symbol and summary" = "symbol", "KEGG: p-values < 0.01" = "kegg"), selected = "kegg"),
               p()
             ),
             mainPanel(
               plotOutput("network_graph"),
               tableOutput("network_table")
             )
           )
         )),

output$network_graph <- renderPlot({
  if (input$tissue3 == 2) load("./network_retina.rdt")
  gdt <- network$gdt
  ggplot(gdt, aes(x = group, y = value, fill = geno)) + 
    geom_bar(stat = "identity") + facet_grid(. ~ module) +
    scale_fill_manual(values = c("red", "blue")) +
    theme_bw() + xlab("") + ylab("") + 
    theme(axis.text.x = element_text(angle = 90))
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
