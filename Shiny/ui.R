library(shiny)
library(ggvis)

shinyUI(
  fluidPage(
    titlePanel("Alzheimer's Disease Project"),
    p(strong("Figure 1:"), "Boxplot of gene expression profiles of mouse brain and retina. Data were collected at 4 time points and 2 genotypes."),
    hr(),
    sidebarLayout(
      sidebarPanel(
	      radioButtons("tissue", label = h3("Choose tissue"), choices = list("Brain" = 1, "Retina" = 2), selected = 1),
        textInput("gene", label = h3("Input gene:"), value = "App"),
	      hr(),
	      h5("Useful documents"),
	      p("Joint lab meeting on Feb 17, 2015"),
	      downloadButton('downloadData', 'Download')
      ),
      mainPanel(
        plotOutput("boxplot"),
	      tableOutput("summary")
      )
    )
))
    
