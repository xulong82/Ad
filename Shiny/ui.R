library(shiny)
library(ggvis)

shinyUI(
  navbarPage("Alzheimer's Disease Project",
  tabPanel("Expression",
  fluidPage(
    p(strong("Figure:"), "Boxplot of gene expression profiles of mouse brain and retina. Data were collected at 4 time points and 2 genotypes."),
    hr(),
    sidebarLayout(
      sidebarPanel(
	radioButtons("tissue", label = h3("Choose tissue"), choices = list("Brain" = 1, "Retina" = 2), selected = 1),
        textInput("gene", label = h3("Input gene:"), value = "App"),
	hr(),
	h5("Useful documents"),
	p("Joint lab meeting on Feb 17, 2015"),
	downloadButton('downloadData', 'Download'),
	hr(),
	p(a("GitHub Repository", href = "https://github.com/xulong82/Ad"))
      ),
      mainPanel(
        plotOutput("boxplot"),
	tableOutput("summary")
      )
    )
  )),
  tabPanel("GLM",
  fluidPage(
    p(strong("GLM:"), "Model gene expression data with a generalized linear model (GLM).", strong(a(href="glm.html", "Click here for the protocol"))),
    hr(),
    sidebarLayout(
      sidebarPanel(
	radioButtons("tissue2", label = h3("Choose tissue"), choices = list("Brain" = 1, "Retina" = 2), selected = 1),
	radioButtons("type", label = h3("App or Age genes?"), 
	             choices = list("App" = "app", "Age" = "age"), selected = "app"),
	selectInput("result", label = h3("Choose a table to display"), 
	             choices = list("Symbol by group" = "symbol_by_pattern", "KEGG: p-values < 0.01" = "go_kegg", 
	                            "BP: p-values < 0.001" = "go_bp", "MF: p-values < 0.001" = "go_mf", "CC: p-values < 0.001" = "go_cc"), selected = "symbol_by_pattern"),
	selectInput("select3", label = h3("Less or more brain genes?"), 
	            choices = list("Less: R2 > 0.5, |effect| > 0.2" = "less", "More: R2 > 0.4, |effect| > 0.1" = "more"), selected = "less"),
	selectInput("select4", label = h3("Less or more retina genes?"), 
	            choices = list("Less: R2 > 0.5, |effect| > 0.5" = "less", "More: R2 > 0.5, |effect| > 0.2" = "more"), selected = "less"),
	p()
      ),
      mainPanel(
        plotOutput("glm_graph"),
        tableOutput("glm_table")
      )
    )
  )),
  tabPanel("Adult APP",
           fluidPage(
             p(strong("Stratify the adult APP samples:"), strong(a(href="five.html", "Click here for the protocol"))),
             hr(),
             sidebarLayout(
               sidebarPanel(
                 selectInput("result3", label = h3("Choose a table to display"), 
                             choices = list("Genes" = "symbol", "KEGG: p-values < 0.01" = "KEGG",
                                            "BP: p-values < 0.001" = "BP", "MF: p-values < 0.001" = "MF", "CC: p-values < 0.001" = "CC"), selected = "symbol"),
                 p()
               ),
               mainPanel(
                 plotOutput("five_graph"),
                 tableOutput("five_table")
               )
             )
           )),
  tabPanel("Board",
           h3("To do list:"),
           tags$ol(
             tags$li("Extensive pathway analysis need to be applied on misc gene lists. iRegulon will be the tool for this purpose."), 
             tags$li("Some comparison between the retinal and cerebral samples in details."))
           )
))
    
