library(shiny)
library(ggvis)

load("edges.rdt")

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
    p(strong("GLM:"), "Model gene expression data with a generalized linear model (GLM).", strong(a(href="glm.html", "Click here for the protocol."))),
    hr(),
    fluidRow(
      column(3,
	radioButtons("tissue2", label = h3("Choose tissue"), choices = list("Brain" = 1, "Retina" = 2), selected = 1),
	radioButtons("type", label = h3("App or Age genes?"), 
	             choices = list("App" = "app", "Age" = "age"), selected = "app"),
	selectInput("select3", label = h3("Less or more brain genes?"), 
	            choices = list("Less: R2 > 0.5, |effect| > 0.2" = "less", "More: R2 > 0.4, |effect| > 0.1" = "more"), selected = "less"),
	selectInput("select4", label = h3("Less or more retina genes?"), 
	            choices = list("Less: R2 > 0.5, |effect| > 0.5" = "less", "More: R2 > 0.5, |effect| > 0.2" = "more"), selected = "less"),
	p()
      ),
      column(9,
        plotOutput("glm_graph"),
        textInput("pattern", label = h4("Input a pattern"), value = 10),
        helpText("First number in X-axis label indicates the group number."),
        h4("Gene symbol"),
        textOutput("glm_table1"),
        h4("KEGG enrichment"),
        helpText("p-value < 0.01"),
        tableOutput("glm_table2"),
        h4("GO_BP enrichment"),
        helpText("p-value < 0.001"),
        tableOutput("glm_table3"),
        h4("GO_MF enrichment"),
        helpText("p-value < 0.001"),
        tableOutput("glm_table4"),
        h4("GO_CC enrichment"),
        helpText("p-value < 0.001"),
        tableOutput("glm_table5")
      )
    )
  )),
  tabPanel("Regulator",
    fluidRow(
      column(2,
        radioButtons("genes_by_pattern", "Choose a GLM gene cohorts", choices = names(edgesList), selected = names(edgesList[4])),
        helpText("The binary codes 0001 correlate with FALSE-FALSE-FALSE-TRUE pattern in the GLM tab."),
        radioButtons("layout", "Choose graph layout", choices = list("Sphere" = 1, "Circle" = 2, "Colony" = 3), selected = 1)
      ),
      column(10,
        tags$div(align = "center", p(strong("Transcription"), "factors in gold dots and red labels. Target genes in green dots and blue labels. Size represents number of connections.", strong(a(href="regulator.html", "Click here for the protocol.")))),
        plotOutput("igraph")
      )
  )),
  tabPanel("Adult APP",
           fluidPage(
             p(strong("Stratify the adult APP samples:"), strong(a(href="five.html", "Click here for the protocol."))),
             hr(),
             fluidRow(
               column(3,
                 radioButtons("result3", label = h3("Choose a table to display"), 
                             choices = list("Genes" = "symbol", "KEGG: p-values < 0.01" = "KEGG",
                                            "BP: p-values < 0.001" = "BP", "MF: p-values < 0.001" = "MF", "CC: p-values < 0.001" = "CC"), selected = "symbol")
               ),
               column(9,
                 plotOutput("five_graph"),
                 tableOutput("five_table")
               )
             )
           )),
  tabPanel("Board",
           h3("Compare the retinal and cerebral samples across different GLM patterns in the perspective of GO and KEGG pathways."),
           h1("To extract the prominent signal in both age and app context, parse all GLM parameters together without separating them into app- and age-related groups."),
           h2("Correlate GLM signals with the single cell RNA-seq data from 3005 cerebral cells.")
  )
))
    
