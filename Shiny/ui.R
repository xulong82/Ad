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
	radioButtons("number", label = h3("Less or more genes?"), 
	             choices = list("Less: R2 > 0.5, |effect| > 0.2" = "less", "More: R2 > 0.4, |effect| > 0.1" = "more"), selected = "less"),
	radioButtons("type", label = h3("App or Age genes?"), 
	             choices = list("App" = "app", "Age" = "age"), selected = "app"),
	radioButtons("result", label = h3("Choose a table to display"), 
	             choices = list("Symbol by group" = "symbol_by_pattern", "Symbol and summary" = "symbol", "KEGG: p-values < 0.01" = "go_kegg", 
	                            "BP: p-values < 0.001" = "go_bp", "MF: p-values < 0.001" = "go_mf", "CC: p-values < 0.001" = "go_cc"), selected = "go_kegg"),
	p()
      ),
      mainPanel(
        tableOutput("glm")
      )
    )
  )),
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
  tabPanel("Adult APP",
           fluidPage(
             p(strong("Stratify the adult APP samples:"), strong(a(href="five.html", "Click here for the protocol"))),
             hr(),
             sidebarLayout(
               sidebarPanel(
                 radioButtons("result3", label = h3("Choose a table to display"), 
                              choices = list("Symbol and summary" = "symbol", "KEGG: p-values < 0.01" = "kegg"), selected = "kegg"),
                 p()
               ),
               mainPanel(
                 tableOutput("five_table")
               )
             )
           )),
  tabPanel("Board",
           h3("To do list:"),
           tags$ol(
             tags$li("There are way more significant genes in retinal samples. So much that it is hard to extract meaningful things. I possibly need to put harder stringencies."), 
             tags$li("Extensive pathway analysis need to be applied on misc gene lists. iRegulon will be the tool for this purpose."), 
             tags$li("Some comparison between the retinal and cereberal samples in details.")),
           hr(),
           h4("Leave XuLong a message?"),
           tags$style(type='text/css', "#message {height: 300px; width:600px}"),
           tags$textarea(id = "message", placeholder = 'Type here', rows = 20, ""),
           submitButton("Submit")
           )
))
    
