library(shiny)
library(shinythemes)
library(shinycssloaders)
library(shinyMatrix)

navbarPage(theme = shinytheme("flatly"),
  "ArrayAnalysis",
  id = "Mainset",
  navbarMenu("Workflow",
             tabPanel("Data Input",
                sidebarLayout(
                  sidebarPanel(
                    width = 2,
                    fileInput(inputId = "DESCFILE", "Please select description file"),
                    checkboxInput(inputId = "HEADER",label =  "Header",value = TRUE),
                    radioButtons(inputId =  "SEP",label = "Separator",choices = c(Comma = ",", Semicolon = ";", Tab = "\t"),selected = "\t"),
                    radioButtons(inputId = "QUOTE",label = "Quote",choices = c(None = "", "Double Quote" = '"', "Single Quote" = "'")),
                    tags$hr(),
                    fileInput(inputId = "DATAFILE", "Please select .zip data file", accept = ".zip"),
                    tags$br(),
                    tags$br(),
                    actionButton(inputId = "next2",label = 'Next',icon = icon("arrow-right"))
                  ),
                  mainPanel(
                    tableOutput("desc"),
                    tableOutput("rawdata"),
                    textOutput("CDFnotify")
                  )
                )  
             ),
             tabPanel("Quality Control",
              tabsetPanel(id = "QCsubset",
                tabPanel("Sample Quality",
                  sidebarLayout(
                    sidebarPanel(
                      width = 2,
                      tags$b("Sample Quality"),
                      uiOutput(outputId = "sqInputs"),
                      actionButton(inputId = "prev1",label = 'Prev',icon = icon("arrow-left")),
                      actionButton(inputId = "next3",label = 'Next',icon = icon("arrow-right"))
                    ),
                    mainPanel(
                      plotOutput(outputId = "sampleprep",height = "100%", width = "100%"),
                      plotOutput(outputId = "ratioplot",height = "100%", width = "100%"),
                      plotOutput(outputId = "ratioplot2",height = "100%", width = "100%"),
                      plotOutput(outputId = "rnadegplot",height = "100%", width = "100%")
                    )
                  )
                ),
                tabPanel("Hybridization and overall signal quality",
                  sidebarLayout(
                    sidebarPanel(
                      width = 2,
                      tags$b("Hybridization and overall signal quality"),
                      uiOutput("hybridInput"),
                      actionButton(inputId = "prev2",label = 'Prev',icon = icon("arrow-left")),
                      actionButton(inputId = "next4",label = 'Next',icon = icon("arrow-right"))
                    ),
                    mainPanel(
                      plotOutput("hybrid",height = "100%", width = "100%"),
                      plotOutput("percentpresent",height = "100%", width = "100%"),
                      plotOutput("posnegdistro",height = "100%", width = "100%"),
                      plotOutput("backgroundintensity",height = "100%", width = "100%")
                    )
                  )
                ),
                tabPanel("Signal comparability and bias diagnostic",
                  sidebarLayout(
                    sidebarPanel(
                      width = 2,
                      tags$b("Signal Distribution"),
                      uiOutput("signalInput"),
                      #tags$b("Intensity-dependent bias"),
                      #checkboxInput(inputId = "MAraw",label = "MA plot",value = FALSE),
                      tags$b("Spatial bias"),
                      checkboxInput(inputId = "layoutplot",label = "Array reference layout",value = T),
                      checkboxInput(inputId = "posnegCOI",label = "Pos/Neg Center of intensity",value = T),
                      #checkboxInput(inputId = "spatialimages",label = "2D Images",value = FALSE),
                      #checkboxInput(inputId = "PLMimage",label = "All PLM based images",value = FALSE),
                      tags$b("Probe-set homogeneity"),
                      checkboxInput(inputId = "Nuse",label = "NUSE plot",value = T),
                      checkboxInput(inputId = "RLE",label = "RLE plot",value = T),
                      actionButton(inputId = "prev3",label = 'Prev',icon = icon("arrow-left")),
                      actionButton(inputId = "next5",label = 'Next',icon = icon("arrow-right"))
                    ),
                    mainPanel(
                      plotOutput("scalefactors",height = "100%", width = "100%"),
                      plotOutput("boxplotRaw",height = "100%", width = "100%"),
                      plotOutput("densityRaw",height = "100%", width = "100%"),
                      plotOutput("controlplotaffx",height = "100%", width = "100%"),
                      plotOutput("controlplotbox",height = "100%", width = "100%"),
                      #plotOutput("rawMAplot",height = "100%", width = "100%"),
                      plotOutput("layoutplot",height = "100%", width = "100%"),
                      plotOutput("posnegCOI",height = "100%", width = "100%"),
                      #plotOutput("spatialImages",height = "100%", width = "100%"),
                      #plotOutput("PLMimages",height = "100%", width = "100%"),
                      plotOutput("rawNuse",height = "100%", width = "100%"),
                      plotOutput("rawRle",height = "100%", width = "100%")
                    )
                  )
                ),
                tabPanel("Array correlation",
                  sidebarLayout(
                    sidebarPanel(
                      width = 2,
                      tags$b("Array correlation"),
                      checkboxInput(inputId = "correlRaw",label = "Correlation plot",value = T),
                      checkboxInput(inputId = "PCAraw",label = "PCA analysis",value = T),
                      checkboxInput(inputId = "clusterRaw",label = "Hierarchical clustering",value = T),
                      selectInput(inputId = "clusteroption1",label = "Distance calculation method",choices = c("Pearson","Spearman","Euclidian")),
                      selectInput(inputId = "clusteroption2",label = "Clustering method",choices = c("ward","single","complete","average","mcquitty","median","centroid")),
                      actionButton(inputId = "prev4",label = 'Prev',icon = icon("arrow-left")),
                      actionButton(inputId = "next6",label = 'Next',icon = icon("arrow-right"))
                    ),
                    mainPanel(
                      plotOutput("correlRaw",height = "100%", width = "100%"),
                      plotOutput("PCAraw",height = "100%", width = "100%"),
                      plotOutput("clusterRaw",height = "100%", width = "100%")
                    )
                  )
                )
              )
             ),
             tabPanel("Preprocessing and evaluation",
              tabsetPanel(id = "Presubset",
                tabPanel("Normalization method and annotation",
                  sidebarLayout(
                    sidebarPanel(
                      width = 2,
                      selectInput(inputId = "normMeth", label = "Normalization method", choices = c("RMA","GCRMA","PLIER","none")),
                      selectInput(inputId = "perGroup", label = "Normalization per experimental group or using all arrays", choices = c("dataset","group")),
                      selectInput(inputId = "annotations",label = "Annotations",choices = c("No annotations","Custom annotations","Upload annotation file")),
                      conditionalPanel(
                        condition = "input.annotations=='Custom annotations'",
                        checkboxInput(inputId = "customCDF", label = "Custom annotation file (CDF)", value = T),
                        selectInput(inputId = "species",label = "Species",choices = c("Anopheles gambiae","Arabidopsis thaliana","Bos taurus","Caenorhabditis elegans","Canis familiaris", "Danio rerio","Drosophila melanogaster","Gallus gallus","Homo sapiens","Macaca mulatta","Mus musculus", "Oryza sativa","Rattus norvegicus","Saccharomyces cerevisiae","Schizosaccharomyces pombe","Sus scrofa"), selected = "Homo sapiens"),
                        selectInput(inputId = "CDFtype",label = "Annotation format",choices = c("ENTREZG","REFSEQ","ENSG","ENSE","ENST","VEGAG","VEGAE","VEGAT","TAIRG","TAIRT","UG","MIRBASEF","MIRBASEG")),
                      ),
                      conditionalPanel(
                        condition = "input.annotations=='Upload annotation file'",
                        fileInput(inputId = "annot_file",label = "Upload annotation file",multiple = FALSE)
                      ),
                      actionButton(inputId = "preprocessing",label = "Run preprocessing", icon = icon("refresh")),
                      tags$br(),
                      tags$br(),
                      actionButton(inputId = "prev5", label = 'Prev',icon = icon("arrow-left")),
                      actionButton(inputId = "next7", label = 'Next',icon = icon("arrow-right"))
                    ),
                    mainPanel(
                      tableOutput("normdatatable")
                    )
                  )
                ),
                tabPanel("Signal comparability and bias of normalized intensities",
                  sidebarLayout(
                    sidebarPanel(
                      width = 2,
                      tags$b("Signal comparability and bias of normalized intensities"),
                      checkboxInput(inputId = "boxplotNorm", label = "Boxplot Norm", value = FALSE),
                      checkboxInput(inputId = "densityNorm", label = "Density Norm", value = FALSE),
                      checkboxInput(inputId = "MAnorm", label = "MA Norm", value = FALSE),
                      actionButton(inputId = "prev6", label = 'Prev',icon = icon("arrow-left")),
                      actionButton(inputId = "next8", label = 'Next',icon = icon("arrow-right"))
                    ),
                    mainPanel(
                      plotOutput("boxplotNorm",height = "100%", width = "100%"),
                      plotOutput("densityNorm",height = "100%", width = "100%"),
                      plotOutput("MAnorm",height = "100%", width = "100%")
                    )
                  )
                ),
                tabPanel("Normalized array correlation",
                  sidebarLayout(
                    sidebarPanel(
                      width = 2,
                      tags$b("Normalized array correlation"),
                      checkboxInput(inputId = "correlNorm", label = "Correlation Norm",  value = FALSE),
                      checkboxInput(inputId = "PCAnorm", label = "PCA Norm", value = FALSE),
                      checkboxInput(inputId = "clusterNorm", label = "Cluster Norm", value = FALSE),
                      actionButton(inputId = "prev7", label = 'Prev',icon = icon("arrow-left")),
                      actionButton(inputId = "next9", label = 'Next',icon = icon("arrow-right"))
                    ),
                    mainPanel(
                      plotOutput("correlNorm", height = "100%", width = "100%"),
                      plotOutput("PCAnorm",height = "100%", width = "100%"),
                      plotOutput("clusterNorm",height = "100%", width = "100%")
                    )
                  ) 
                )
              )
             ),
             tabPanel("Differential Analysis",
              sidebarLayout(
                sidebarPanel(
                 # width = 2,
                  uiOutput("contrasts"),
                  tags$b("Draw histograms"),
                  #checkboxInput(inputId = "pvalHist", label = "Plot P-Value histogram for each comparison", value = FALSE),
                  #checkboxInput(inputId = "logFCHist", label = "Plot adapted Fold Change historams for each comparison", value = FALSE),
                  actionButton(inputId = "startStat",label = "Calculate statistics",icon = icon("refresh")),
                  tags$br(),
                  tags$br(),
                  actionButton(inputId = "prev8",label = 'Prev',icon = icon("arrow-left")),
                  actionButton(inputId = "next10",label = 'Next',icon = icon("arrow-right"))
                ),
                mainPanel(
                  uiOutput("statOutput"),
                )
              )
             ),
             tabPanel("GO Enrichment",
             sidebarLayout(
               sidebarPanel(
                 width = 2,
                 tags$b("P-value cuttoff to create a list of significant genes"),
                 numericInput(inputId = "Pvalcutoff",label = "P value cutoff for significant genes",value = 0.05,min = 0,max = 2,step = 0.005),
                 tags$hr(),
                 tags$b("GO Enrichment Options"),
                 numericInput(inputId = "gopvalCutoff",label = "P value cutoff for enrichment test report",value = 0.05,min = 0,max = 1,step = 0.005),
                 selectInput(inputId = "ont",label = "Select ontology",choices = c("ALL","BP","MF","CC")),
                 actionButton(inputId = "goAnalysis",label = "Start GO enrichment"),
                 tags$br(),
                 tags$br(),
                 actionButton(inputId = "prev9",label = 'Prev',icon = icon("arrow-left")),
                 actionButton(inputId = "next11",label = 'Next',icon = icon("arrow-right"))
               ),
               mainPanel(
                 uiOutput("goOutput")
               )
             )
             ),
             tabPanel("Pathway Analysis",
                      sidebarLayout(
                        sidebarPanel(
                          uiOutput("pathwayControls")
                        ),
                        mainPanel(
                          uiOutput("pathwayOutput")
                        )
                      )
                      ),
             tabPanel("Network Analysis",
                      sidebarLayout(
                        sidebarPanel(),
                        mainPanel()
                      )
                      ),
             tabPanel("Downloads",
                      sidebarLayout(
                        sidebarPanel(
                          checkboxInput(inputId = "rawQCplot_download", label = "Download all selected plots", value = T),
                          checkboxInput(inputId = "Normtable_download", label = "Download all generated tables", value = T)
                        ),
                        mainPanel()
                      )
                      )
  ),
  tabPanel("Documentation",
           tags$h3("TODO")
           )
)