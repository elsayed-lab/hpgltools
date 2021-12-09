#source("./expt_app_helpers.R")

#' UI side Shiny App
#'
#'
#' @import shiny 
#'

ui <- fluidPage(
  # App title ----
  headerPanel("Summarized Experiment Plots"),
  #h4(paste("Notes:" , dataset$notes)),
  br(),
  
  #fluidRow(
   # splitLayout(cellWidths = c("50%", "50%"),
    #            textInput(inputId = "data_dir", label = "Data Directory To Use:", 
     #                     value = "/mnt/cbcb/fs00_reesyxan/cbcb-lab/nelsayed/scratch/reesyxan/mmusculus_iprgc_2019/expt.rda"),
      #          actionButton("button", "Load Data"))),
  
  textInput(inputId = "data_dir", label = "Data Directory To Use:", value = "/mnt/cbcb/fs00_reesyxan/cbcb-lab/nelsayed/scratch/reesyxan/mmusculus_iprgc_2019/expt.rda"),
  actionButton("button", "Load Data"),
  
  #h4(paste("Dataset Found:" , dataset$notes)),
  #h4(paste("Notes:" , dataset$notes)),
  
  # Sidebar panel for inputs ----
  sidebarPanel(
    width = 2,
    # Input: Selector for variable to plot against mpg ----
    # selectInput(inputId = "data",
    #            label = "Dataset:",
    #           c("TMRC3" = "TMRC3",
    #            "TMRC2" = "TMRC2"
    #            )),
    
    h5("Normalization Options"),
    selectInput(inputId = "convert",
                label = "Convert:",
                c("FALSE" = "none",
                  "cpm" = "cpm",
                  "rpkm" = "rpkm",
                  "cbcb" = "cbcb"
                )),
    selectInput(inputId = "normalization",
                label = "Normalization:",
                c("FALSE" = "FALSE",
                  "quant" = "quant",
                  "cbcb" = "cbcb"
                )),
    selectInput(inputId = "transform",
                label = "Transform:",
                c("FALSE" = "FALSE",
                  "log2" = "log2"
                )),
    
    checkboxInput(inputId = "filter",
                  label = "Filter",
                  TRUE),
    br(),
    h5("Plotting Options"),
    uiOutput('pca_plotting_color')
  ),
  
  mainPanel(
    width = 10,
    # Output: Tabset w/ plot, summary, and table ----
    tabsetPanel(type = "tabs",
                tabPanel("Table",
                         h4("Meta Data to Display"),
                         uiOutput("metadata_display"),
                         verbatimTextOutput("dfStr"),
                         br(),
                         dataTableOutput("table")),
                tabPanel("2D PCA Plots",
                         br(),
                         fluidRow(
                           splitLayout(cellWidths = c("50%", "50%"),
                                       h4("PCA Plot 1"),
                                       h4("PCA Plot 2"))),
                         fluidRow(
                           splitLayout(cellWidths = c("50%", "50%"),
                                       selectInput(inputId = "xaxis1",
                                                   label = "X-Axis",
                                                   c(paste("PC", 1:10))),
                                       selectInput(inputId = "xaxis2",
                                                   label = "X-Axis",
                                                   c(paste("PC", 1:10)),
                                                   selected = "PC 2")),
                           tags$head(tags$style(HTML("
                              .shiny-split-layout > div {
                                overflow: visible;
                              }
                              ")))),
                         fluidRow(
                           splitLayout(cellWidths = c("50%", "50%"),
                                       selectInput(inputId = "yaxis1",
                                                   label = "Y-Axis",
                                                   c(paste("PC", 1:10)),
                                                   selected = "PC 2"),
                                       selectInput(inputId = "yaxis2",
                                                   label = "Y-Axis",
                                                   c(paste("PC", 1:10)),
                                                   selected = "PC 3"))),
                         br(),
                         
                         fluidRow(
                           splitLayout(cellWidths = c("48%", "52%"),
                                       plotOutput("PCA1"), plotOutput("PCA2"))),
                         br(),
                         br(),
                ),
                
                tabPanel("FOR NAJIB: 3D PCA Plot",
                         # selectInput(inputId = "plotting_color",
                         #            label = "Plotting Color Labels",
                         #           names(dataset[["initial_metadata"]]),
                         #          selected = "typeofcells"),
                         br(),
                         plotlyOutput("PCA"))
                )
    )
  )

#EOF
