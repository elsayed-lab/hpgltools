dataset <- get_data(data = "TMRC3")

#' Set up the Main siny page.
#'
#' @import shiny
#' @import plotly
ui <- fluidPage(
    # App title ----
    headerPanel("TMRC3 Plots"),
    h4(paste("Notes:" , dataset$notes)),
    br(),

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
        selectInput(inputId = "pca_plotting_color",
                    label = "Plotting Color Labels",
                    c("Default", names(dataset[["initial_metadata"]])))
    ),

    mainPanel(
        width = 10,
        # Output: Tabset w/ plot, summary, and table ----
        tabsetPanel(type = "tabs",
                    tabPanel("Table",
                             h4("Meta Data to Display"),
                             selectInput("columns", "Columns",
                                         names(dataset[["initial_metadata"]]), multiple = TRUE),
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
                    )))

server <- function(input, output) {
  output$table <- renderDataTable({
    dataset[["initial_metadata"]][, input$columns]
  })

  pca <- reactive({
    all_norm <- sm(normalize_expt(dataset,
                                  transform = input$transform,
                                  norm = input$normalization,
                                  convert = input$convert,
                                  filter = input$filter))
    pca <- plot_pca(all_norm, plot_labels=FALSE)
    return(pca)
  })

  output$PCA <- renderPlotly({
    if(input$pca_plotting_color == "Default"){
      fig <- plot_ly(data = pca()$table,
                     x = ~ pc_1,
                     y = ~ pc_2,
                     z = ~ pc_3)

      fig <- fig %>% add_markers()
      fig <- fig %>% layout(scene = list(
                                xaxis = list(title = paste0('PC1: ', pca()$prop_var[1], "% variance")),
                                yaxis = list(title = paste0('PC2: ', pca()$prop_var[2], "% variance")),
                                zaxis = list(title = paste0('PC3: ', pca()$prop_var[3], "% variance"))))
      fig

    } else {
      fig <- plot_ly(data = pca()$table,
                     x = ~ pc_1,
                     y = ~ pc_2,
                     z = ~ pc_3,
                     color = ~ as.factor(dataset[["initial_metadata"]][, input$pca_plotting_color]),
                     colors = brewer.pal(length(unique(dataset[["initial_metadata"]][, input$pca_plotting_color])), "RdYlBu"))

      fig <- fig %>% add_markers()
      fig <- fig %>% layout(scene = list(
                                xaxis = list(title = paste0('PC1: ', pca()$prop_var[1], "% variance")),
                                yaxis = list(title = paste0('PC2: ', pca()$prop_var[2], "% variance")),
                                zaxis = list(title = paste0('PC3: ', pca()$prop_var[3], "% variance"))))
      fig
    }
  })

  output$PCA1 <- renderPlot({
    if (input$pca_plotting_color == "Default") {
      fig <- ggplot(data = pca()$table,
                    aes(x = pca()$table[, paste0("pc_", substrRight(input$xaxis1, 1))],
                        y = pca()$table[, paste0("pc_", substrRight(input$yaxis1, 1))])) +
        geom_point(color = "Blue") +
        theme_classic() +
        xlab(paste0(input$xaxis1, ": ", pca()$prop_var[as.numeric(substrRight(input$xaxis1, 1))], "% variance")) +
        ylab(paste0(input$yaxis1, ": ", pca()$prop_var[as.numeric(substrRight(input$yaxis1, 1))], "% variance")) +
        theme(legend.position="none")
      fig
    } else {
      fig <- ggplot(data = pca()$table,
                    aes(x = pca()$table[, paste0("pc_", substrRight(input$xaxis1, 1))],
                        y = pca()$table[, paste0("pc_", substrRight(input$yaxis1, 1))],
                        color = as.factor(dataset[["initial_metadata"]][, input$pca_plotting_color]))) +
        geom_point() +
        theme_classic() +
        xlab(paste0(input$xaxis1, ": ", pca()$prop_var[as.numeric(substrRight(input$xaxis1, 1))], "% variance")) +
        ylab(paste0(input$yaxis1, ": ", pca()$prop_var[as.numeric(substrRight(input$yaxis1, 1))], "% variance")) +
        theme(legend.position = "none")
      fig
    }
  })

  output$PCA2 <- renderPlot({
    if (input$pca_plotting_color == "Default") {
      fig <- ggplot(data = pca()$table,
                    aes(x = pca()$table[, paste0("pc_", substrRight(input$xaxis2, 1))],
                        y = pca()$table[, paste0("pc_", substrRight(input$yaxis2, 1))])) +
        geom_point(color = "Blue") +
        theme_classic() +
        xlab(paste0(input$xaxis2, ": ", pca()$prop_var[as.numeric(substrRight(input$xaxis2, 1))], "% variance")) +
        ylab(paste0(input$yaxis2, ": ", pca()$prop_var[as.numeric(substrRight(input$yaxis2, 1))], "% variance")) +
        labs(colour = input$pca_plotting_color)
      fig
    } else {
      fig <- ggplot(data = pca()$table,
                    aes(x = pca()$table[, paste0("pc_", substrRight(input$xaxis2, 1))],
                        y = pca()$table[, paste0("pc_", substrRight(input$yaxis2, 1))],
                        color = as.factor(dataset[["initial_metadata"]][, input$pca_plotting_color]))) +
        geom_point() +
        theme_classic() +
        xlab(paste0(input$xaxis2, ": ", pca()$prop_var[as.numeric(substrRight(input$xaxis2, 1))], "% variance")) +
        ylab(paste0(input$yaxis2, ": ", pca()$prop_var[as.numeric(substrRight(input$yaxis2, 1))], "% variance")) +
        labs(colour = input$pca_plotting_color)
      fig
    }
  })
}

# EOF
