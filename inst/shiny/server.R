#' Server side Shiny App
#'
#' Identify rows in a data frame with high sclaed dispersion. This rule was
#' taken from the Seurat package.
#'
#' @param input standard input for shiny app
#' @param output standard output for shiny app
#'
#' @importFrom plotly plot_ly
#' @importFrom ggplot2 ggplot
#'

server <- function(input, output) {

   observeEvent(input$button, { 
     
     dataset <- reactive({
       data <- get_data(input$data_dir)
       return(data)
       })
     
     output$table <- renderDataTable({
       dataset()[["initial_metadata"]][, input$columns]
  })
     
     
     output$pca_plotting_color <- renderUI({
       selectInput(inputId = "pca_plotting_color",
                   label = "Plotting Color Labels",
                   c("Default", output$metadata_columns))
     })
     
     
     output$metadata_display <- renderUI({
       selectInput("columns", "Columns",
                   output$metadata_columns, 
                   multiple = TRUE)
     })
     
     
     
  pca <- reactive({
    all_norm <- sm(normalize_expt(dataset(),
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
                     color = ~ as.factor(dataset()[["initial_metadata"]][, input$pca_plotting_color]),
                     colors = brewer.pal(length(unique(dataset()[["initial_metadata"]][, input$pca_plotting_color])), "RdYlBu"))
      
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
                        color = as.factor(dataset()[["initial_metadata"]][, input$pca_plotting_color]))) +
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
                        color = as.factor(dataset()[["initial_metadata"]][, input$pca_plotting_color]))) +
        geom_point() +
        theme_classic() +
        xlab(paste0(input$xaxis2, ": ", pca()$prop_var[as.numeric(substrRight(input$xaxis2, 1))], "% variance")) +
        ylab(paste0(input$yaxis2, ": ", pca()$prop_var[as.numeric(substrRight(input$yaxis2, 1))], "% variance")) +
        labs(colour = input$pca_plotting_color)
      fig
    }
  })
   })
}

#EOF