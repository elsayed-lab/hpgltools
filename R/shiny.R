## A place to put shiny things

#' make a slideable threshold for cutting off significant genes in a DE table.
#'
#' @import shiny
slide_de_threshold <- function(de_table, contrast = 1, lfc = 1.0, p = 0.05, according_to = "deseq") {
  df <- de_table[["data"]][[contrast]][, c("deseq_logfc", "deseq_adjp")]
  ui <- fluidPage(
      titlePanel(title = h4("Test", align = "center")),
      sidebarPanel(
          sliderInput("lfc", "logFC cutoff: ", 1.0, min = 0.0, max = 12.0, step = 0.1),
          sliderInput("p", "p-value cutoff: ", 0.05, min = 0.0, max = 0.2, step = 0.001)),
      mainPanel(
          plotOutput("volcano"),
          tableOutput("volcano_summary")
      ))

  server <- function(input, output) {
    myProcess <- reactiveVal(NULL)
    observeEvent(input$stopSession, {
      cat(sprintf("Closing session %s\n", session$token))
      session$close()
    })
    dat <- reactive({
      p_idx <- df[["deseq_adjp"]] <= input[["p"]]
      lfc_idx <- abs(df[["deseq_logfc"]]) >= input[["lfc"]]
      both_idx <- p_idx & lfc_idx
      test <- df[both_idx, ]
      df[["lp"]] <- -1.0 * log10(df[["deseq_adjp"]])
      df[["lf_sig"]] <- "insig"
      df[lfc_idx, "lf_sig"] <- "sig"
      df[["p_sig"]] <- "insig"

      df[p_idx, "p_sig"] <- "sig"
      df[["both_sig"]] <- "insig"
      df[both_idx, "both_sig"] <- "sig"

      df[["lf_sig"]] <- as.factor(df[["lf_sig"]])
      df[["p_sig"]] <- as.factor(df[["p_sig"]])
      df[["both_sig"]] <- as.factor(df[["both_sig"]])
      ret <- list(
          "p" = p_idx,
          "lfc" = lfc_idx,
          "both" = both_idx,
          "df" = df)
      return(ret)
    })
    output$volcano <- renderPlot({
      stuff <- dat()
      df <- stuff[["df"]]
      ggplot(df, aes_string(x = "deseq_logfc", y = "lp", color = "both_sig")) +
        geom_point()
    }, height = 400,width = 600)
    output$volcano_summary <- renderTable({
      stuff <- dat()
      summary_df <- data.frame("lfc" = c(0, 0),
                               "p" = c(0, 0),
                               "both" = c(0, 0))
      rownames(summary_df) <- c("sig", "insig")
      summary_df["sig", "lfc"] <- sum(stuff[["lfc"]])
      summary_df["insig", "lfc"] <- sum(!stuff[["lfc"]])
      summary_df["sig", "p"] <- sum(stuff[["p"]])
      summary_df["insig", "p"] <- sum(!stuff[["p"]])
      summary_df["sig", "both"] <- sum(stuff[["both"]])
      summary_df["insig", "both"] <- sum(!stuff[["both"]])
    })
  }
  app <- shinyApp(ui, server)
  return(app)
}
