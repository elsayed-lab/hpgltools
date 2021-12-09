#!/usr/bin/env Rscript

devtools::load_all("../../")
library(shiny)
library(hpgltools)
library(ggplot2)
library(plotly)
library(dplyr)
library(RColorBrewer)
library(stringr)

ui <- hpgltools:::ui
server <- hpgltools:::server
shinyApp(ui = ui, server = server)
