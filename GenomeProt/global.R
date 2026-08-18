suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(shiny)
  library(shinyjs)
  library(shinydashboard)
  library(reticulate)
})

options(shiny.maxRequestSize = 21474836480) # 20 GB

options(scipen = 999)