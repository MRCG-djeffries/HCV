library(shinyjs)
library(shinycssloaders)
library(shinyFeedback)
library(ggplot2)
library(RSQLite)
library(shinyWidgets)
library(dbplyr)
library(shinydashboard)
library(cowplot)
ui <- fluidPage(
  shinyFeedback::useShinyFeedback(),
  shinyjs::useShinyjs(),
  # Application Title
  tabsetPanel(type = "tabs",
  tabPanel("Model parameters",
    parameters_module_ui("param_table")
  ),
  tabPanel("Epidemic data",fit_module_ui("plot_incidence")),
  tabPanel("Fit model",fit_to_data_module_ui("plot_model")),
  tabPanel("Intervention model",simple_intervention_module_ui("sim_model"))
)
)
      
