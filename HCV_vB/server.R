server <- function(input, output, session) {

  # user session$userData to store user data that will be needed throughout
  # the Shiny application
  #session$userData$email <- 'tycho.brahe@tychobra.com'
  session$userData$db_trigger <- reactiveVal(0)

  callModule(parameters_module,id="param_table")
  callModule(fit_module,id="plot_incidence")
  callModule(fit_to_data_module,id="plot_model")
  callModule(simple_intervention_module,id="sim_model")
  daty <- reactive ({NULL})


}
