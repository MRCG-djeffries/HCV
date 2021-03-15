fit_to_data_module_ui <- function(id) {
  ns <- NS(id)
  
  # tagList(
  dashboardPage(
    dashboardHeader(title = "choose parameter set",titleWidth = 250),
    dashboardSidebar(width = 250, shinyjs::useShinyjs(),
                     actionBttn(ns("model_fit1"), label="Fit model 1",size="xs"),
                     actionBttn(ns("model_fit2"), label="Fit model 2",size="xs"),
                     tags$hr(),
                     tags$h5("These are live model fits dynamically running the compartment model for each parameter set"),
                     tags$hr(),
                     tags$h4(id="wordy2","Note the parameter sets were obtained by minimising the errors to the two fitting points on the epidemic data page (50% and 45%). This is an intenisve multi-variate optimization process and is not suitable for server sided applications. Additionally there can be
                             subtle differences in fit, depending which parameters are actually optimized for",
                             tags$style(HTML("#wordy2{color: yellow;}"))),
                     actionBttn(ns("opt_detail"), label="Optimization details",size="xs"),
                     tags$hr(),
                     actionBttn(ns("reset_input"), label="Reset",size="xs")
    ),
    dashboardBody(shinyjs::useShinyjs(),tags$head(tags$style(HTML('
      .content-wrapper {
        background-color: #fff;
      }
    '
    ))), 
    h5(id="wordy",htmlOutput(ns("plot_text"), container = span)),tags$style(HTML("#wordy{color: red;}")),
    fluidRow(
      column(width = 12,
             mainPanel(plotOutput(ns("plot_model"),height="600px") %>% withSpinner() ),
      )
    )
    )
  )
  
}

fit_to_data_module <- function(input, output, session) {
  output$plot_model =  renderPlot(NULL) # so no spinner until button pushed
                     
  observeEvent(input$model_fit1, {
    output$plot_text <- renderUI({HTML(paste0("Fit for 3 parameters: probability of transmission, 
    F0 to F1 rate for current and former PWIDs.","<br>", 
    "The fitted rates of F0 to F1 for current and former are 0.077 and 0.078","<br>",
    "Compared to the literature estimated rates of 0.106 and 0.116","<br>",
    "Note in the parameter database, there is no value for probability of transmission (estimated as 13.9%)","<br>",
    "GENERALLY AT LEAST ONE PARAMETER IS USED TO CALIBRATE THE MODEL TO KNOWN EPIDEMIC DATA"
    ))})
    output$plot_model =  renderPlot(runmodel(2))
    ns =session$ns
    output$ui =renderUI(   plotOutput(ns("plot_model"))  )

  })
  observeEvent(input$reset_input,{
    output$plot_model =  renderPlot(NULL)
    output$plot_text = renderUI({HTML("The model is fit to the epideic data, by varying 3 (fit 1) or 5 (fit 2) parameters")})
  })
  observeEvent(input$model_fit2, {
    output$plot_text <- renderUI({HTML(paste0("Fit for 5 parameters: probability of transmission, 
    F0 to F1 rate for current and former PWIDs and relapse rate, injecting career duration","<br>", 
    "The fitted rates of F0 to F1 for current and former are 0.104 and 0.078","<br>",
    "Compared to the literature estimated rates of 0.106 and 0.116","<br>",
    "The fitted rates of relapse rate, injecting career duration are 0.105 and 18.39 years","<br>",
    "Compared to the literature estimated rates of 0.027 and 17 years","<br>",
    "Note in the parameter database, there is no value for probability of transmission (estimated as 10.4%)","<br>",
     "GENERALLY AT LEAST ONE PARAMETER IS USED TO CALIBRATE THE MODEL TO KNOWN EPIDEMIC DATA"
    ))})
    output$plot_model =  renderPlot(runmodel(1))
    ns =session$ns
    output$ui =renderUI(   plotOutput(ns("plot_model"))  )
  })
  output$plot_text = renderUI({HTML("The model is fit to the epideic data, by varying 3 (fit 1) or 5 (fit 2) parameters")})
  #output$plot_text <- renderText({"The model is fit to the epideic data, by varying 3 (fit 1) or 5 (fit 2) parameters"})
  observeEvent(input$opt_detail, {
    ns =session$ns
    #addResourcePath("picy", "~/HCV/www")
    
    showModal(modalDialog(
      size="l",
      HTML('<img src="optofig.PNG">')))
   })
  
}
