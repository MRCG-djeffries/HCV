simple_intervention_module_ui <- function(id) {
  ns <- NS(id)
 
  # tagList(
  dashboardPage(
    dashboardHeader(title = "model interventions",titleWidth = 250),
    dashboardSidebar(width = 250, shinyjs::useShinyjs(),    chooseSliderSkin("Flat", color = "blue"),
                     tags$h5("Simulate an intervention from 2020"),
                     sliderInput(ns("cessation"), "Decrease injecting career by (%):",
                                 min = 0, max = 50,
                                 value = 0, step = 1),
                     sliderInput(ns("relapse"), "Reduce use relapse by (%):",
                                 min = 0, max = 50,
                                 value = 0, step = 1),
                     tags$br(),
                     actionBttn(ns("run_sim"), label="Run simulation",size="xs"),
                     tags$hr(),
                     actionBttn(ns("reset_input1"), label="Reset",size="xs")
    ),
    dashboardBody(shinyjs::useShinyjs(),tags$head(tags$style(HTML('
      .content-wrapper {
        background-color: #fff;
      }
    '
    ))), 
    h5(id="wordy",htmlOutput(ns("plot_simtext"), container = span)),tags$style(HTML("#wordy{color: red;}")),
    fluidRow(
      column(width = 12,
             mainPanel(plotOutput(ns("plot_simmodel"),height="600px") %>% withSpinner() ),
      )
    )
    )
  )
  
}

simple_intervention_module <- function(input, output, session) {
 
   output$plot_simmodel =  renderPlot(NULL) # so no spinner until button pushed
  # show("plot_simmodel")
  # data =list(cess=-1,rela=-1)

  data = eventReactive(input$run_sim,{
    list(cess=input$cessation,rela=input$relapse)
  })
  observeEvent(input$reset_input1,{
    reset("cessation")
    reset("relapse")
    output$plot_simmodel =  renderPlot(NULL)
  })  
  observeEvent(input$run_sim,{
    output$plot_simmodel =  renderPlot(
        runmodel_intervention(data()$cess,data()$rela)
      )
    ns =session$ns
    output$ui =renderUI(   plotOutput(ns("plot_simmodel"))  )
  })
  
  output$plot_simtext = renderUI({HTML("Move sliders to simulate intervention (3 parameter model fitting) and click button to run")})
  
 
  
  
}
