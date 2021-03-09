fit_module_ui <-function(id){
  ns = NS(id)
  
  ui <- fluidPage(
    titlePanel(title=h4("Epidemic data for fitting", align="center")),
    tags$h5(id="head1","1. Incidence scaling due to drug market activity fitted as a multiplier for probability of infection"),
    tags$style(HTML("#head1{color: red;}")),
    fluidRow(
      column(width = 12,
             mainPanel(plotOutput(ns("plot_incidence"))),
      )),
    tags$h5(id="head2","2. Fitting points"),
    tags$style(HTML("#head2{color: red;}")),
    tags$ul(
      tags$li("50% of current PWID are chronically infected with HCV in 2015"), 
      tags$li("45% of current and former chronically infected with HCV are in stages F0 & F1 in 2015"), 
    ),
    tags$h5(id="head3","3. Scaling to a population"),
    tags$style(HTML("#head3{color: red;}")),
    tags$h5(id="The model assumes a standardised population of 1,000 subjects"),
    tags$ul(
      tags$li("In 2015 there were 40,000 HCV infected current PWID"), 
      tags$li("In 2015 there were 144,000 HCV infected former PWID"), 
    ),
  )
}

fit_module <- function(input, output, session) {

  output$plot_incidence<-renderPlot({
    plot_inci_trend()
})

}