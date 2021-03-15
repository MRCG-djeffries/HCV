parameters_module_ui <- function(id) {
  ns <- NS(id)

  # tagList(
  dashboardPage(title="HCV model",
    dashboardHeader(title = "choose parameter set",titleWidth = 250),
    dashboardSidebar(width = 250, shinyjs::useShinyjs(),
                     selectInput(ns("category"), "Category:",
                                 c("Transition" = "Transition",
                                   "Liver related mortality" = "Liver related mortality",
                                   "Non liver related mortality" = "Non liver related mortality",
                                   "Drug use"="Drug use",
                                   "One off costs"="One off costs",
                                   "Costs associated with treatment"="Costs associated with treatment",
                                   "Other annual costs"="Other annual costs",
                                   "Health utilities"="Health utilities"
                                   ),
                                 multiple = TRUE,size=8,selectize=FALSE),
                     actionBttn(ns("reset_cat"), label="Reset",size="xs"),
                     tags$hr()
                   
    ),
    dashboardBody(shinyjs::useShinyjs(),tags$head(tags$style(HTML('
      .content-wrapper {
        background-color: #fff;
      }
    '
    ))), 
    fluidRow(
      column(
        width = 12,
        title = "Parameters in model",
        DTOutput(ns("param_table")) %>%
          withSpinner(),
        tags$br(),
        tags$br()
      )
    )
   )
  )
 
}

parameters_module <- function(input, output, session) {

 
  dats <- reactive({
    session$userData$db_trigger()

    out <- NULL
    tryCatch({
      out <- conn %>%
        tbl('parameter_data') %>%
        collect()  %>%
         arrange(Vnum)

    }, error = function(err) {

      print(err)
      showToast("error", "Database Connection Error")

    })
    format(out$value, big.mark=',', scientific=FALSE) 
    
    out

  })

  observeEvent(input$reset_cat, {
    reset("category")
  })
 
  output$param_table <- renderDT({
    if (length(input$category)>0){
      out=dats()%>% filter(Variable %in% input$category)
      # for ( i in 1 : length(input$category)){
      #     out =out %>% filter(Variable == input$category[i])
      # }
    headerCallback = JS(
      "function( thead, data, start, end, display ) {
      $(thead).closest('thead').find('th').eq(1).css('color', 'white');
      $(thead).closest('thead').find('th').eq(2).css('color', 'white');
      $(thead).closest('thead').find('th').eq(3).css('color', 'red');
      $(thead).closest('thead').find('th').eq(5).css('color', 'white');
      $(thead).closest('thead').find('th').eq(6).css('color', 'red');
              }"
    )
    datatable(
      out,
      rownames = FALSE,
       colnames = c('Vnum', 'Variable', 'Stratum', 'Parameter', 'Units', 'Type', 'Value', 'Ref'),
      selection = "none",
      class = "compact stripe row-border nowrap",
     
      escape = -1,
      extensions = c("Buttons"),
      options = list(
        scrollX = TRUE,
        dom = 'Blftip',
        pageLength = 15,
        buttons = list(
          list(
            extend = "excel",
            text = "Download",
            title = paste0("parameters-", Sys.Date()),
            exportOptions = list(
              columns = 1:(length(out) - 1)
            )
          )
        ),
        columnDefs = list(
          list(targets = c(0,1,3,4,5,6), orderable = FALSE),
          list(className = 'dt-left', targets = 0:4),
          list(className = 'dt-center', targets = 6)
        ),
        headerCallback = JS(headerCallback),
        initComplete = JS(
          "function(settings, json) {",
          "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff',' text-align': 'left;'});",
          "}")
      )
    )%>%
      formatStyle('Vnum', `text-align` = 'left')%>%
      formatStyle('Value', `text-align` = 'right',color = 'red',fontWeight = 'bold')%>%formatCurrency('Value', '')%>%
      formatStyle('Parameter', `text-align` = 'left',color = 'red',fontWeight = 'bold')
    }
    
  })








}
