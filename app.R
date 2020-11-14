# Define UI for app that draws a histogram ----
library("shiny")
library("tidyverse")

countries = read.csv("population_by_country_2020.csv")
ui <- fluidPage(
  
  # App title ----
  titlePanel("What are the parameters for your event?"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      selectInput(inputId = "country",
                  label = "Which country do you live in?",
                  choices = countries$Country,
                  selected = "United Kingdom"),
      dateInput(inputId = "date_event",
                label ="When will the event occur?", 
                value ="2020-09-01",
                min = NULL,
                max = NULL,
                format = "yyyy-mm-dd",
                startview = "month",
                weekstart = 0,
                language = "en",
                width = NULL,
                autoclose = TRUE,
                datesdisabled = NULL,
                daysofweekdisabled = NULL),
      
      numericInput(inputId = "Duration",
                   label = "Duration of the event (in minutes) ",
                   value = 90,
                   min=0),
      # Horizontal line ----
      tags$hr(),
      radioButtons(inputId = "unit",
                   label="What measurement unit will you be using",
                   choices = c("ft" = 1,
                               "m" = 0),
                   selected = 0,
                   inline = TRUE),
      numericInput(inputId = "Length",
                   label = "Length (in your selected metric). ",
                   value = 1,
                   min=0),
      numericInput(inputId = "Width (in your selected metric).",
                   label = "Width ",
                   value = 1,
                   min=0),
      numericInput(inputId = "Height (in your selected metric).",
                   label = "Height ",
                   value = 1,
                   min=0),
      numericInput(inputId = "Pressure (in atm)",
                   label = "Pressure ",
                   value = 0.95,
                   min=0),
      numericInput(inputId = "Temperature",
                   label = "Temperature (in Celsius)",
                   value = 20,
                   min=0),
      # Horizontal line ----
      tags$hr(),
      numericInput(inputId = "Ventilation_out",
                   label = "Ventilation with outside air",
                   value = 0.7,
                   min=0),
      numericInput(inputId = "Decay_rate",
                   label = "Decay rate of the virus",
                   value = 0.62,
                   min=0),
      numericInput(inputId = "Deposition",
                   label = "Deposition to surfaces",
                   value = 0.3,
                   min=0),
      numericInput(inputId = "controls",
                   label = "Additional control measures",
                   value = 0,
                   min=0),
      numericInput(inputId = "N",
                   label = "Total of people present",
                   value = 0,
                   min=0),
      numericInput(inputId = "N_inf",
                   label = "Nb of infected individuals",
                   value = 0,
                   min=0),
      numericInput(inputId = "p_susc",
                   label =  "Percentage of susceptible people",
                   value = 100.0,
                   min=0, max=100),
      radioButtons(inputId = "activity",
                   label="What activity will the participants be performing?",
                   choices = c("singing" = 0,
                               "eating" = 1,
                               "walking" = 2),
                   selected = 0,
                   inline = TRUE),
      radioButtons(inputId = "mask",
                   label="Will the participants be required to wear any mask",
                   choices = c("no" = 0,
                               "yes" = 1),
                   selected = 0,
                   inline = TRUE),
      
      # Input: Select a file ----
      fileInput("file1", "Upload participants' information (choose CSV File)",
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      
      # Horizontal line ----
      tags$hr(),
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      tabsetPanel(
        )
      
    )

  )
)


# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  
  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  # ------------------ App virtualenv setup (Do not edit) ------------------- #
  
  
  #output$probs <- renderDataTable({
  #output$region = renderUI({
  #  country_regions = dplyr::filter(regions, country == input$country)
  #  selectInput('region2', 'What region do you live in?', country_regions$region)
  #})
  
  dataInput <- reactive({
    volume = extract_volume(input$length, input$width, input$height)
    first_order_loss_rate = extract_first_order(input$ventilation,
                                                input$control,
                                                input$decay,
                                                input$deposition)
    ventilation_rate_per_person = extract_ventilation_rate_per_person(volume, 
                                                                      input$ventilation,
                                                                      input$control,
                                                                      input$n)
    nb_infective_people <-  compute_number_infective_people()
    quanta_emission_rate <- compute_quanta_emission_rate(input$quanta_exhalation_rate,
                                                         input$mask_efficiency ,
                                                         input$prop_mask,
                                                         nb_infective_people)
    quanta_concentation <- compute_quanta_concentation(quanta_emission_rate,
                                                       first_order_loss_rate,
                                                       volume,  
                                                       input$duration)
    quanta_inhaled_per_person <- compute_quanta_inhaled_per_person(quanta_concentration,
                                                                   breathing_rate,
                                                                   input$duration,
                                                                   input$inhalation_mask_efficiency,
                                                                   input$prop_mask)
    
    p = 1-exp(-quanta_inhaled_per_person)
    return(list(p_infection= p, p_hosp = input$hosp_rate *p, p_death = p * input$death_rate))
    
  })
  
  output$distPlot <- renderPlot({
    x =  dataInput()
    
    hist(x, breaks = seq(from=0, to=1, by=0.025), col = "#75AADB", border = "white",
         xlab = "Probability",
         main ="Your probability distribution" )
  })

  
  


 
}


shinyApp(ui, server)
