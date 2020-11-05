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
      numericInput(inputId = "Duration",
                   label = "Duration of the event (in minutes) ",
                   value = 90,
                   min=0),
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
                   label="Will the participants be required to wer any mask",
                   choices = c("no" = 0,
                               "yes" = 1),
                   selected = 0,
                   inline = TRUE),
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
  
  observe({
    x <- input$country
    
  })
  
  
  #output$probs <- renderDataTable({
  #output$region = renderUI({
  #  country_regions = dplyr::filter(regions, country == input$country)
  #  selectInput('region2', 'What region do you live in?', country_regions$region)
  #})
  
  
  dataInput <- reactive({
    testMethod()
    
  })
  
  
  


 
}


shinyApp(ui, server)
