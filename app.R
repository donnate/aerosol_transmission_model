# Define UI for app that draws a histogram ----
library("shiny")
library("tidyverse")
setwd("~/Dropbox/aerosol_transmission_model/")
source("individual_probabilities.R")
source("helper_functions.R")
source("aerosol_functions.R")
source("preprocessing_functions.R")

countries = read.csv("population_by_country_2020.csv")
SENSITIVITY = c(0,0,0.019,0.0327,0.560,0.653,0.718,0.746,0.737,0.718,0.7,0.68,0.662,0.644,0.625)
RELATIVE_INFECTIOUSNESS = c(0,0.01,0.05,0.2,0.6,0.88,0.98,1,1,1,0.95,0.8,0.4,0.2,0.1,0.01)
TAU = 0.1

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
      numericInput(inputId = "n",
                   label = "Total of people present",
                   value = 0,
                   min=0),
      radioButtons(inputId = "activity",
                   label="What activity will the participants be performing?",
                   choices = c("singing" = 0,
                               "eating" = 1,
                               "walking" = 2,
                               "light_exercise" = 3,
                               "strenuous_exercise" = 4),
                   selected = 0,
                   inline = TRUE),
      radioButtons(inputId = "mask",
                   label="Will the participants be required to wear any mask",
                   choices = c("no" = 0,
                               "yes" = 1),
                   selected = 0,
                   inline = TRUE),
      radioButtons(inputId = "mixing",
                   label="Are the people going to be mixing in that event?",
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
      
      
      # Horizontal line ----
      tags$hr(),
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Data file ----
      tableOutput("contents")
      
    )
    
  )
)


# Define server logic required to draw a histogram ----
server <- function(input, output, session) {

  dataInput <- reactive({
    
    
    #### Step 1: load the participants data + location data. In the absence of data, assume homogeneity.
    ####### Step 1.a: Load participants
    input = list(country = "France",
                 date_event= "2020-09-01",
                 duration = 90,
                 unit = "m",
                 length=100,
                 width = 100,
                 height = 5,
                 pressure= 0.95,
                 yemperature = 20,
                 ventilation_out = 0.7,
                 decay_rate = 0.62,
                 deposition = 0.3,
                 controls = 0,
                 n=60,
                 activity = 0,
                 mask = 0,
                 mixing = 0
                 )
    df <- read.csv(input$file1$datapath,
                   header = TRUE,
                   sep = ",")
    #df <- read.csv("~/Dropbox/mock_data.csv",
    #                             header = TRUE,
    #                               sep = ",")
    ####### Enrich the dataset by computing the prevalence of the virus up
    ####### to 14 days before the event
    prevalence =  rep(0.005,  length(SENSITIVITY)) #extract_prevalence()
    
    ###### Step 1.b: compute room parameters for aerosolization
    volume = extract_volume(input$length, input$width, input$height)
    first_order_loss_rate = extract_first_order(input$ventilation,
                                                input$control,
                                                input$decay,
                                                input$deposition)
    ventilation_rate_per_person = extract_ventilation_rate_per_person(volume, 
                                                                      input$ventilation,
                                                                      input$control,
                                                                      input$n)
    
    ##########################################
    #### Step 2: compute probability of being infectious, hospitalized,and dying
    
   df$p =  unlist(sapply(1:5, function(x){
      compute_infectiousness_probability(sensitivity = SENSITIVITY,
                                         prevalence=prevalence,
                                         df$profession[x],
                                         df$high_risk_contact[x],
                                         df$mask_wearing[x],
                                         df$nb_people_hh[x])
      }))
    
    df$p_hosp =  unlist(sapply(1:5, function(x){
      compute_hospitalization_probability(c(df$age[x], df$Pregnant[x],
                                            df$Chronic_Renal_Insufficiency[x],
                                            df$Diabetes[x], df$Immunosuppression[x],
                                            df$COPD[x], df$Obesity[x], 
                                            df$Hypertension[x],
                                            df$Tobacco[x], df$Cardiovascular_Disease[x],
                                            df$Asthma[x], df$Sex[x]))
    }))
    
    df = df %>% rowwise() %>% mutate(p_death = death_probability(age, Pregnant,
                                                                 Chronic_Renal_Insufficiency,
                                                                 Diabetes, Immunosuppression, COPD,
                                                                 Obesity, Hypertension, Tobacco,
                                                                 Cardiovascular_Disease,Asthma, Gender))

    ##########################################
    nb_infective_people = rep(0,B)
    nb_infections = rep(0,B)
    p = rep(0,B)
    nb_hospitalizations = rep(0,B)
    nb_deaths = rep(0,B)
    #### Step 3: compute probability of infecting people and adverse outcomes using MCMC simulations
    for (b in 1:B){
      ####### draw infected people (their infectivity bucket)
      Z = apply(pi,axis= 1, function(x){dirichlet(1,x)})
      Z = Z[(Z>1)]  #### keep infected people only
      nb_infective_people[b] = length(Z)
      
      ####### Step 3.a Direct contacts
      contacts  = rnorm(length(Z), mu = input$mu, sd = input$sd)
      nb_infections[b] = sapply(1:length(Z),
                            function(x){rbinom(contacts[x], TAU *  RELATIVE_INFECTIOUSNESS[x])})
      
      ####### Step 3.b Aerosolization
      quanta_emission_rate <- compute_quanta_emission_rate(input$quanta_exhalation_rate,
                                                           input$mask_efficiency ,
                                                           input$prop_mask,
                                                           nb_infective_people[b])
      quanta_concentation <- compute_quanta_concentation(quanta_emission_rate,
                                                         first_order_loss_rate,
                                                         volume,  
                                                         input$duration)
      quanta_inhaled_per_person <- compute_quanta_inhaled_per_person(quanta_concentration,
                                                                     breathing_rate,
                                                                     input$duration,
                                                                     input$inhalation_mask_efficiency,
                                                                     input$prop_mask)
      
      p[b] = 1-exp(-quanta_inhaled_per_person)
      nb_infections[b] =  nb_infections[b] + rbinom(N - nb_infections[b] - nb_infective_people[b],
                                                    p[b])
      
    }
    
    #### Step 4: Compute (by MCMC simulation) the number of people with adverse outcomes
    for (b in 1:B){
      ###### Sample from the list
      nb_hospitalizations[b] = sum(sapply(sample(df$p_hosp, nb_infections[b]), function(x){rbinom(1,x)}))
      nb_deaths[b] = sum(sapply(sample(df$p_death, nb_infections[b]), function(x){rbinom(1,x)}))
    }
    return(list(nb_infections=nb_infections, nb_deaths=nb_deaths, nb_infections=nb_infections))
  })
  
  output$contents <- renderTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file1)
    
    df <- read.csv(input$file1$datapath,
                   header = TRUE,
                   sep = ",")
    return(head(df))
    
    
  })
  
  
  
  
  
}


shinyApp(ui, server)
