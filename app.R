# Define UI for app that draws a histogram ----
library("shiny")
library("tidyverse")
library(MCMCpack)
setwd("~/Dropbox/aerosol_transmission_model/")
source("individual_probabilities.R")
source("helper_functions.R")
source("aerosol_functions.R")
source("preprocessing_functions.R")

countries = read.csv("population_by_country_2020.csv")
BREATHING_RATE = 0.7
SENSITIVITY = c(0,0,0.019,0.0327,0.560,0.653,0.718,0.746,0.737,0.718,0.7,0.68,0.662,0.644,0.625)
RELATIVE_INFECTIOUSNESS = c(0, 0.01,0.05,0.2,0.6,0.88,0.98,1,1,1,0.95,0.8,0.4,0.2,0.1,0.01)
MASK_EFFICIENCY = 0.5  ### 50% is the recommended value
MASK_INHALATION_EFFICIENCY =0.3
TAU = 0.2
MU = 20
SD = 4

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
                value ="2020-11-15",
                min = NULL,
                max = "2020-11-29",
                format = "yyyy-mm-dd",
                startview = "month",
                weekstart = 0,
                language = "en",
                width = NULL,
                autoclose = TRUE,
                datesdisabled = NULL,
                daysofweekdisabled = NULL),
      
      numericInput(inputId = "duration",
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
      numericInput(inputId = "length",
                   label = "Length (in your selected metric). ",
                   value = 10,
                   min=0),
      numericInput(inputId = "width",
                   label = "Width (in your selected metric).",
                   value = 10,
                   min=0),
      numericInput(inputId = "height",
                   label = "Height (in your selected metric).",
                   value = 10,
                   min=0),
      numericInput(inputId = "pressure",
                   label = "Pressure (in atm) ",
                   value = 0.95,
                   min=0),
      numericInput(inputId = "temperature",
                   label = "Temperature (in Celsius)",
                   value = 20,
                   min=0),
      numericInput(inputId = "RH",
                   label = "Relative Humidity (from 20 to 70%)",
                   value = 60,
                   min=20,
                   max=70),
      numericInput(inputId = "UV",
                   label = "UV Index from 0 (indoors) to 10 (sunny, outside)",
                   value = 0,
                   min=0,
                   max=10),
      # Horizontal line ----
      tags$hr(),
      numericInput(inputId = "ventilation",
                   label = "Ventilation with outside air",
                   value = 0.7,
                   min=0),

      numericInput(inputId = "deposition",
                   label = "Deposition to surfaces",
                   value = 0.3,
                   min=0),
      numericInput(inputId = "control",
                   label = "Additional control measures",
                   value = 0,
                   min=0),
      numericInput(inputId = "n",
                   label = "Total of people present",
                   value = 5,
                   min=1),
      numericInput(inputId = "time2event",
                   label = "Number of days between Antigen Testing and Event",
                   value = 0,
                   min=0,
                   max=10),
      selectInput(inputId = "activity",
                   label="What activity will the participants be performing?",
                   choices = read_csv("quanta_emission_rates.csv")$Activity,
                   selected = NULL,
                   multiple=TRUE),
      radioButtons(inputId = "mask",
                   label="Will the participants be required to wear any mask",
                   choices = c("no" = 0,
                               "yes" = 1),
                   selected = 0,
                   inline = TRUE),
      numericInput(inputId = "prop_mask",
                   label="What proportion (%) of the participants do you expect to wear any mask?",
                   value=0,
                   min=0,
                   max=100),
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
      
      # Output: Histogram ----
      tabsetPanel(
        tabPanel("Disclaimer", htmlOutput("disclaimer")),
        tabPanel("Plot", htmlOutput("distRes"), plotOutput("distPlot")),
        #tabPanel("Table", tableOutput("probs")),
        tabPanel("Report", h3(htmlOutput("Report"))))
      
    )
  )
    
  
)


# Define server logic required to draw a histogram ----
server <- function(input, output, session) {

  dataInput <- reactive({
    
    
    #### Step 1: load the participants data + location data. In the absence of data, assume homogeneity.
    ####### Step 1.a: Load participants
    #df <- read.csv(input$file1$datapath,
    #               header = TRUE,
    #               sep = ",")
    DECAY = (7.56923714795655+1.41125518824508*(input$temperature-20.54)/10.66 +
            0.02175703466389*(input$RH-45.235)/28.665+7.55272292970083*((input$UV*0.185)-50) / 50 +
            (input$temperature-20.54)/10.66*(input$UV*0.185-50)/50*1.3973422174602)*60  #https://www.dhs.gov/science-and-technology/sars-airborne-calculator
    print(paste0("Decay:",  DECAY))
    df <- read.csv("mock_data.csv",
                                  header = TRUE,
                                    sep = ",")  ### for debugging
    ####### Enrich the dataset by computing the prevalence of the virus up
    ####### to 14 days before the event
<<<<<<< HEAD
    prevalence_df =  read_csv("chosen_prevalence_data.csv")#rep(0.005,  length(SENSITIVITY)) #extract_prevalence()
    filtered_prevalence_df = filter(prevalence_df, Date_of_infection<=as.Date(input$date_event) & Date_of_infection>=as.Date(input$date_event)-14)
    PREVALENCE = filtered_prevalence_df$infection_prevalence
    filtered_prevalence_df2 = filter(prevalence_df, Date_of_infection==as.Date(input$date_event))
    Baseline_infections_on_event_day<-rnorm(n=1000, mean=filtered_prevalence_df2$Infection_prevalence, sd=filtered_prevalence_df2$sd_Infection_prevalence)
    
    
=======
    prevalence_df =  read_csv("chosen_prevalence_data.csv") 
    filtered_prevalence_df = prevalence_df%>% filter(Date_of_infection<=input$date_event & Date_of_infection>=as.Date(input$date_event)-14)
    PREVALENCE = filtered_prevalence_df$Infection_prevalence

>>>>>>> 6eb4ea282e3414df0b58843231d7a05935f6ba9c
    ###### Step 1.b: compute room parameters for aerosolization
    volume = extract_volume(input$length, input$width, input$height, input$unit)
    first_order_loss_rate = extract_first_order(input$ventilation,
                                                input$control,
                                                DECAY,
                                                input$deposition)
    #print(paste0("first_order_loss_rate: ",first_order_loss_rate))
    ventilation_rate_per_person = extract_ventilation_rate_per_person(volume, 
                                                                      input$ventilation,
                                                                      input$control,
                                                                      input$n)
    #print(paste0("ventilation_rate_per_person : ",ventilation_rate_per_person ))
    ##########################################
    #### Step 2: compute probability of being infectious, hospitalized,and dying
    group_assignment = sapply(0:(length(SENSITIVITY)), function(x){paste0("p",x)})
    df[sapply(1:(length(SENSITIVITY)), function(x){paste0("p",x)})] =  t(sapply(1:nrow(df), function(x){
      compute_infectiousness_probability(sensitivity = SENSITIVITY,
                                         prevalence=PREVALENCE,
                                         df$profession[x],
                                         df$high_risk_contact[x],
                                         df$nb_people_hh[x]
                                         )
      }))
    df$p0 = 1 - apply(df[sapply(1:(length(SENSITIVITY)), function(x){paste0("p",x)})],1,sum)
    #print(df[group_assignment])
    df$p_hosp =  unlist(sapply(1:nrow(df), function(x){
      compute_hospitalization_probability(df$Age[x], df$Pregnant[x],
                                            df$Chronic_Renal_Insufficiency[x],
                                            df$Diabetes[x], df$Immunosuppression[x],
                                            df$COPD[x], df$Obesity[x], 
                                            df$Hypertension[x],
                                            df$Tobacco[x], df$Cardiovascular_Disease[x],
                                            df$Asthma[x], df$Sex[x])
    }))
    
    df$p_death =  unlist(sapply(1:nrow(df), function(x){
      compute_death_probability(df$Age[x], df$Pregnant[x],
                                            df$Chronic_Renal_Insufficiency[x],
                                            df$Diabetes[x], df$Immunosuppression[x],
                                            df$COPD[x], df$Obesity[x], 
                                            df$Hypertension[x],
                                            df$Tobacco[x], df$Cardiovascular_Disease[x],
                                            df$Asthma[x], df$Sex[x])
    }))
    
    

    ##########################################
    B = 1000
    nb_infective_people = rep(0,B)
    nb_infections = rep(0,B)
    p = rep(0,B)
    nb_hospitalizations = rep(0,B)
    nb_deaths = rep(0,B)
    #### Step 3: compute probability of infecting people and adverse outcomes using MCMC simulations
    for (b in 1:B){
      ####### draw infected people (their infectivity bucket)
      Z = apply(df[group_assignment], MARGIN = 1, function(x){which(rmultinom(1,1,x)>0)-1})
      Z = Z[(Z>0)]  #### keep infected people only
      nb_infective_people[b] = length(Z)
      if (nb_infective_people[b] > 0){
        ####### Step 3.a Direct contacts
        contacts  = rnorm(length(Z), mean = MU, sd = SD)
        nb_infections[b] = sapply(1:length(Z),
                                  function(x){rbinom(1,round(abs(contacts[x])), TAU *  RELATIVE_INFECTIOUSNESS[Z[x] + input$time2event])})
        
        ####### Step 3.b Aerosolization
        quanta_emission_rate <- compute_quanta_emission_rate(activity,
                                                             MASK_EFFICIENCY ,
                                                             input$prop_mask,
                                                             nb_infective_people[b])
        #print(paste0("quanta_emission_rate: ",quanta_emission_rate))

        quanta_concentration <- compute_quanta_concentation(quanta_emission_rate,
                                                           first_order_loss_rate,
                                                           volume,  
                                                           input$duration,
                                                           nb_infective_people[b])
        #print(paste0("quanta_concentration: ",quanta_concentration))
        quanta_inhaled_per_person <- compute_quanta_inhaled_per_person(quanta_concentration,
                                                                       BREATHING_RATE,
                                                                       input$duration,
                                                                       MASK_INHALATION_EFFICIENCY,
                                                                       0.01 * input$prop_mask)
        #print(paste0("quanta_inhaled_per_person: ",quanta_inhaled_per_person))
        
        p[b] = 1.-exp(-quanta_inhaled_per_person[[1]])
        nb_infections[b] =  nb_infections[b] + rbinom(1, input$n - nb_infections[b] - nb_infective_people[b],
                                                      p[b])
        
      }

      
    }
    
    #### Step 4: Compute (by MCMC simulation) the number of people with adverse outcomes
    for (b in 1:B){
      ###### Sample from the list
      #print(nb_infections[b])
      if (nb_infections[b] >0){
           nb_hospitalizations[b] = sum(sapply(sample(df$p_hosp, nb_infections[b]), function(x){rbinom(1,1,x)}))
           nb_deaths[b] = sum(sapply(sample(df$p_death, nb_infections[b]), function(x){rbinom(1,1,x)}))
      }
    }
    return(list(nb_infections=nb_infections, nb_deaths=nb_deaths, nb_hospitalizations=nb_hospitalizations))
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
  
  output$disclaimer = renderUI({
    tags$div(
      tags$p("This calculator uses the predicted number of infections in a country, the characteristics of the event, the participants and the screening protocol to estimate the number of infections, hospitalisations and deaths that are likely to result from holding the event"), 
      tags$p("Our method has been designed in collaboration with medical experts. However, we are not medical experts ourselves. Do not rely on this tool for medical advice."), 
      tags$p("We do not save any of the information you input or upload"),
      tags$p("To use this calculator, please fill out all of the questions on the left hand side, then click on the 'plot' or 'report' tabs at the top of the screen to see our predictions")
    )
  })
  

  
  pt1 <- reactive({
    ggplot(prevalence_df,aes(x=Date_of_infection, y=Infection_prevalence))+
    geom_point()+
    geom_line()+
    geom_errorbar(aes(ymin=Infection_prevalence-sd_Infection_prevalence,ymin=Infection_prevalence+sd_Infection_prevalence))+
    theme_classic()
  })
  
  pt2 <- reactive({
    ggplot(Baseline_infections_on_event_day)+
      geom_histogram()+
      theme_classic()
    
  })
  
  pt3 <- reactive({
    hist(x$nb_infections, col = "#75AADB", border = "white",
         xlab = "N",
         main ="Nb of predicted infections",
         na.rm=TRUE)
  })
  
  output$plotgraph = renderPlot({
    ptlist <- list(pt1(),pt2(),pt3())
    if (length(ptlist)==0) return(NULL)
    
    grid.arrange(grobs=ptlist,nrow=length(ptlist))
  })
  
  output$distRes <- renderText({
    x =  dataInput()
    #### Write the report
    conclusion = ""
    conclusion = paste0(conclusion, "</div>")
    HTML(conclusion)
  })
  
  
  
  
  
}


shinyApp(ui, server)
