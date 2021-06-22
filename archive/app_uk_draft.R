# Define UI for app that draws a histogram ----
library("shiny")
library("tidyverse")
library(gridExtra)
library(lme4)
library(locfit)
source("proba_adverse_outcome.R")
source("individual_probabilities.R")
source("helper_functions.R")
source("aerosol_functions.R")
source("covid_case_predictions.R")
#setwd("~/Dropbox/aerosol_transmission_model/")

BREATHING_RATE = 0.012 * 60 
DEPOSITION = 0.24
MASK_EFFICIENCY = 0.5  ### 50% is the recommended value
MASK_INHALATION_EFFICIENCY = 0.3
PRESSURE = 0.95
RELATIVE_INFECTIOUSNESS = rev(c(0, 0.01,0.05,0.2,0.6,0.88,0.98,1,1,1,0.95,0.8,0.4,0.2,0.1,0.01))
SENSITIVITY = rev(c(0,0,0.019,0.0327,0.560,0.653,0.718,0.746,0.737,0.718,0.7,0.68,0.662,0.644,0.625))
VACCINATION_EFFICACY_FIRST_DOSE = c(0,0, 0,rep(60, 10))/100 #### ASK JACK FOR THE DATA
VACCINATION_EFFICACY_SECOND_DOSE = c(60,70, 80,rep(95, 10))/100 #### ASK JACK FOR THE DATA
HOUSEHOLD_TRANSMISSION = 0.5
PERIOD_FOR_FITTING = 28

TAU = 0.06/4
MU = 5
SD = 2


#### Includes regional vaccination data
VACCINATIONS <- read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations.csv")
VACCINATIONS_US <- read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/us_state_vaccinations.csv")
intersection_col = intersect(colnames(VACCINATIONS), colnames(VACCINATIONS_US))
VACCINATIONS <- rbind(VACCINATIONS[,intersection_col], VACCINATIONS_US[,intersection_col] %>% filter(location != "United States"))




##### Extract regional prevalence rates
#english_prev <- readxl::read_xlsx("covid19infectionsurveydatasets20210319.xlsx", sheet="1f", trim_ws = TRUE, col_types = c("date",rep("numeric",81 )), skip=4)
english_prev <- readxl::read_xlsx("covid19infectionsurveydatasets20210319.xlsx", sheet="1b", trim_ws = TRUE, col_types = c("date",rep("numeric",6 ),rep("text",3 )), skip=4)
english_prev = english_prev[2:nrow(english_prev),]
english_prev["Region"] = "England"
irish_prev <- readxl::read_xlsx("covid19infectionsurveydatasets20210319.xlsx", sheet="4b", trim_ws = TRUE, col_types = c("date",rep("numeric",6 ),rep("text",3 )), skip=4)
irish_prev = irish_prev[2:nrow(irish_prev),]
irish_prev["Region"] = "Northern Ireland"
wales_prev <- readxl::read_xlsx("covid19infectionsurveydatasets20210319.xlsx", sheet="3b", trim_ws = TRUE, col_types = c("date",rep("numeric",6 ),rep("text",3 )), skip=4)
wales_prev = wales_prev[2:nrow(wales_prev),]
wales_prev["Region"] = "Wales"
scotland_prev <- readxl::read_xlsx("covid19infectionsurveydatasets20210319.xlsx", sheet="5b", trim_ws = TRUE, col_types = c("date",rep("numeric",6 ),rep("text",3 )), skip=4)
scotland_prev =  scotland_prev[2:nrow(scotland_prev),]
scotland_prev["Region"] = "Scotland"
british_prev <- rbind(english_prev, irish_prev, scotland_prev, wales_prev)
COUNTRY_DATA <- read.csv(file="https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv", header=T, sep=",")
COUNTRY_DATA$date <- (as.Date(COUNTRY_DATA$date, "%Y-%m-%d"))

SENSITIVITY_FUNCTION <- read_csv("lancaster_pcr.csv")
PROBABILITY_SYMPTOMS <- c(0.0,0.0,1.5, 26.3, 45.0, 52.5, 57.8, 60.0,60.0, 60.0,60.0, 60.0, 60.0,60.0, 60.0)/100
PROBABILITY_INFECTIOUS <- c(0.0, 0.0, 1.4, 26, 67, 100, 100, 100 ,100, 100, 100, 86, 53, 31, 17, 8.42, 4.21, 2.11, 1.05)/100  #### AVERAGE OF TWO-DAY OFFSET
# regions_list =list ("England" =c("North East","North West","Yorkshire and The Humber",
#                               "East Midlands","West Midlands","East of England",
#                               "London","South East","South West"),
#                     "Wales" = c("Isle of Anglesey", "Gwynedd","Conwy", "Denbighshire",
#                                 "Flintshire","Wrexham",
#                                 "Ceredigion", "Pembrokeshire", "Carmarthenshire",
#                                 "Powys", "Caerphilly", "Blaenau Gwent", "Torfaen",
#                                 "Monmouthshire", "Newport",
#                                 "Swansea", "Neath Port Talbot",
#                                "Bridgend", "Rhondda Cynon Taf", "Merthyr Tydfil",
#                                 "Vale of Glamorgan", "Cardiff"),
#                     "Northern Ireland" = c("Northern", "Western","Belfast", "South Eastern",
#                                 "Southern"),
#                     "Scotland" = c("Na h-Eileanan Siar; Highland; Moray; Orkney Islands; Shetland Islands; Aberdeen City; Aberdeenshire; Argyll and Bute
#                                    Clackmannanshire; Falkirk; Stirling; Angus; Dundee City; Fife; Perth and Kinross
#                                    East Renfrewshire; Inverclyde; Renfrewshire; West Dunbartonshire; East Dunbartonshire; Glasgow City
#                                    East Lothian; Midlothian; City of Edinburgh; West Lothian
#                                    South Lanarkshire; North Lanarkshire
#                                    Dumfries and Galloway; East Ayrshire; North Ayrshire; Scottish Borders; South Ayrshire")
# )

ui <- fluidPage(
  
  # App title ----
  titlePanel("What are the parameters for your event?"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      numericInput(inputId = "n_england",
                   label = "Total number of people from England present at the event",
                   value = 1000,
                   min=0),
      numericInput(inputId = "n_wales",
                   label = "Total number of people from Wales present at the event",
                   value = 1000,
                   min=0),
      numericInput(inputId = "n_scotland",
                   label = "Total number of people from Scotland present at the event",
                   value = 1000,
                   min=0),
      numericInput(inputId = "n_ni",
                   label = "Total number of people from Northern Ireland present at the event",
                   value = 1000,
                   min=0),
      # selectInput(inputId = "region",
      #             label = "Which region do you live in?",
      #             choices = regions_list,
      #             selected="London"),
      # uiOutput("region"),
      dateInput(inputId = "date_event",
                label ="When will the event occur?", 
                value =as.character(Sys.Date()+1),
                min = "2020-03-01",
                max = as.character(Sys.Date() + 75),
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
      numericInput(inputId = "p_symptoms",
                   label = "Probability of developing symptoms ",
                   value = 0.6,
                   min=0,
                   max=1),
      numericInput(inputId = "p_lie",
                   label = "Probability of lying: pretending to feel fine even with symptoms",
                   value = 0.4,
                   min=0,
                   max=1),
      # Horizontal line ----
      tags$hr(),
      radioButtons(inputId = "unit",
                   label="What measurement unit will you be using?",
                   choices = c("ft" = 1,
                               "m" = 0),
                   selected = 0,
                   inline = TRUE),
      numericInput(inputId = "length",
                   label = "Length (in your selected metric). ",
                   value = 50,
                   min=0),
      numericInput(inputId = "width",
                   label = "Width (in your selected metric).",
                   value = 50,
                   min=0),
      numericInput(inputId = "height",
                   label = "Height (in your selected metric).",
                   value = 5,
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
                   label = "UV Index from 0 (indoors) to 10 (noon, sunny, outside)",
                   value = 0,
                   min=0,
                   max=10),
      # Horizontal line ----
      tags$hr(),
      selectInput(inputId = "ventilation",
                  label = "Which category will be the ventilation at the event be like?",
                  choices = read_csv("ventilation_parameters.csv")$category,
                  selected = "Daycare",
                  multiple=FALSE),
      radioButtons(inputId = "control",
                   label = "Additional control measures?",
                   choices = c("None"= 0,
                               "HEPA filter"= 440),
                   selected = 0), inline=TRUE,
      numericInput(inputId = "time2event",
                   label = "Number of days between Antigen Testing and Event",
                   value = 2,
                   min=0,
                   max=10),
      selectInput(inputId = "activity",
                  label="What activity will the participants be performing?",
                  choices = read_csv("quanta_emission_rates.csv")$Activity,
                  selected = "Standing:Loudly speaking",
                  multiple=TRUE),
      # radioButtons(inputId = "mask",
      #              label="Will the participants be required to wear any mask",
      #              choices = c("no" = 0,
      #                          "yes" = 1),
      #              selected = 0,
      #              inline = TRUE),
      numericInput(inputId = "prop_mask",
                   label="What proportion (%) of the participants do you expect to wear any mask?",
                   value=0,
                   min=0,
                   max=100),
      # radioButtons(inputId = "mixing",
      #              label="Are the people going to be mixing in that event?",
      #              choices = c("no" = 0,
      #                          "yes" = 1),
      #              selected = 0,
      #              inline = TRUE),
      
      # Input: Select a file ----
      # fileInput("file1", "Optional: Upload participants' information (choose CSV File)",
      #           multiple = FALSE,
      #           accept = c("text/csv",
      #                      "text/comma-separated-values,text/plain",
      #                      ".csv")),
      
      # Horizontal line ----
      tags$hr(),
      
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      tabsetPanel(
        tabPanel("Disclaimer", htmlOutput("disclaimer"), plotOutput("plotinput"))
        # tabPanel("Plot", htmlOutput("distRes"), plotOutput("plotgraph"), tableOutput("contents")),
        #tabPanel("Table", tableOutput("probs")),
        # tabPanel("Report", h3(htmlOutput("Report")), tableOutput("summary_participants_properties"),
        #       # tableOutput("summary_comparisons"))
      )
      
    )
  )
  
  
)



# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  
  dataInput <- reactive({
    
    
    #### Step 1: load the participants data + location data. In the absence of data, assume homogeneity.
    
    MAX_DATE = as.Date(max(british_prev$Date, na.rm = TRUE), fmt="%Y-%m-%d")
    PERIOD_FOR_PREDICTING = as.numeric(as.Date(input$date_event)- as.Date(max(british_prev$Date, na.rm = TRUE), fmt="%Y-%m-%d"))
    N = list("England"=input$n_england,"Wales"=input$n_wales,"Scotland"=input$n_scotland,"Northern Ireland"=input$n_ni)
    REGIONS = c("England", "Wales", "Scotland", "Northern Ireland")
    DECAY = max(0, (7.56923714795655 + 
                      1.41125518824508 * (input$temperature-20.54)/10.66 +
                      0.02175703466389*(input$RH-45.235)/28.665 + 
                      7.55272292970083*((input$UV*0.185)-50) / 50 +
                      (input$temperature-20.54)/10.66*(input$UV*0.185-50)/50*1.3973422174602) *60)  #https://www.dhs.gov/science-and-technology/sars-airborne-calculator
    
    ####### Enrich the dataset by computing the prevalence of the virus up
    ####### to 14 days before the event
    a = 0
    for (region in REGIONS){
      data2fit = british_prev %>% filter(Region == region, Date > MAX_DATE - PERIOD_FOR_FITTING) %>% select(`Modelled % testing positive for COVID-19`)
      res_df = compute_prevalence(as.Date(as.character(input$date_event)), data2fit,
                                  country_data =COUNTRY_DATA,nb_curves=100, 
                                  distance=Difference_function, 
                                  period4predicting=PERIOD_FOR_PREDICTING, period4fitting =28, 
                                  div=2, distance_type="MSE")
      if (a==0){
        prevalence_df = res_df$prediction_summary
        prevalence_df["region"] = region
        samples_prev = res_df$samples
        samples_prev["region"] = region
      }else{
        temp = res_df$prediction_summary
        temp["region"] = region
        temp_samples = res_df$samples
        temp_samples["region"] = region
        prevalence_df = rbind(prevalence_df, temp)
        samples_prev = rbind(samples_prev, temp_samples)
      }
      a = a + 1
    }
    
    ### Compute number of susceptible people
    #prevalence_df =  read_csv("prevalence_", input$country, "_data.csv")  #rep(0.005,  length(SENSITIVITY)) #extract_prevalence()
    future_prevalence_df = prevalence_df %>%
      dplyr::filter(time > 0)
    
    ###### Predict Vaccination Rates
    
    a = 0
    for (region in REGIONS){
      data2fit = VACCINATIONS %>% filter(location== region) %>% mutate(vaccinated = people_vaccinated_per_hundred - people_fully_vaccinated_per_hundred)
      ### Just run a simple logistic regression with regression splines (order 1)
      out0 = locfit(people_fully_vaccinated_per_hundred ~ lp(date, h=10, deg=1), data=data2fit %>%filter(is.na(vaccinated) == FALSE))
      #plot(unlist(data2fit%>%filter(is.na(vaccinated) == FALSE)%>% select(date)), unlist(data2fit%>%filter(is.na(vaccinated) == FALSE)%>% select(vaccinated)), 'l', col='blue',lwd=2)
      #lines(unlist(data2fit%>%filter(is.na(vaccinated) == FALSE)%>% select(date)), fitted(out0), 'l', col='red',lwd=2)
      ##temp_full_df = compute_vaccinations(as.character(input$date_event), data2fit,
      # country_data = VACCINATIONS,nb_curves=100, 
      # distance=Difference_function, 
      # period4predicting=NULL, period4fitting =NULL, 
      # div=2, distance_type="MSE")
      temp_full_df = data.frame(
        vaccinated = predict(out0, newdata = data.frame(date = seq(from=MAX_DATE+1, to= MAX_DATE +PERIOD_FOR_PREDICTING, by=1)))/100,
        time = 1:PERIOD_FOR_PREDICTING,
        "type"=rep("fully vaccinated", PERIOD_FOR_PREDICTING))
      temp_full_df["region"] = region
      model1 <- lm(vaccinated ~ date, data = data2fit %>% filter(date> as.Date("2021-02-01")))
      library(locfit)
      out1 = locfit(vaccinated ~ lp(date, h=10, deg=1), data=data2fit %>%filter(is.na(vaccinated) == FALSE))
      temp_df = data.frame(
        vaccinated = predict(out1, newdata = data.frame(date = seq(from=MAX_DATE+1, to= MAX_DATE +PERIOD_FOR_PREDICTING, by=1)))/100,
        time = 1:PERIOD_FOR_PREDICTING,
        "type"=rep("vaccinated", PERIOD_FOR_PREDICTING))
      temp_df["region"] = region
      
      
      ##temp_df = compute_vaccinations(as.character(input$date_event), data2fit, country="United Kingdom",
      # country_data = COUNTRY_DATA,nb_curves=100, 
      # distance=Difference_function, 
      # period4predicting=NULL, period4fitting =NULL, 
      # div=2, distance_type="MSE")
      temp_df["type"] = "vaccinated"
      if (a==0){
        vaccinations_df = rbind(temp_full_df, temp_df)
      }else{
        vaccinations_df = rbind(vaccinations_df, temp_full_df, temp_df)
      }
      a = a + 1
    }
    
    vaccinations_df = vaccinations_df %>% group_by(region, type) %>% 
      mutate(
        difference = vaccinated - lag(vaccinated)
      )
    vaccinations_df$difference = apply(vaccinations_df,1, function(x){as.numeric(ifelse(is.na(x["difference"]), x["vaccinated"], x["difference"]))})
    vac_df = dcast(vaccinations_df %>% select(region, time, difference, type), region +type  ~ time, value.var = "difference")
    ##########################################
    #### Step 2: compute probability of being infectious, and vulnerable
    group_assignment = c(sapply(1:PERIOD_FOR_PREDICTING,
                                function(x){paste0(x)}))
    df = dcast(filtered_prevalence_df %>% select(region, time, prevalence), region ~ time, value.var = "prevalence")
    row.names(df) = df$region
    df = df[,group_assignment]
    nb_people_infected <- sapply(REGIONS, function(region){sum(N[region] * df[region,group_assignment ] )})
    df["p_inf"] = 1 - apply(df[,2:ncol(df)],1,sum)
    df[,group_assignment] = df[, group_assignment] * ((1-PROBABILITY_SYMPTOMS) + input$p_symptoms * input$p_lie) ### these people won't know that they are infected
    ##
    #### Effect of the test
    for (i in 1:PERIOD_FOR_PREDICTING){
      df[,  as.character(i)] = df[, as.character(i)] * SENSITIVITY[min(PERIOD_FOR_PREDICTING-i+1, length(SENSITIVITY))] * RELATIVE_INFECTIOUSNESS[min(PERIOD_FOR_PREDICTING-i+1, length(RELATIVE_INFECTIOUSNESS))]
    }
    df["p_comp"] = 1 - apply(df[,as.character(1:PERIOD_FOR_PREDICTING)],1,sum)
    nb_people_infectious_at_the_event <- sapply(REGIONS, function(region){sum(N[region] * df[region,group_assignment ] )})
    nb_people_detected <- nb_people_infected - nb_people_infectious_at_the_event 
    
    nb_people_vulnerable_at_the_event <- sapply(REGIONS, function(r){
      
      
    })
    ######### COMPUTE ASSOCIATED DISTRIBUTIONS
    dist_infected<- sapply(1:B, function(b){
      
    })
    ######### 
    
    
    ###### Step 3: compute room parameters for aerosolization
    volume = extract_volume(input$length, input$width, input$height, input$unit)
    surface = extract_surface(input$length, input$width, input$unit)
    occupant_density = input$n /surface
    ventilation <- lookup_ventilation(input$ventilation, input$n, occupant_density, surface, volume)
    first_order_loss_rate = extract_first_order(ventilation,
                                                as.numeric(input$control),
                                                DECAY,
                                                DEPOSITION)
    
    ventilation_rate_per_person = extract_ventilation_rate_per_person(volume, 
                                                                      ventilation,
                                                                      as.numeric(input$control),
                                                                      input$n)
    
    quanta_emission_rate0 <- compute_quanta_emission_rate(input$activity,
                                                          MASK_EFFICIENCY ,
                                                          input$prop_mask,
                                                          1)
    
    quanta_concentration0 <- compute_quanta_concentation(quanta_emission_rate0,
                                                         first_order_loss_rate,
                                                         volume,  
                                                         input$duration,
                                                         1)
    print(paste0("ventilation_rate_per_person : ",ventilation_rate_per_person ))
    
    ##########################################
    B = 5000
    nb_infective_people = rep(0,B)
    nb_infections = rep(0,B)
    nb_secondary_infections = rep(0,B)
    p = rep(0,B)
    nb_hospitalizations = rep(0,B)
    nb_secondary_hospitalizations = rep(0,B)
    nb_deaths = rep(0,B)
    nb_secondary_deaths = rep(0,B)
    Baseline_infections_on_event_day <- rep(0, B)
    Baseline_hosp_on_event_day <- rep(0, B)
    Baseline_deaths_on_event_day <- rep(0, B)
    ##### Maybe I should add a negative binomoal
    nb_infective_people = rpois(B,input$n_england * sum(df["England", as.character(1:PERIOD_FOR_PREDICTING)]) + input$n_ni * sum(df["Northern Ireland", as.character(1:PERIOD_FOR_PREDICTING)]) + 
                                  input$n_scotland * sum(df["Scotland", as.character(1:PERIOD_FOR_PREDICTING)]) + input$n_wales * sum(df["Wales", as.character(1:PERIOD_FOR_PREDICTING)]))
    
    #### Step 3: compute probability of infecting people and adverse outcomes using MCMC simulations
    withProgress(message = paste0('Running ', B, ' simulations'), value = 0, {
      for (b in 1:B){
        ####### draw infected people (their infectivity bucket)
        # Z = c(rpois(input$n_england,df["England", group_assignment]),rmultinom(1,input$n_ni,df["Northern Ireland", group_assignment]),
        #       rmultinom(1,input$n_scotland,df["Scotland", group_assignment]),rmultinom(1,input$n_wales,df["Wales", group_assignment]))
        # nb_infective_people[b] = sum((Z>0))
        # Z = Z[(Z>0)]  #### keep infected people only
        if (nb_infective_people[b] > 0){
          ####### Step 3.a Direct contacts
          #contacts  =rnorm(length(Z), mean = MU, sd = SD)
          #contacts[contacts<1]=1
          #### Attempt at spatial modelling -- perhaps best to keep the two models separate for now as we are developping the tool
          #nb_infections[b] = sum(sapply(1:length(Z),
          #                              function(x){rbinom(1,round(abs(contacts[x])), 
          #                                                 TAU *  RELATIVE_INFECTIOUSNESS[min(16,Z[x] + input$time2event)])}))
          
          ####### Step 3.b Aerosolization
          #quanta_emission_rate <- quanta_emission_rate0 * nb_infective_people[b]
          #print(paste0("quanta_emission_rate: ",quanta_emission_rate))
          
          #quanta_concentration <- abs(quanta_concentration0) * sum(sapply(1:length(Z), function(x){RELATIVE_INFECTIOUSNESS[min(16,Z[x] + input$time2event)]}))
          #print(paste0("quanta_concentration: ",quanta_concentration))
          quanta_concentration <- abs(quanta_concentration0) * sum(sapply(1:length(Z), function(x){nb_infections}))
          quanta_inhaled_per_person <- compute_quanta_inhaled_per_person(quanta_concentration,
                                                                         BREATHING_RATE,
                                                                         input$duration,
                                                                         MASK_INHALATION_EFFICIENCY,
                                                                         0.01 * input$prop_mask)
          #print(paste0("quanta_inhaled_per_person: ",quanta_inhaled_per_person))
          
          p[b] = 1.-exp(-quanta_inhaled_per_person[[1]])
          #nb_infections[b] =  nb_infections[b] + rbinom(1, input$n - nb_infections[b] - nb_infective_people[b],
          #                                              p[b])
          nb_infections[b] =  rbinom(1, input$n - nb_infective_people[b],
                                     p[b])
          little_p = rnorm(n=1, mean=filtered_prevalence_df$Infection_prevalence,
                           sd=filtered_prevalence_df$sd_Infection_prevalence)
          Baseline_infections_on_event_day[b] <- nb_infective_people[b] + rbinom(1, input$n, max(min(little_p, 1), 0))
        }
        
        
        
        
        #### Step 4: Compute (by MCMC simulation) the number of people with adverse outcomes and secondary attack rates
        #for (b in which(is.nan(nb_infections) == FALSE)){
        ###### Sample from the list
        if (nb_infections[b] >0){
          nb_secondary_infections[b] = sum(sapply(sample(df$nb_people_hh, nb_infections[b]),function(x){
            ifelse(x>0, rbinom(1, x, HOUSEHOLD_TRANSMISSION),0)}))
          nb_hospitalizations[b] = sum(sapply(sample(df$p_hosp, nb_infections[b]), function(x){rbinom(1,1,x)}))
          nb_deaths[b] = sum(sapply(sample(df$p_death, nb_infections[b]), function(x){rbinom(1,1,x)}))
          if (nb_secondary_infections[b]>0){
            print("here")
            nb_secondary_hospitalizations[b] = sum(sapply(sample(df$p_hosp, nb_secondary_infections[b]), function(x){rbinom(1,1,x)}))
            nb_secondary_deaths[b] = sum(sapply(sample(df$p_death, nb_secondary_infections[b]), function(x){rbinom(1,1,x)}))
          }
        }
        if (Baseline_infections_on_event_day[b] >0){
          Baseline_hosp_on_event_day[b] = sum(sapply(sample(df$p_hosp, Baseline_infections_on_event_day[b]), function(x){rbinom(1,1,x)}))
          Baseline_deaths_on_event_day[b] = sum(sapply(sample(df$p_death, Baseline_infections_on_event_day[b]), function(x){rbinom(1,1,x)}))
        }
        incProgress(1/B, detail = paste("*"))
      }
    })
    
    
    return(list(nb_infections=nb_infections, nb_deaths=nb_deaths, nb_hospitalizations=nb_hospitalizations,
                Baseline_infections_on_event_day = Baseline_infections_on_event_day,
                Baseline_hosp_on_event_day = Baseline_hosp_on_event_day,
                Baseline_deaths_on_event_day = Baseline_deaths_on_event_day,
                prevalence_df=prevalence_df, n = input$n,
                nb_secondary_infections = nb_secondary_infections,
                nb_secondary_deaths=nb_secondary_deaths,
                nb_secondary_hospitalizations=nb_secondary_hospitalizations, 
                df=df))
  })
  
  output$contents <- renderTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    
    # or all rows if selected, will be shown.
    x =  dataInput()
    # set up cut-off values 
    breaks <- c(0,2,4,6,8,10,15,20,25,30,40,50,60,70,80,90,100,200, x$n)
    # specify interval/bin labels
    tags_x <- c("[0-2)","[2-4)", "[4-6)", "[6-8)", "[8-10)", "[10-15)", "[15-20)","[16-18)", "[18-20)",
                "[20-25)","[25-30)", "[30-40)", "[40-50)", "[50-60)", "[60-70)","[70-80)", "[80-90)","[90-100)")
    # bucketing values into bins
    df <- data.frame("Setting"= c("Event", "Baseline"))
    group_tags <- cut(x$nb_infections, 
                      breaks=breaks, 
                      include.lowest=TRUE, 
                      right=FALSE)
    gp2 = cut(x$Baseline_infections_on_event_day, 
              breaks=breaks, 
              include.lowest=TRUE, 
              right=FALSE)
    options(digits=7)
    df <- cbind(df, rbind(summary(group_tags)/length(x$nb_infections) ,summary(gp2)/length(x$nb_infections)))
    return(df)
    
    
  })
  
  output$disclaimer = renderUI({
    tags$div(
      tags$p("This calculator uses the predicted number of infections in a country, the characteristics of the event, the participants and the screening protocol to estimate the number of infections, hospitalisations and deaths that are likely to result from holding the event"), 
      tags$p("This risk estimator is not exact!!! It should not be understood as exact, but rather, to understand order of magnitude for the risk."), 
      tags$p("We do not save any of the information you input or upload"),
      tags$p("To use this calculator, please fill out all of the questions on the left hand side, then click on the 'plot' or 'report' tabs at the top of the screen to see our predictions")
    )
  })
  
  
  
  pt1 <- reactive({
    x =  dataInput()
    ggplot(x$prevalence_df,aes(x=Date_of_infection, y=Infection_prevalence))+
      geom_point() +
      geom_line() +
      geom_errorbar(aes(ymin=Infection_prevalence-sd_Infection_prevalence,ymax=Infection_prevalence+sd_Infection_prevalence))+
      theme_classic()
  })
  
  output$plotinput =renderPlot({
    
    x =  dataInput()
    temp  =data.frame(x = c(0,5,10,15), 
                      y=c(20,20,20,20), 
                      n = c(x$n_england, x$n_wales, x$n_scotland, x$n_ni), 
                      text = c(paste0("England \n ", x$n_england),  paste0("Wales \n ", x$n_wales), 
                               paste0("Scotland \n", x$n_scotland),paste0("Northern Ireland \n ", x$n_ni))
    )
    pt1<- ggplot(temp, aes(x, y, fill = n, label = text)) +
      geom_tile(width = 4, height = 4) + # make square tiles
      geom_text(color = "white") + # add white text in the middle
      coord_fixed() + # make sure tiles are square
      theme_void() +
      ggtitle("Number of participants per region")
    
    pt2<- ggplot(temp, aes(x, y, fill = n, label = text)) +
      geom_tile(width = 4, height = 4) + # make square tiles
      geom_text(color = "white") + # add white text in the middle
      coord_fixed() + # make sure tiles are square
      theme_void() +
      ggtitle("Number of infectious and vulnerable participants per region")
    
    grid.arrange(pt1,pt2,pt3,pt4,nrow=4)
    
  })
  
  output$plotinfected =renderPlot({
    
    x =  dataInput()
    temp  =data.frame(x = c(0,5,10,15), 
                      y=c(20,20,20,20), 
                      n = x$n_infected, 
                      text = c(paste0("England \n ", x$n_england),  paste0("Wales \n ", x$n_wales), 
                               paste0("Scotland \n", x$n_scotland),paste0("Northern Ireland \n ", x$n_ni))
    )
    ggplot(temp, aes(x, y, fill = n, label = text)) +
      geom_tile(width = 4, height = 4) + # make square tiles
      geom_text(color = "white") + # add white text in the middle
      coord_fixed() + # make sure tiles are square
      theme_void() +
      ggtitle("Average Number of infected participants per region")
    
  })
  
  output$plotinfectious=renderPlot({
    
    x =  dataInput()
    temp  =data.frame(x = c(0,5,10,15), 
                      y=c(20,20,20,20), 
                      x_v=c(2.5,7.5,12.5,17.5), 
                      n = x$n_infectious, 
                      v = x$n_vulnerable,
                      text_i = c(paste0("England \n ", x$n_england),  paste0("Wales \n ", x$n_wales), 
                                 paste0("Scotland \n", x$n_scotland),paste0("Northern Ireland \n ", x$n_ni))
    )
    ggplot(temp, aes(x, y, fill = n, label = text)) +
      geom_tile(width = 4, height = 4) + # make square tiles
      geom_text(color = "white") + # add white text in the middle
      coord_fixed() + # make sure tiles are square
      theme_void() +
      ggtitle("Average Number of infectious and susceptible participants per region")
    
  })
  
  output$plotgraph = renderPlot({
    x =  dataInput()
    print(x$prevalence_df)
    pt1 <- ggplot(x$prevalence_df, aes(x=Date_of_infection, y=Infection_prevalence*10^6))+
      geom_point()+
      geom_line()+
      geom_errorbar(aes(ymin=10^6*(Infection_prevalence-sd_Infection_prevalence),
                        ymax=10^6*(Infection_prevalence+sd_Infection_prevalence)))+
      theme_classic()+
      scale_y_continuous(limits=c(0,NA), labels = scales::comma)+
      labs(title = "Projected Prevalence", x="Time (Days)", y="COVID cases per million")
    
    pt2 <-  ggplot(data.frame(N=x$Baseline_infections_on_event_day))+
      geom_histogram(aes(N,y = (..count..)/sum(..count..)))+
      theme_classic() +
      labs(title="Distribution of the number of infections \n (baseline, no event)",
           x ="Number of cases", y = "Probability")
    
    
    pt3 <- ggplot(data.frame(N=x$nb_infections))+
      geom_histogram(aes(x=N,y = (..count..)/sum(..count..)))+
      theme_classic() + 
      labs(title="Distribution of the number of (primary) infections \n (with event)",
           x ="Number of cases", y = "Probability")
    
    pt4 <- ggplot(data.frame(N=x$nb_infections + x$nb_secondary_infections))+
      geom_histogram(aes(x=N,y = (..count..)/sum(..count..)))+
      theme_classic() + 
      labs(title="Distribution of the number of total infections \n (primary and secondary, with event)",
           x ="Number of cases", y = "Probability")
    
    ptlist <- list(pt1,pt2,pt3, pt4)
    if (length(ptlist)==0) return(NULL)
    
    grid.arrange(pt1,pt2,pt3,pt4,nrow=2)
  })
  
  output$distRes <- renderText({
    x =  dataInput()
    #### Write the report
    conclusion = ""
    conclusion = paste0(conclusion, "</div>")
    HTML(conclusion)
  })
  
  
  output$summary_participants_properties<-renderTable({
    x =  dataInput()
    #### Renders table with summary statistics for the participants (age and gender, nb of people in household)
    df = data.frame(rbind(
      c(paste0(nrow(input$n), ' participants'), paste0(100 *mean(unlist(x$df["Sex"]) == "Male"), " % Male"), 
        paste0(100 *mean(unlist(x$df["Sex"]) == "Female"), " % Female")),
      c(mean(unlist(x$df["Age"])), quantile(unlist(x$df["Age"]), 0.025), quantile(unlist(x$df["Age"]),0.975)),
      c(mean(unlist(x$df["nb_people_hh"])), quantile(unlist(x$df["nb_people_hh"]), 0.025),
        quantile(unlist(x$df["nb_people_hh"]), 0.975))
    ))
    colnames(df) <- c("Mean", "2.5th Quantile", "97.5th Quantile")
    row.names(df) <- c("General", "Age", "Nb of people in household")
    print(df)
  })
  
  output$summary_comparisons <-renderTable({
    #### Renders table with summary statistics for the participants (age and gender, nb of people in household)
    x = dataInput()
    df = data.frame(rbind(
      c(mean(x$nb_infections), quantile(x$nb_infections, 0.025), quantile(x$nb_infections, 0.975)),
      c(mean(x$Baseline_infections_on_event_day), quantile(x$Baseline_infections_on_event_day, 0.025), quantile(x$Baseline_infections_on_event_day, 0.975)),
      c(mean(x$nb_hospitalizations), quantile(x$nb_hospitalizations, 0.025), quantile(x$nb_hospitalizations, 0.975)),
      c(mean(x$Baseline_hosp_on_event_day), quantile(x$Baseline_hosp_on_event_day, 0.025), quantile(x$Baseline_hosp_on_event_day, 0.975)),
      c(mean(x$nb_deaths), quantile(x$nb_deaths, 0.025), quantile(x$nb_deaths, 0.975)),
      c(mean(x$Baseline_deaths_on_event_day), quantile(x$Baseline_deaths_on_event_day, 0.025), quantile(x$Baseline_deaths_on_event_day, 0.975))
    ))
    colnames(df) <- c("Mean", "2.5th Quantile", "97.5th Quantile")
    row.names(df) <- c("Event: Infections", "Baseline: Infections",
                       "Event: Hospitalizations", "Baseline: Hospitalizations", 
                       "Event: Deaths", "Baseline: Deaths")
    print(df)
  })
  
  
  comparison_with_other_diseases<-reactive({
    #### Renders comparison of the fatality of COVID with respect to other diseases
    "By "
  })
  
  comparison_with_car_accidents <-reactive({
    #### Renders comparison of fatalities due to car accidents
  })
  
  comparison_with_H0 <-reactive({
    #### Renders comparison with H0 (table)
    
  })
  
  comparison_with_H0_text<-reactive({
    dat =  dataInput()
    
    
    #### Renders comparison with H0 (text ---- this is an important comparison, so it perhaps deserves more explanations than the rest).
    i = round(mean(dat$nb_infections),2)
    ii =round( 0.5 * (quantile(x = dat$nb_infections, 0.975) - quantile(x = dat$nb_infections, 0.025)),2)
    j = round(mean(dat$nb_hospitalizations),2)
    jj =round( 0.5 * (quantile(x = dat$nb_hospitalizations, 0.975) - quantile(x = dat$nb_hospitalizations, 0.025)),2)
    k = round(mean(dat$nb_deaths),2)
    kk =round( 0.5 * (quantile(x = dat$nb_deaths, 0.975) - quantile(x = dat$nb_deaths, 0.025)),2)
    
    x = round(mean(dat$Baseline_infections_on_event_day),2)
    xx = round(0.5 * (quantile(x = dat$Baseline_infections_on_event_day, 0.975) - quantile(x = dat$Baseline_infections_on_event_day, 0.025)),2)
    y = round(mean(dat$Baseline_hosp_on_event_day),2)
    yy =round( 0.5 * (quantile(x = dat$Baseline_hosp_on_event_day, 0.975) - quantile(x = dat$Baseline_hosp_on_event_day, 0.025)),2)
    z = round(mean(dat$Baseline_deaths_on_event_day),2)
    zz = round(0.5 * (quantile(x = dat$Baseline_deaths_on_event_day, 0.975) - quantile(x = dat$Baseline_deaths_on_event_day, 0.025)),2)
    a = round(i/x,4)
    b=  round(j/y,4)
    d = round(ii/ xx,4)
    p = paste0("The event is projected to yield a total of ", i, " new infections (+/-", ii , "), ", j, " hospitalizations (+/-", jj , "), and ", k, " deaths (+/-", kk , ").  \n")
    #print(p)
    if (input$date_event > Sys.Date()){
      p = paste0(p, "By comparison, if the event were not to take place, based on the projected prevalence of the disease, we could expect ", x, "(+/-", xx , ") new infections among the participants, yielding ",
                 y, "(+/-", yy , ") hospitalizations and ", z, "(+/-", zz , ") deaths \n. As such, the average effect of the event would be an increase of ",
                 ifelse( a>1, "an increase of ", "a decrease of "),
                 100 * abs(a-1), " % new infections among participants (", ifelse( d>1, "an increase of ", "a decrease of "),
                 100 * abs(d-1), " % in the 95th quantile), and ", ifelse( b>1, "an increase of ", "a decrease of "), abs(b -1)*100,
                 " % new hospitalizations (primary). Check the table below for a more complete comparison of primary (and secondary) infections, hosptializations and deaths.")
    }else{
      p = paste0(p, "By comparison, if the event had not taken place, based on the observed prevalence of the disease, we could expect ", x, "(+/-", xx , ") new infections among the participants, yielding ",
                 y, "(+/-", yy , ") hospitalizations and ", z, "(+/-", zz , ") deaths. \n As such, the average effect of the event is ",
                 ifelse(a>1, "an increase of ", "a decrease of "),
                 abs(a-1)*100, " new infections among participants (", ifelse(d>1, "an increase of ", "a decrease of "),
                 abs(d-1)*100, " in the 95th quantile), and ", ifelse(b>1, "an increase of ", "a decrease of "), abs(b -1)*100,
                 " new hospitalizations (primary). Check the table below for a more complete comparison of primary (and secondary) infections, hosptializations and deaths.")
    }
    #print(p)
  })
  
  
  output$Report <- renderText({
    ### Report that needs to quantify the risk of holding the event with respect to other known diseases/ prevelence
    ### We need: (1) the H0 (better visaulization)
    ###          (2) Comparison with the risk of having car accidents
    ###          (3) Other diseases to compare the fatality with
    gsub(pattern = "\\n", replacement = "<br/>" ,paste0('<p style="font-family:verdana;font-size:14px">',
                                                        toString(comparison_with_H0_text()),
                                                        toString(comparison_with_other_diseases()),  
                                                        "</p>"))
    
  })
  
  
  
  
  
  
}


shinyApp(ui, server)
