# Define UI for app that draws a histogram ----
library("shiny")
library("tidyverse")
library(gridExtra)
library(lme4)
library(locfit)
library(data.table)
library(reshape2)
library(tibbletime)
source("proba_adverse_outcome.R")
source("individual_probabilities.R")
source("helper_functions.R")
source("vaccination.R")
source("aerosol_functions.R")
source("beta_params.R")
source("covid_case_predictions.R")
source("under_ascertainment_bias.R")
#setwd("~/Dropbox/aerosol_transmission_model/")


rolling_mean <- rollify(mean, window = 7)

BREATHING_RATE = 0.012 * 60 
DEPOSITION = 0.24
MASK_EFFICIENCY = 0.5  ### 50% is the recommended value
MASK_INHALATION_EFFICIENCY = 0.3
PRESSURE = 0.95
RELATIVE_INFECTIOUSNESS = rev(c(0, 0.01,0.05,0.2,0.6,0.88,0.98,1,1,1,0.95,0.8,0.4,0.2,0.1,0.01))
SENSITIVITY = rev(c(0,0,0.019,0.0327,0.560,0.653,0.718,0.746,0.737,0.718,0.7,0.68,0.662,0.644,0.625))
VACCINATION_EFFICACY_FIRST_DOSE = c(0,0, 0,rep(60, 25))/100 #### ASK JACK FOR THE DATA
VACCINATION_EFFICACY_SECOND_DOSE = c(60,70, 80,rep(95, 25))/100 #### ASK JACK FOR THE DATA
HOUSEHOLD_TRANSMISSION = 0.5
PERIOD_FOR_FITTING = 28
REGIONS = c("England", "Wales", "Scotland", "Northern Ireland")
TAU = 0.06/4
MU = 5
SD = 2
B = 1000

#### Includes regional vaccination data
VACCINATIONS <- read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations.csv")
VACCINATIONS_US <- read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/us_state_vaccinations.csv")
intersection_col = intersect(colnames(VACCINATIONS), colnames(VACCINATIONS_US))
VACCINATIONS <- rbind(VACCINATIONS[,intersection_col], VACCINATIONS_US[,intersection_col] %>% filter(location != "United States"))


f <- "cases_data"
if (!file.exists(f)){
  british_prev = read_csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&metric=newCasesByPublishDate&metric=newPeopleVaccinatedCompleteByPublishDate&metric=newPeopleVaccinatedFirstDoseByPublishDate&metric=newPeopleVaccinatedSecondDoseByPublishDate&format=csv")
  write_csv(data, file=f)
}else{
  british_prev = read_csv(f)
}

pop <- list("England" = 56286961, "Wales" = 3152879, "Scotland" = 5463300, "Northern Ireland" = 1893667)
british_prev["pop"] = sapply(british_prev$areaName, function(r){pop[[r]]})
british_prev  = british_prev %>% mutate(CasesPerMillion = newCasesByPublishDate/pop *1e6,
                                        FirstDosePerMillion  = newPeopleVaccinatedFirstDoseByPublishDate/pop *1e6,
                                        SecondDosePerMillion  = newPeopleVaccinatedCompleteByPublishDate/pop *1e6) %>%
  filter(date<Sys.Date())
colnames(british_prev)[4] = "location"
##### Extract regional prevalence rates
#english_prev <- readxl::read_xlsx("covid19infectionsurveydatasets20210319.xlsx", sheet="1f", trim_ws = TRUE, col_types = c("date",rep("numeric",81 )), skip=4)
# english_prev <- readxl::read_xlsx("covid19infectionsurveydatasets20210319.xlsx", sheet="1b", trim_ws = TRUE, col_types = c("date",rep("numeric",6 ),rep("text",3 )), skip=4)
# english_prev = english_prev[2:nrow(english_prev),]
# english_prev["Region"] = "England"
# irish_prev <- readxl::read_xlsx("covid19infectionsurveydatasets20210319.xlsx", sheet="4b", trim_ws = TRUE, col_types = c("date",rep("numeric",6 ),rep("text",3 )), skip=4)
# irish_prev = irish_prev[2:nrow(irish_prev),]
# irish_prev["Region"] = "Northern Ireland"
# wales_prev <- readxl::read_xlsx("covid19infectionsurveydatasets20210319.xlsx", sheet="3b", trim_ws = TRUE, col_types = c("date",rep("numeric",6 ),rep("text",3 )), skip=4)
# wales_prev = wales_prev[2:nrow(wales_prev),]
# wales_prev["Region"] = "Wales"
# scotland_prev <- readxl::read_xlsx("covid19infectionsurveydatasets20210319.xlsx", sheet="5b", trim_ws = TRUE, col_types = c("date",rep("numeric",6 ),rep("text",3 )), skip=4)
# scotland_prev =  scotland_prev[2:nrow(scotland_prev),]
# scotland_prev["Region"] = "Scotland"
# british_prev <- rbind(english_prev, irish_prev, scotland_prev, wales_prev)%>% separate("Ratio of estimated number of people testing positive for COVID-19", c("nb", "ratio"), sep = " in ")%>%
#   mutate("CasesPerMillion" = as.numeric(nb)/as.numeric(ratio) * 1e6)
# british_prev$date <- (as.Date(british_prev$date, "%Y-%m-%d"))


COUNTRY_DATA <- read.csv(file="https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv", header=T, sep=",")
COUNTRY_DATA$date <- (as.Date(COUNTRY_DATA$date, "%Y-%m-%d"))
bias<- compute_underascertainment_bias(Sys.Date() - (21+14), "United Kingdom", COUNTRY_DATA, Case_to_death_delay= 21, plot=FALSE)
bias_corr  = mean(bias$value)
bias_sd  = sd(bias$value)
##### COMPUTE UNDER-ASCERTAINMENT BIAS:

MU_INCUBATION = c(1.63,1.51, 1.75) 
SIGMA_INCUBATION = c( 0.50, 0.45, 0.55)
PROBA_ASYMPT = 0.24
PROBA_SYMPTOMS = as.numeric(table(factor(sapply(1:10000,function(b){mu = rnorm(1,1.63, 0.12); s = rnorm(1,0.5, 0.05); return(round(rlnorm(1,mu,s)))}), levels = 1:40))/(10000))* (1-PROBA_ASYMPT)

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
                   min=10,
                   max=30),
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
      numericInput(inputId = "prop_mask",
                   label="What proportion (%) of the participants do you expect to wear any mask?",
                   value=50,
                   min=0,
                   max=100),
      # Horizontal line ----
      tags$hr(),
      
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      tabsetPanel(
        tabPanel("Efficiency Screening Protocol", plotOutput("plotinput",  width = 800, height = 600)),
        tabPanel("Plot Results", plotOutput("plotgraph"), tableOutput("contents")),
      tabPanel("Disclaimer", htmlOutput("disclaimer"))
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
    
    ##### A BUNCH OF PRELIMINARY PARAMETERS
    MAX_DATE = as.Date(max(british_prev$date, na.rm = TRUE), fmt="%Y-%m-%d")
    #PERIOD_FOR_PREDICTING = as.numeric(as.Date(input$date_event)- as.Date(max(british_prev$date, na.rm = TRUE), fmt="%Y-%m-%d"))
    PERIOD_FOR_PREDICTING = as.numeric(as.Date(input$date_event)- MAX_DATE)
    N = list("England"=input$n_england,"Wales"=input$n_wales,"Scotland"=input$n_scotland,"Northern Ireland"=input$n_ni)
    N_TOT = sum(unlist(N))
    DECAY = max(0, (7.57+ 
                      1.41* (input$temperature-20.54)/10.66 +
                      0.0218 *(input$RH-45.24)/28.67 + 
                      7.55 *((input$UV*0.185)-50) / 50 +
                      (input$temperature-20.54)/10.66*(input$UV*0.185-50)/50*1.40) *60)  #https://www.dhs.gov/science-and-technology/sars-airborne-calculator
    
    
    participant_df = data.frame( "p_hosp" = 0.05, "p_death" = 0.01, "p_hh_transmission" = 0.3) #### to be made more tailored to the event at ahnd
    ##########################################
    ##########################################
    #### Step 1: COMPUTE NUMBER OF INFECTIOUS AND SUSCEPTIBLE PEOPLE AT THE EVENT
    ##########################################
    ##########################################
    
    
    ##########################################
    ####### 1a. Predict the prevalence using K-NN
    ##########################################
    D = 10
    withProgress(message = 'Computing distribution of the number of infected partipants', value = 0, {
      
    incProgress(1/D, detail = "Predicting the prevalence")
    a = 0
    for (region in REGIONS){
      data2fit = british_prev %>% filter(location == region, date > MAX_DATE - PERIOD_FOR_FITTING) %>% select(`CasesPerMillion`)
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
    #### Correct for all of these under-ascertainment issues
    future_prevalence_df = prevalence_df %>%
      dplyr::filter(time > 0) %>%
      mutate(prevalence = 1/bias_corr * prevalence,
             sd_prevalence = 1/bias_corr * sd_prevalence)
    
    ##########################################
    ####### 1b. Predict vaccination rate using loess fit.
    ##########################################
    incProgress(2/D, detail = "Predicting the vaccination rates")


    
    ##########################################
    #### Step 1.c : compute probability of being infectious, and vulnerable
    ##########################################
    group_assignment = c(sapply(1:PERIOD_FOR_PREDICTING,
                                function(x){paste0(x)}))
    df = reshape2::dcast(future_prevalence_df %>% select(region, time, prevalence), region ~ time, value.var = "prevalence")
    row.names(df) = df$region
    df = df[,group_assignment]
    nb_people_infected <- sapply(REGIONS, function(region){sum(N[region] * df[region,group_assignment ] )})
    df["p_inf"] = 1 - apply(df[,2:ncol(df)],1,sum)
    ##
    #### Effect of the test + symptoms screening
    SENSITIVITY = c(sapply(1: min(40, PERIOD_FOR_PREDICTING - input$time2event ), function(i){
      sum(PROBA_SYMPTOMS[1:(PERIOD_FOR_PREDICTING - input$time2event - i)]) * 0.2 * input$p_lie +
        0.4 * input$p_lie  * sum(PROBA_SYMPTOMS[(PERIOD_FOR_PREDICTING - input$time2event - i + 1): (PERIOD_FOR_PREDICTING - i)]) +
        (1-  sum(PROBA_SYMPTOMS[1: (PERIOD_FOR_PREDICTING - i)]) + PROBA_ASYMPT) * 0.4
    }),  1- sum(PROBA_SYMPTOMS[1:2]) * (1-input$p_lie) + PROBA_ASYMPT, 1- PROBA_SYMPTOMS[1] * (1-input$p_lie) + PROBA_ASYMPT, rep(1.0, max(0, PERIOD_FOR_PREDICTING -40) ))
    
    for (i in 1:PERIOD_FOR_PREDICTING){
      df[,  as.character(i)] = df[, as.character(i)] * SENSITIVITY[min(PERIOD_FOR_PREDICTING-i+1, length(SENSITIVITY))] * RELATIVE_INFECTIOUSNESS[min(PERIOD_FOR_PREDICTING-i+1, length(RELATIVE_INFECTIOUSNESS))]
    }
    
    df["p_comp"] = 1 - apply(df[,as.character(1:PERIOD_FOR_PREDICTING)],1,sum)
    nb_people_infectious_at_the_event <- sapply(REGIONS, function(region){sum(N[[region]] * df[region,group_assignment ] )})
    nb_people_detected <- nb_people_infected - nb_people_infectious_at_the_event 

    incProgress(4/D, detail = "Computing the number of Susceptible participants ")
    
    p_immunity = immunity("United Kingdom", input$date_event, VACCINATIONS, ODDS_RATIO_VACCINE_DOSE1, ODDS_RATIO_VACCINE_DOSE2, plot=FALSE)
    p_immunity_worst = immunity("United Kingdom", input$date_event, VACCINATIONS, ODDS_RATIO_VACCINE_DOSE1_WORST, ODDS_RATIO_VACCINE_DOSE2_WORST, plot=FALSE)
    p_immunity_best = immunity("United Kingdom", input$date_event, VACCINATIONS, ODDS_RATIO_VACCINE_DOSE1_BEST, ODDS_RATIO_VACCINE_DOSE2_BEST, plot=FALSE)
    beta_immunity_params = beta.parms.from.quantiles(c(p_immunity_worst,  p_immunity_best))
    nb_people_vulnerable_at_the_event <- sapply(REGIONS, function(r){mean(sapply(1:B, function(b){
      #### ASSUME FIRST DOSE BECOMES FULLY VACCINATED AFTER FOUR WEEKS
        sum(rbeta(N[[r]], beta_immunity_params$a, beta_immunity_params$b))
    }))})
    ######### COMPUTE ASSOCIATED DISTRIBUTIONS
    # dist_infectious<- bind_rows(lapply(1:B, function(b){
    #   t = data.frame(sapply(REGIONS, function(region){sapply(1:100, function(y){dpois(y, sum(N[[region]] * df[region,group_assignment ] ))})}))
    #   t["n"] = 1:100
    #   return(t)
    # }))
    ######### 
    
    
    

    
    
    ##########################################
    ##########################################
    ###### Step 2: TRANSMISSION DYNAMICS
    ##########################################
    ##########################################
    
    ##########################################
    ###### Step 2.a: compute room parameters for aerosolization
    ##########################################
    volume = extract_volume(input$length, input$width, input$height, input$unit)
    surface = extract_surface(input$length, input$width, input$unit)
    occupant_density = N_TOT /surface
    ventilation <- lookup_ventilation(input$ventilation, N_TOT, occupant_density, surface, volume)
    first_order_loss_rate = extract_first_order(ventilation,
                                                as.numeric(input$control),
                                                DECAY,
                                                DEPOSITION)
    
    ventilation_rate_per_person = extract_ventilation_rate_per_person(volume, 
                                                                      ventilation,
                                                                      as.numeric(input$control),
                                                                      N_TOT)
    
    quanta_emission_rate0 <- compute_quanta_emission_rate(input$activity,
                                                          MASK_EFFICIENCY ,
                                                          input$prop_mask/100,
                                                          1)
    
    print(paste0("ventilation_rate_per_person : ",ventilation_rate_per_person ))

    ##########################################
    #### Step 2.b: MC simulations
    ##########################################
    
    
    ##### Maybe I should add a negative binomial 
    nb_infective_people = rpois(B,input$n_england * sum(df["England", as.character(1:PERIOD_FOR_PREDICTING)]) + input$n_ni * sum(df["Northern Ireland", as.character(1:PERIOD_FOR_PREDICTING)]) + 
                                  input$n_scotland * sum(df["Scotland", as.character(1:PERIOD_FOR_PREDICTING)]) + input$n_wales * sum(df["Wales", as.character(1:PERIOD_FOR_PREDICTING)]))
    
    incProgress(5/D, detail = "Simulating the number of infectious participants")
    proba_dist_infectiousness <- sapply(1:min(100, N_TOT), function(n){dpois(n, input$n_england * sum(df["England", as.character(1:PERIOD_FOR_PREDICTING)]) + input$n_ni * sum(df["Northern Ireland", as.character(1:PERIOD_FOR_PREDICTING)]) + 
                                                                               input$n_scotland * sum(df["Scotland", as.character(1:PERIOD_FOR_PREDICTING)]) + input$n_wales * sum(df["Wales", as.character(1:PERIOD_FOR_PREDICTING)]))})
    
    
    
    ##########################################
    ##########################################
    ###### Step 3: RESULTS
    ##########################################
    ##########################################
    
    #### Step 3: compute probability of infecting people and adverse outcomes using MCMC simulations
    incProgress(6/D, detail = "Simulating the quanta concentrations")
    proba_infection <- sapply(1:min(100, N_TOT), function(n){
        
        q_e = n * as.numeric(quanta_emission_rate0)
        
        q_c <- as.numeric(compute_quanta_concentation(q_e, first_order_loss_rate = first_order_loss_rate, 
                                                      volume =volume,
                                                      duration =input$duration, 
                                                      nb_infective_people = n))
        
        quanta_inhaled_per_person <- compute_quanta_inhaled_per_person(quanta_concentration = q_c,  
                                                                       breathing_rate = BREATHING_RATE,
                                                                       duration=input$duration, 
                                                                       inhalation_mask_efficiency = MASK_INHALATION_EFFICIENCY,
                                                                       prop_mask = 0.01 * input$prop_mask)
        #print(c(q_e, q_c,quanta_inhaled_per_person))
        #print(paste0("quanta_inhaled_per_person: ",quanta_inhaled_per_person))
        ###### Now there is a lot of uncertainty around that parameter
        ###### We model potential super spreader events by using a pareto distirbution
        
        return(quanta_inhaled_per_person)
      })
    # proba_infection <-t(sapply(1:B, function(b){
    #     sapply(1:min(100, N_TOT), function(n){
    #       
    #       q_e = n * as.numeric(quanta_emission_rate0)
    #       
    #       q_c <- as.numeric(compute_quanta_concentation(q_e, first_order_loss_rate = first_order_loss_rate, 
    #                                          volume =volume,
    #                                          duration =input$duration, 
    #                                          nb_infective_people = n))
    #       
    #       quanta_inhaled_per_person <- compute_quanta_inhaled_per_person(quanta_concentration = q_c,  
    #                                                      breathing_rate = BREATHING_RATE,
    #                                                      duration=input$duration, 
    #                                                      inhalation_mask_efficiency = MASK_INHALATION_EFFICIENCY,
    #                                                      prop_mask = 0.01 * input$prop_mask)
    #       #print(c(q_e, q_c,quanta_inhaled_per_person))
    #       #print(paste0("quanta_inhaled_per_person: ",quanta_inhaled_per_person))
    #       ###### Now there is a lot of uncertainty around that parameter
    #       ###### We model potential super spreader events by using a pareto distirbution
    #       
    #       return(1.-exp(-rgamma(1, quanta_inhaled_per_person*100, rate = 100)))
    #     })
    #   }))
    
    n_susc <-  apply(sapply(1:N_TOT, function(b){rbinom(B, 1, 
                                              sum(rbeta(1, beta_immunity_params$a, beta_immunity_params$b)))}),1,sum)
    incProgress(7/D, detail = "Computing infections (this takes a while)")

    L  = min(max(which(proba_dist_infectiousness > 1e-10)), N_TOT)
    n_infections <- t(sapply(1:B, function(b){
      sapply(1:L, function(n){
        min(n_susc[b], rpareto(1,n_susc[b] * proba_infection[n] *0.16/1.16, 1.16))
      })
    }))
    res = data.frame("Average Number of Transmissions" = apply(as.matrix(n_infections),2,mean)%*% as.matrix(proba_dist_infectiousness[1:L]),
                     "Q97.5 Number of Transmissions" = apply(as.matrix(n_infections),2,quantile, 0.975)%*% as.matrix(proba_dist_infectiousness[1:L]),
                     "Q2.5 Number of Transmissions" = apply(as.matrix(n_infections),2,quantile, 0.025)%*% as.matrix(proba_dist_infectiousness[1:L]),
                     "Q99 Number of Transmissions" = apply(as.matrix(n_infections),2,quantile, 0.99)%*% as.matrix(proba_dist_infectiousness[1:L]),
                     "Q1 Number of Transmissions" = apply(as.matrix(n_infections),2,quantile, 0.01)%*% as.matrix(proba_dist_infectiousness[1:L]))

    
    incProgress(8/D, detail = "Computing summaries for the number of people infected")
    # N_SUSC = round(sum(unlist(nb_people_vulnerable_at_the_event)))
    # n_infections_summarized <- apply(as.matrix(n_infections),2,mean)%*% as.matrix(proba_dist_infectiousness[1:L])
    # n_infections_q975 <- apply(n_infections, 1, function(x){quantile(x, 0.975)}) %*% as.matrix(proba_dist_infectiousness)[1:L]
    # n_infections_q25 <- as.matrix(t(sapply(1:B, function(b){
    #   sapply(1:min(100, n_susc[b]), function(n){
    #     qbinom(0.025, rpois(1,n_susc[b]), proba_infection[b,n])
    #   })})) ) %*% as.matrix(proba_dist_infectiousness)
    # n_infections_conditional <- apply(as.matrix(n_infections)%*% as.matrix(proba_dist_infectiousness),1,sum)/sum(proba_dist_infectiousness)
    # n_infections_conditional_q975 <- as.matrix(t(sapply(1:B, function(b){
    #   sapply(1:min(100, N_SUSC), function(n){
    #     qbinom(0.975, rpois(1,N_SUSC), proba_infection[b,n])
    #   })})) ) %*% as.matrix(proba_dist_infectiousness)/sum(proba_dist_infectiousness)
    # n_infections_conditional_q25 <- as.matrix(t(sapply(1:B, function(b){
    #   sapply(1:min(100, N_SUSC), function(n){
    #     qbinom(0.025, rpois(1,N_SUSC), proba_infection[b,n])
    #   })})) ) %*% as.matrix(proba_dist_infectiousness)/sum(proba_dist_infectiousness)
    
    #n_infections_summarized = data.frame(infections = res$Average.Number.of.Transmissions ,
    #                                     p_infection = sum(proba_infection * proba_dist_infectiousness),
    #                                     p_conditional_infection = sum(proba_infection * proba_dist_infectiousness)/sum(proba_dist_infectiousness),
    #                                     q975_infections = res$Q97.5.Number.of.Transmissions, 
    #                                     q25_infections = res$Q02.5.Number.of.Transmissions, 
    #                                     sim=1:B) %>% mutate(type="Event")
    
    
        #### Step 4: Compute (by MCMC simulation) the number of people with adverse outcomes and secondary attack rates
        #for (b in which(is.nan(nb_infections) == FALSE)){
        ###### Sample from the list
    incProgress(8/D, detail = "Computing Baseline rates")   
    prev_day_event = filter(prevalence_df, time ==max(prevalence_df$time))
    proba_baseline <- t(sapply(1:B, function(b){
      a = sum(sapply(REGIONS, function(r){
        temp = prev_day_event %>% dplyr::filter(region==r)
        little_p = rnorm(n=1, mean=temp$prevalence,
                         sd=temp$sd_prevalence)
        max(min(little_p, 1), 0)
      }))
      return(a)
    }))
    
    infections_baseline <- data.frame(rbindlist(lapply(1:B, function(b){
      aa = data.frame(t(sapply(1:min(100, N_TOT), function(n){
        a  = proba_dist_infectiousness[n] 
        nn = min(rpois(1, sum(unlist(nb_people_vulnerable_at_the_event))), N_TOT - n ) ##round(min(sum(unlist(nb_people_vulnerable_at_the_event)), (N_TOT - n)))
        return(c(b, n, a, nn, proba_baseline[b], 
                 a * nn * proba_baseline[b],  qbinom(0.975,nn, a),qbinom(0.025,nn, a) ))
        
      })))
      colnames(aa) <- c("sim", "n", "p_infectious", "susceptible", "p_transmission", "av", "q975", "q25")
      return(aa)
    })))
    infections_baseline$sim = as.factor(infections_baseline$sim)
    n_infections_baseline = as_tibble(infections_baseline) %>% dplyr::group_by(sim) %>% dplyr::summarise(infections = sum(av),
                                                                                                         p_infection = sum(av/(susceptible)),
                                                                                                         p_conditional_infection = sum(av/(susceptible))/ sum(p_transmission),
                                                                                                         conditional_infections = sum(av)/ sum(p_transmission),
                                                                                                         q975_infections = sum(q975 *p_transmission ),
                                                                                                         q975_conditional_infections = sum(q975 *p_transmission)/sum(p_transmission),
                                                                                                         q25_infections = sum(q25 *p_transmission ),
                                                                                                         q25_conditional_infections = sum(q25 *p_transmission)/sum(p_transmission),
                                                                                                         ) %>% select(infections, conditional_infections, q975_infections, q975_conditional_infections,
                                                                                                                       q25_infections,q25_conditional_infections,
                                                                                                                      p_infection,
                                                                                                                      p_conditional_infection,sim) %>% mutate(type="Null")


    proba_damage <- t(sapply(1:B, function(b){
         ##### Need to introduce some heterogeneity there
         c(participant_df$p_hosp, participant_df$p_death, participant_df$p_hh_transmission)
    }))
    })
    
    

    
    
    return(list(res=res,
                proba_baseline = proba_baseline,
                proba_damage = proba_damage,
                prevalence_df=prevalence_df, n = input$n,
                n_infections_baseline  = n_infections_baseline ,
                participant_df = participant_df,
                samples_prev=samples_prev,
                nb_people_infected = nb_people_infected,
                nb_people_infectious_at_the_event =nb_people_infectious_at_the_event ,
                nb_people_detected = nb_people_detected,
                nb_people_vulnerable_at_the_event=nb_people_vulnerable_at_the_event,
                N=N))
  })
  
  output$contents <- renderTable({
    
    x =  dataInput()
    # set up cut-off values 
    breaks <- c(0,2,4,6,8,10,15,20,25,30,40,50,60,70,80,90,100,200, x$n)
    # specify interval/bin labels
    tags_x <- c("[0-2)","[2-4)", "[4-6)", "[6-8)", "[8-10)", "[10-15)", "[15-20)","[16-18)", "[18-20)",
                "[20-25)","[25-30)", "[30-40)", "[40-50)", "[50-60)", "[60-70)","[70-80)", "[80-90)","[90-100)","[100-)")
    # bucketing values into bins
    
    df <- data.frame("Setting"= c("Event", "Baseline", "Event, conditional on at least one infectious participant"))
    group_tags <- cut(unlist(x$n_infections_summarized%>% filter(type=="Event") %>% select(infections)), 
                      breaks=breaks, 
                      include.lowest=TRUE, 
                      right=FALSE)
    gp2 = cut(unlist(x$n_infections_summarized%>% filter(type=="Null") %>% select(infections)), 
              breaks=breaks, 
              include.lowest=TRUE, 
              right=FALSE)
    gp3 = cut(unlist(x$n_infections_summarized%>% filter(type=="Event") %>% select(conditional_infections)), 
              breaks=breaks, 
              include.lowest=TRUE, 
              right=FALSE)
    options(digits=7)
    df <- cbind(df, rbind(summary(group_tags)/(nrow(x$n_infections_summarized)/2), summary(gp2)/ (nrow(x$n_infections_summarized)/2),summary(gp3)/(nrow(x$n_infections_summarized)/2)))
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
    print(x)
    temp  =data.frame(x = c(0,50,100,150), 
                      y=c(20,20,20,20), 
                      n = sapply(REGIONS, function(r){x$N[[r]]}), 
                      text = sapply(REGIONS, function(r){paste0(r, "\n ", x$N[[r]])})
                      )
    
    pt1<- ggplot(temp, aes(x, y, fill = n, label = text)) +
      geom_tile(width = 20, height=20) + # make square tiles
      geom_text(color = "white", size=3) + # add white text in the middle
      coord_fixed() + # make sure tiles are square
      theme_void() +
      ggtitle("Number of participants per region")
    
    x$samples_prev$region = factor( x$samples_prev$region , levels=REGIONS)
    pt2 <- ggplot(x$samples_prev)+
      geom_line(aes(x=time, y=value, group = variable))+
      facet_wrap(.~region, nrow=1) + theme_bw()+
      ggtitle("Predicted Prevalence Per Region")+scale_y_log10()
      
    
    temp2  =data.frame(x = c(0,50,100,150), 
                       y=c(20,20,20,20),
                      n = c(x$nb_people_infected[["England"]], x$nb_people_infected[["Wales"]], x$nb_people_infected[["Scotland"]], x$nb_people_infected[["Northern Ireland"]]), 
                      text = sapply(REGIONS, function(r){paste0("Infected \n", r, " \n ",round(x$nb_people_infected[[r]],4) )}),
                      n_infectious = sapply(REGIONS, function(r){x$nb_people_infectious_at_the_event[[r]]}),
                      yi=c(0,0,0,0), xi = c(-10,40,90,140),
                      xv = c(10,60,110,160),
                      text_i = sapply(REGIONS, function(r){paste0("Infectious \n ", r, " \n ",round(as.numeric(x$nb_people_infectious_at_the_event[[r]]), 4))}),
                      text_v = sapply(REGIONS, function(r){paste0( "Susceptible \n", r, " \n ",round(x$nb_people_vulnerable_at_the_event[[r]],4) )}),
                      n_v =  sapply(REGIONS, function(r){x$nb_people_vulnerable_at_the_event[[r]]})
                      )
    
    
    pt3<- ggplot(temp2, aes(x, y, fill = n, label = text)) +
      geom_tile(width = 20, height = 20) + # make square tiles
      geom_text(color = "white", size=3) + # add white text in the middle
      coord_fixed() + # make sure tiles are square
      theme_void() +
      geom_tile(data = temp2, aes(xi, yi, fill = n_infectious), width = 20, height = 20) + # make square tiles+
      geom_text(data=temp2,aes(x=xi, y=yi,label=text_i), color="white", size=3) + 
      geom_tile(aes(xv, yi, fill = n_v, fill = n_v, label = text_v), width = 20, height = 20) + # make square tiles+
      geom_text(data=temp2,aes(x=xv, y=yi,label=text_v), color="white", size=3) + 
      ggtitle("Number of infectious and vulnerable participants per region")
    a=0
    for (r in REGIONS){
      if (a==0){
        temp3 = data.frame("Nb infectious"= rpois(1e6,x$nb_people_infectious_at_the_event[[r]]),
                           region=r)
      }else{
        temp3 = rbind(temp3,data.frame("Nb infectious"= rpois(1e6,x$nb_people_infectious_at_the_event[[r]]),
                           region=r))
      }
      a=a+1
    }
    pt4 <- ggplot(temp3)+
      geom_histogram(aes(x=`Nb.infectious`,y=..count../sum(..count..))) +
      facet_wrap(.~region, nrow=1) + theme_bw() + 
      scale_y_log10() + ylab("Log density") + xlab("Number infectious people at the event")
    
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
                      text_i = sapply(REGIONS, function(r)(paste0(r, "\n ", x$N[[r]])))
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
    pt1 <- ggplot(x$prevalence_df %>% filter(region=="England") %>% group_by(Date_of_cases) %>% summarise(prevalence=mean(prevalence), sd_prevalence=sum(sd_prevalence)), 
                  aes(x=Date_of_cases, y=prevalence*10^6))+
      geom_point()+
      geom_line() +
      geom_errorbar(aes(ymin= 10^6*(prevalence-2*sd_prevalence),
                        ymax=10^6*(prevalence+2*sd_prevalence)))+
      theme_classic()+
      scale_y_continuous(limits=c(0,NA), labels = scales::comma)+
      labs(title = "Projected Prevalence", x="Time (Days)", y="COVID cases per million")
    
    pt2 <-  ggplot(x$n_infections_summarized %>% filter(type=="Null"))+
      geom_density(aes(x=infections,y = (..count..)/sum(..count..), fill="Average number of infections"))+
      geom_density(aes(x=q975_infections,y = (..count..)/sum(..count..), fill="97.5th Quantile for the number of infections"), alpha=0.5)+
      geom_density(aes(x=q25_infections,y = (..count..)/sum(..count..), fill="2.5th Quantile for the number of infections"), alpha=0.5)+
      theme_classic() + xlim(0, 0.4)+
      labs(title="Distribution of the average number of infections \n (baseline, no event)",
           x ="Number of cases", y = "Probability")
    
    
    pt3 <-  ggplot(x$n_infections_summarized %>% filter(type=="Event"))+
      geom_density(aes(x=infections,y = (..count..)/sum(..count..), fill="Average number of infections"))+
      geom_density(aes(x=q975_infections,y = (..count..)/sum(..count..), fill="97.5th Quantile for the number of infections"), alpha=0.5)+
      geom_density(aes(x=q25_infections,y = (..count..)/sum(..count..), fill="2.5th Quantile for the number of infections"), alpha=0.5)+
      theme_classic() + xlim(0, 10) + 
      labs(title="Distribution of the average number of infections \n (With event)",
           x ="Number of cases", y = "Probability")
    
    pt4 <-  ggplot(x$n_infections_summarized %>% filter(type=="Event"))+
      geom_density(aes(x=conditional_infections,y = (..count..)/sum(..count..), fill="Average number of infections"))+
      geom_density(aes(x=q975_infections,y = (..count..)/sum(..count..), fill="97.5th Quantile for the number of infections"), alpha=0.5)+
      geom_density(aes(x=q25_infections,y = (..count..)/sum(..count..), fill="2.5th Quantile for the number of infections"), alpha=0.5)+
      theme_classic() + xlim(0, 25) + 
      labs(title="Distribution of the average number of infections, conditional on n_i>0 \n (With event)",
           x ="Number of cases", y = "Probability")
    
    # pt4 <- ggplot(data.frame(N=x$nb_infections + x$nb_secondary_infections))+
    #   geom_histogram(aes(x=N,y = (..count..)/sum(..count..)))+
    #   theme_classic() + 
    #   labs(title="Distribution of the number of total infections \n (primary and secondary, with event)",
    #        x ="Number of cases", y = "Probability")
    
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
