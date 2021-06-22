# Define UI for app that draws a histogram ----
library("shiny")
library("tidyverse")
library(gridExtra)
library(lme4)
library(locfit)
library(data.table)
library(reshape2)
library(tibbletime)
library(EnvStats)


source("helper_functions.R")
source("vaccination.R")
source("aerosol_functions.R")
source("beta_params.R")
source("covid_case_predictions.R")
source("under_ascertainment_bias.R")
source("screening_efficiency.R")
source("relative_infectiousness.R")
#setwd("~/Dropbox/aerosol_transmission_model/")


url <- "https://twitter.com/intent/tweet?text=Event%20Risk%20Estimation%20Tool&url=https://homecovidtests.shinyapps.io/CERTIFIC_risk/"
url_fb <- "http://www.facebook.com/share.php?u=https://homecovidtests.shinyapps.io/CERTIFIC_risk/&quote=Event%20Risk%20Estimation%20Tool"
share <- list(
  title = "Event Risk Estimation Tool for the Management of Live Events",
  url = "http://www.facebook.com/share.php?u=https://homecovidtests.shinyapps.io/CERTIFIC_risk/",
  description = "Given the parameters for an event, display the associated transmission risk. For reference and information purpose. Not for medical use."
)

###### START BY DEFINING GLOBAL PARAMETERS
B = 1000  ### Nb of simulations we want to run
BREATHING_RATE = 0.012 * 60    
DEPOSITION = 0.24
MASK_EFFICIENCY = 0.8  ### 50% is the recommended value
MASK_INHALATION_EFFICIENCY = 0.5
PRESSURE = 0.95
#### GLOBAL FUNCTIONS
rolling_mean <- rollify(mean, window = 7)

#### CASES AND VACCINATION DATA
COUNTRY_DATA <- read.csv(file="https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv", header=T, sep=",")
COUNTRY_DATA$date <- (as.Date(COUNTRY_DATA$date, "%Y-%m-%d"))
ifr_data=read.csv("country_ifr_data.csv", header=T)

COUNTRY_LIST =  intersect(unique(COUNTRY_DATA$location),ifr_data$location)


##### COMPUTE UNDER-ASCERTAINMENT BIAS



ui <- fluidPage(
  
  # App title ----
  titlePanel("What are the parameters for the event?"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      selectInput(inputId = "country",
                  label = "Which country will the event be in",
                  choices = COUNTRY_LIST,
                  selected = "United Kingdom"),
      numericInput(inputId = "n_people",
                   label = "Total number of people present at the event",
                   value = 1000,
                   min=0),
      dateInput(inputId = "date_event",
                label ="When will the event occur?", 
                value =as.character(Sys.Date()+14),
                min = "2020-03-01",
                max = as.character(Sys.Date() + 60),
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
      numericInput(inputId = "p_lie",
                   label = "Probability of lying: pretending to feel fine even with symptoms",
                   value = 0.5,
                   min=0,
                   max=1),
      numericInput(inputId = "prop_mask",
                   label="What proportion (%) of the participants do you expect to wear any mask?",
                   value=100,
                   min=0,
                   max=100),
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
                  label = "Which category will the ventilation at the event be like?",
                  choices = read_csv("ventilation_parameters.csv")$category,
                  selected = "Music/Theater/Dance",
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
      # Horizontal line ----
      # Horizontal line ----
      tags$hr(),
      # Create url with the 'twitter-share-button' class
      tags$a(href=url, "Tweet", class="twitter-share-button"),
      # Copy the script from https://dev.twitter.com/web/javascript/loading into your app
      # You can source it from the URL below. It must go after you've created your link
      includeScript("http://platform.twitter.com/widgets.js"),
      #submitButton("Calculate", icon("refresh")),
      tags$a(href=url_fb, "Share on Facebook", class="fb-share-button")
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      # Output: Histogram ----
      tabsetPanel(
        tabPanel("Disclaimer", htmlOutput("disclaimer")),
        tabPanel("Results", plotOutput("plotinput",  width = 800, height = 1000),tableOutput("contents"))
       )
      )
  )
  
)



# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  dataInput <- reactive({
    
    
    bias<- compute_underascertainment_bias(min(as.Date(as.character(input$date_event)),  Sys.Date()) - (21+14),
                                           input$country, COUNTRY_DATA, 
                                           Case_to_death_delay= 21, 
                                           date_max=min(as.Date(as.character(input$date_event)) +14, Sys.Date()),
                                           plot=FALSE)
    
    bias_corr  = mean(bias$value)
    bias_corr = ifelse(is.na(bias_corr), 1, bias_corr )
    bias_sd  = sd(bias$value)
    
    

    ##### A BUNCH OF PRELIMINARY PARAMETERS
    MAX_DATE = as.Date(Sys.Date() -1, fmt="%Y-%m-%d")
    #PERIOD_FOR_PREDICTING = as.numeric(as.Date(input$date_event)- as.Date(max(british_prev$date, na.rm = TRUE), fmt="%Y-%m-%d"))
    PERIOD_FOR_PREDICTING = max(28, (as.numeric(as.Date(input$date_event)- MAX_DATE-1)))
    PERIOD_FOR_FITTING = 28
    infectiousness_all = compute_relative_infectiousness(input, period4predicting = PERIOD_FOR_PREDICTING, plot=FALSE) 


    N_TOT = input$n_people
    DECAY = max(0, (7.57+ 
                      1.41* (input$temperature-20.54)/10.66 +
                      0.0218 *(input$RH-45.24)/28.67 + 
                      7.55 *((input$UV*0.185)-50) / 50 +
                      (input$temperature-20.54)/10.66*(input$UV*0.185-50)/50*1.40) *60)  #https://www.dhs.gov/science-and-technology/sars-airborne-calculator
    
  
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
      
      if(as.Date(input$date_event) >=Sys.Date()+ 2){
          data2fit = COUNTRY_DATA %>% filter(date >= Sys.Date() - PERIOD_FOR_FITTING , location == input$country) %>% select(new_cases_smoothed_per_million)
          if (nrow(data2fit) == 0){
            data2fit = COUNTRY_DATA %>% filter(date >= Sys.Date() - PERIOD_FOR_FITTING , location == "World") %>% select(new_cases_smoothed_per_million)
          }
          res_df = compute_prevalence(as.Date(as.character(input$date_event)), data2fit,
                                        country_data =COUNTRY_DATA, nb_curves=100, 
                                        distance=Difference_function, 
                                        period4predicting=PERIOD_FOR_PREDICTING + 1, period4fitting = min(PERIOD_FOR_FITTING, nrow(data2fit)), 
                                        div=2, distance_type="MSE")
          prevalence_df = res_df$res
          samples_prev = res_df$output
    
          #### Correct for all of these under-ascertainment issues
          future_prevalence_df = prevalence_df %>%
            dplyr::filter(time > 0) %>%
            mutate(prevalence = 1/bias_corr * prevalence,
                   sd_prevalence = 1/bias_corr * sd_prevalence)
      }else{
        data2fit = COUNTRY_DATA %>% filter(date >= as.Date(as.character(input$date_event)) - PERIOD_FOR_PREDICTING,
                                           date <= as.Date(as.character(input$date_event)), location == input$country) %>%  select(new_cases_smoothed_per_million)
        #### Extrapolate if the data is 
        if (nrow(data2fit) < PERIOD_FOR_PREDICTING +1){
          data2fit= rbind(data2fit, data.frame("new_cases_smoothed_per_million" = rep(data2fit$new_cases_smoothed_per_million[nrow(data2fit)], PERIOD_FOR_PREDICTING +1 -nrow(data2fit) ) ))
        }
        future_prevalence_df  = data.frame("time" = 1:(PERIOD_FOR_PREDICTING + 1),
                                           "prevalence" =  1/bias_corr * 1e-6 * data2fit$new_cases_smoothed_per_million,
                                           "sd_prevalence" = rep(0,PERIOD_FOR_PREDICTING+1),
                                           "Date_of_cases" =seq(from=as.Date(as.character(input$date_event)) - PERIOD_FOR_PREDICTING, to = as.Date(as.character(input$date_event)), by="day"))
      }
      
      future_prevalence_df = future_prevalence_df %>% mutate(ymin = prevalence -2*sd_prevalence,
                                                   ymax = prevalence +2*sd_prevalence)
      future_prevalence_df["ymin"] = sapply(future_prevalence_df["ymin"], function(x){ifelse(x<0,0,x)})
      future_prevalence_df["ymax"] = sapply(future_prevalence_df["ymax"], function(x){ifelse(x>1,1,x)})
      
      

      
      
      ##########################################
      #### Step 1.b : compute probability of being infectious, and vulnerable
      ##########################################
      incProgress(2/D, detail = "Predicting the infectiousness (This might take a while).")
      group_assignment = c(sapply(1:PERIOD_FOR_PREDICTING, function(x){paste0(x)})) #### STOP one day before the event
      
      df = pivot_wider(future_prevalence_df %>% select(time, prevalence), names_from = c("time"), values_from = "prevalence")
      #proba_null <- future_prevalence_df[PERIOD_FOR_PREDICTING,"prevalence"]
      nb_people_infected <- sum(N_TOT * df)
      
      #### Effect of the test + symptoms screening
      df_sample = data.frame(matrix(0, B,PERIOD_FOR_PREDICTING))
      colnames(df_sample) = group_assignment
      
      for (i in 1:(PERIOD_FOR_PREDICTING)){
        ii = PERIOD_FOR_PREDICTING - i
        ind = which(infectiousness_all$Date.of.Infection == -(PERIOD_FOR_PREDICTING - i))
        
        df_sample[as.character(i)]= sapply(1:B, function(b){
          min(1,max(0,ifelse(future_prevalence_df$sd_prevalence[which(future_prevalence_df$time == i)] >0, rnorm(1, future_prevalence_df$prevalence[which(future_prevalence_df$time == i)], future_prevalence_df$sd_prevalence[which(future_prevalence_df$time == i)]),
                             future_prevalence_df$prevalence[which(future_prevalence_df$time == i)] )))* max(min(1, 0.01* ifelse( (infectiousness_all$infectiousness_event_sd[ind] > 0),
                                                                                                                                rnorm(1, infectiousness_all$infectiousness_event[ind], infectiousness_all$infectiousness_event_sd[ind]),
                                                                                                                                infectiousness_all$infectiousness_event[ind])),0)
        
                    })
      }
      
      
      
      nb_people_infectious_at_the_event <- mean(N_TOT* apply(df_sample[group_assignment ],1, sum))
      nb_infective_people <- sapply(1:B, function(b){rbinom(1,N_TOT, sum(df_sample[b,group_assignment ]))})
      nb_people_detected <- nb_people_infected - nb_people_infectious_at_the_event 

      
      incProgress(3/D, detail = "Computing the number of Susceptible participants ")
      
      if (as.Date(input$date_event) > as.Date("2021-02-15")){
        p_immunity = immunity(input$country, as.Date(input$date_event), ifelse(as.Date(input$date_event)<=Sys.Date(), as.Date(input$date_event), Sys.Date()-1), VACCINATIONS, ODDS_RATIO_VACCINE_DOSE1, ODDS_RATIO_VACCINE_DOSE2, plot=FALSE)
        p_immunity_worst = immunity(input$country, as.Date(input$date_event), ifelse(as.Date(input$date_event)<=Sys.Date(), as.Date(input$date_event), Sys.Date()-1), VACCINATIONS, ODDS_RATIO_VACCINE_DOSE1_WORST, ODDS_RATIO_VACCINE_DOSE2_WORST, plot=FALSE)
        p_immunity_best = immunity(input$country, as.Date(input$date_event), ifelse(as.Date(input$date_event)<=Sys.Date(), as.Date(input$date_event), Sys.Date()-1), VACCINATIONS, ODDS_RATIO_VACCINE_DOSE1_BEST, ODDS_RATIO_VACCINE_DOSE2_BEST, plot=FALSE)
        beta_immunity_params = beta.parms.from.quantiles(c(p_immunity_worst,  p_immunity_best))
        n_susc  <- (sapply(1:B, function(b){
        #### ASSUME FIRST DOSE BECOMES FULLY VACCINATED AFTER FOUR WEEKS
            if(p_immunity<0.5){
              max(0,N_TOT - rpois(1,sum(rbeta(N_TOT, beta_immunity_params$a, beta_immunity_params$b))))
            }else{
              max(0, rpois(1,sum(rbeta(N_TOT, beta_immunity_params$b, beta_immunity_params$a))))
            }
            #sum(sapply(1:N_TOT, function(u){rbinom(1, 1,1-rbeta(1, beta_immunity_params$a, beta_immunity_params$b))}))
          }))
      }else{
        n_susc = rep(N_TOT,B)
      }
      
      nb_people_vulnerable_at_the_event <- mean(n_susc)
      ##########################################
      ##########################################
      ###### Step 2: TRANSMISSION DYNAMICS
      ##########################################
      ##########################################
      
      ##########################################
      ###### Step 2.a: compute room parameters for aerosolization
      ##########################################
      incProgress(3/D, detail = "Computing Room parameters for aerosolization ")
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
      
     

      ##########################################
      ##########################################
      ###### Step 3: RESULTS
      ##########################################
      ##########################################
      
      #### Step 3: compute probability of infecting people and adverse outcomes using MCMC simulations
      
      
      incProgress(4/D, detail = "Computing infection distribution")
      nb_infective_people = apply(df_sample, 1, function(x){rbinom(1,N_TOT, sum(x ))})
      
      proba_dist_infectiousness <- sapply(1:min(100, N_TOT), function(n){dpois(n,  N_TOT * mean(apply(df_sample,1,sum)) )})
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
      
      L  = min(max(which(proba_dist_infectiousness > 1e-15)), N_TOT)
      
      
      n_infections <- sapply(1:B, function(b){
        t(sapply(1:L, function(n){
          min(n_susc[b], rpareto(1,n_susc[b] * proba_infection[n] *0.16/1.16, 1.16))
        }))
      })
      

                                      
      
      incProgress(5/D, detail = "Computing summaries for the number of people infected")
      res = data.frame("Average Number of Transmissions" = t(as.matrix(apply(n_infections,1,mean)))%*% as.matrix(proba_dist_infectiousness[1:L]),
                       "Q97.5 Number of Transmissions" = t(as.matrix(apply(n_infections,1,quantile, 0.975)))%*% as.matrix(proba_dist_infectiousness[1:L]),
                       "Q2.5 Number of Transmissions" = t(as.matrix(apply(n_infections,1,quantile, 0.025)))%*% as.matrix(proba_dist_infectiousness[1:L]),
                       "Q99 Number of Transmissions" = t(as.matrix(apply(n_infections,1,quantile, 0.99)))%*% as.matrix(proba_dist_infectiousness[1:L]),
                       "Q1 Number of Transmissions" = t(as.matrix(apply(n_infections,1,quantile, 0.01)))%*% as.matrix(proba_dist_infectiousness[1:L]),
                       "Q50 Number of Transmissions" = t(as.matrix(apply(n_infections,1,quantile, 0.5)))%*% as.matrix(proba_dist_infectiousness[1:L]))
      
 
    
      ##########################################
      ##########################################
      ###### Step 4: COMPARISON WITH BASELINE
      ##########################################
      ##########################################
      incProgress(6/D, detail = "Simulating Infections")   
      n_simulated_infections  = data.frame("n" = (sapply(1:B, function(b){
        ifelse(nb_infective_people[b] == 0, 0, min(n_susc[b], round(rpareto(1, (n_susc[b] - nb_infective_people[b]) * proba_infection[nb_infective_people[b]] *0.16/1.16, 1.16))))
      })), type="Event", date=input$date_event)
      n_simulated_infections = n_simulated_infections %>% filter( n <= res$Q99.Number.of.Transmissions)
      
      incProgress(7/D, detail = "Computing Baseline rates") 
      
      prev_day_event =  future_prevalence_df %>% filter(Date_of_cases == input$date_event)
      proba_baseline <- t(sapply(1:B, function(b){
        max(0,rnorm(n=1, mean=prev_day_event$prevalence,
                    sd=prev_day_event$sd_prevalence))
      }))
      
      n_simulated_infections <- rbind(n_simulated_infections,
      data.frame("n" = (sapply(1:B, function(b){
        ifelse(nb_infective_people[b] == 0, 0, rbinom(1,n_susc[b],proba_baseline[b]))
      })), type="No Event", date=input$date_event))
      
      
     
      
      quantiles= c(0.99, 0.975, 0.5, 0.025, 0.01)
      infections_baseline <- data.frame(rbindlist(lapply(1:B, function(b){
         a = data.frame(t(c(sapply(quantiles, function(x){qbinom(x, n_susc[b],proba_baseline[b])}), proba_baseline[b] *  n_susc[b])))
         colnames(a)= c("Q99.Number.of.Transmissions","Q97.5.Number.of.Transmissions", "Q50.Number.of.Transmissions", "Q2.5.Number.of.Transmissions" , "Q1.Number.of.Transmissions" , "Average.Number.of.Transmissions")
         return(a)
        })))
      
      infections_baseline = infections_baseline[colnames(res)]
      res = rbind(res, apply(infections_baseline , 2, mean))
      res["Type"] = c("Event", "No Event")
      

    })
    

    return(list(res=res,
                proba_baseline = proba_baseline,
                summary_baseline = apply(infections_baseline , 2, mean),
                prevalence_df=future_prevalence_df,
                infectiousness_all = infectiousness_all,
                n = input$n_people,
                n_simulated_infections  =  n_simulated_infections,
                nb_people_infected = nb_people_infected,
                nb_people_infectious_at_the_event =nb_people_infectious_at_the_event ,
                nb_people_detected = nb_people_detected,
                nb_people_vulnerable_at_the_event=nb_people_vulnerable_at_the_event,
                country=input$country
                ))
  })
  
  
  output$contents <- renderTable({
    
    x =  dataInput()
    # set up cut-off values 
    breaks <- c(0,2,4,6,8,10,15,20,25,30,40,50,60,70,80,90,100,200, x$n)
    # specify interval/bin labels
    tags_x <- c("[0-2)","[2-4)", "[4-6)", "[6-8)", "[8-10)", "[10-15)", "[15-20)","[16-18)", "[18-20)",
                "[20-25)","[25-30)", "[30-40)", "[40-50)", "[50-60)", "[60-70)","[70-80)", "[80-90)","[90-100)","[100-)")
    # bucketing values into bins
    
    df <- data.frame("Setting"= c("Event", "Baseline"))
    group_tags <- cut(unlist(x$n_simulated_infections%>% filter(type=="Event") %>% select(n)), 
                      breaks=breaks, 
                      include.lowest=TRUE, 
                      right=FALSE)
    gp2 = cut(unlist(x$n_simulated_infections%>% filter(type=="No Event") %>% select(n)), 
              breaks=breaks, 
              include.lowest=TRUE, 
              right=FALSE)
    options(digits=7)
    df <- cbind(df, rbind(summary(group_tags)/sum(summary(group_tags)),
                          summary(gp2)/sum(summary(gp2))))
    return(df)
    
    
  })
  
  output$disclaimer = renderUI({
    tags$div(
      tags$p("This prototype calculator estimates the number of infections for an event using information regarding the country where the event occurs, 
             the characteristics of the event, the participants and the screening protocol. We choose to output here the distribution of probable number of transmissions, because it allows us to capture more information (average expected number of transmissions, worst case scenarios, etc.).
             In particular, our calculator allows predictions for an event ahead of time, with a screening via antigen D days before the event and a symptom check at entry. This allows us to compute the probability that an infected participant falls through the screening protocol and is actually infectious at the eent (see plot below)."), 
      tags$img(src = "final_plot_infectiousness.png", height = 300, width =600),
      tags$p("You can find more information on the algorithms we use to come up with these numbers in our white paper at the following link."), 
      tags$br(),
      tags$b("Note: This risk estimator is not exact!!! Because COVID-19's transmission is still ill-understood and characterized, this prototype is based on a number of assumptions and simplifications. It should not be understood as exact, but rather, as an estimate of the risk's magnitude, and as a tool to compare the efficiency of mask wearing and vaccination rates."), 
      tags$br(),
      tags$p("To use this calculator, please fill out all of the questions on the left hand side, then click on the 'Result' tab at the top of the screen to see our predictions."),
      tags$b("We are using simulations to compute the uncertainty around our prediction, so the computations might take a little while."),
      tags$br(),
      tags$p("We are a group of researchers around the world, still currently in the process of developing this tool. As a consequence, we are not professional developpers --- so please do forgive any bugs that you might encounter and address your questions and comments to riskmanagement4capacity@gmail.com."),
      tags$p("We do not save any of the information you input."),
      tags$br()
        )
  })
  
  
  
  
  output$plotinput =renderPlot({
    x =  dataInput()

    pt2 <- ggplot(x$prevalence_df)+
      geom_line(aes(x=Date_of_cases, y=1e6 * prevalence), colour="red")+
      geom_ribbon(aes(x=Date_of_cases, ymin=1e6 * ymin, ymax=1e6*ymax), colour="grey", alpha=0.5)+
      theme_bw()+ xlab("Date") + ylab("Incidence (per Million)")+
      ggtitle("Predicted Incidence (per Million) Per Region")+scale_y_log10()


      
    
    pt3<-ggplot(x$infectiousness_all, aes(x=`Date.of.Infection`)) +
      geom_ribbon(aes(x = `Date.of.Infection`, ymin=100 * `p2`, ymax=100 * p97, fill="Escapes Screening"),alpha=0.3)+
      geom_ribbon(aes(x = `Date.of.Infection`, ymin=`infectiousness_q025`, ymax=infectiousness_q975, fill="Infectious"),alpha=0.3)+
      geom_ribbon(aes(x = `Date.of.Infection`, ymin=infectiousness_event_low , ymax=infectiousness_event_high, fill="Infectious at the event"),alpha=0.3)+
      scale_fill_manual(breaks = c("Escapes Screening", "Infectious", "Infectious at the event"), values=c( "gray16", "gray61", "red")) +
      geom_line(aes(x = `Date.of.Infection`, y=100 * `p` , colour="Escapes Screening"), size=1.2)+theme_bw() +
      geom_line(aes(x = `Date.of.Infection`, y=infectiousness, colour = "Infectious"), size=1.2)+
      geom_line(aes(x = `Date.of.Infection`, y=infectiousness * p, colour = "Infectious at the event"), size=1.2)+
      scale_colour_manual(breaks = c("Escapes Screening","Infectious", "Infectious at the event"),  values=c("black", "gray60","red"))+
      theme(legend.text=element_text(size=16))+labs(colour="Average Probability", fill="95% CI Range", size=14) + 
      theme(text = element_text(size=20),
            axis.text.x = element_text(angle=90, hjust=1)) + 
      xlab("Date of Infection\n (with respect to Event Date)") + 
      ylab("Probability (%)")
  

    pt4 <- ggplot(x$n_simulated_infections)+
       geom_density(aes(x=n, fill=type), adjust = 2, alpha=0.4)+
       theme_bw() + 
       ylab("Density") + xlab("Number of new cases/ transmission at the event")
       
    pt5 <- ggplot(x$n_simulated_infections)+
         geom_boxplot(aes(x=type,y=n, fill=type)) +
       theme_bw() + 
         ylab("Number of new cases\n/ transmission at the event") + xlab("Situation")
   
    grid.arrange(pt2,pt3,pt4,pt5, nrow=4)
    
  })
  
  
  

}
shinyApp(ui, server)
