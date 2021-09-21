#### Albert Hall
setwd("~/Dropbox/aerosol_transmission_model/")
source("code_experiments_paper/covid_projectons.R")
#source("proba_symptoms.R")
source("helper_functions.R")
source("vaccination.R")
source("aerosol_functions.R")
source("beta_params.R")
source("covid_case_predictions.R")
source("under_ascertainment_bias.R")
source("screening_efficiency.R")
source("relative_infectiousness.R")
library(ggplot2)

MASK_EFFICIENCY = 0.5  ### 50% is the recommended value
##### Set parameters
input = list("country" =  "United Kingdom",
             "date_event" = as.Date("2021-03-15"),
             "temperature" = 23,
             "UV"=0,
             "RH"=50, 
             "time2event"=2,
             "duration" = 180/60, #### time must be in hours
             "unit" =0,
             "length" =50,
             "width" =50,
             "height"=5,
             "temperature"=20,
             "RH"=40,
             "UV" =3,
             "ventilation"="Daycare",
             "control"=0,
             "activity"="Standing:Loudly speaking",
             "mask"=0.0,
             "prop_mask"=100,
             "mixing"=0,
             "p_symptoms" =0.4,
             "p_lie" = 0.5,
             file1=NULL 
)
B=100000
N_TOT = 5000
volume = 86650
surface = 20000

country = input$country
PERIOD_FOR_FITTING = 7 
NCURVES = 100
MAX_DATE = input$date_event - 28
COUNTRY_DATA <- read.csv(file="https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv", header=T, sep=",")
COUNTRY_DATA$date <- (as.Date(COUNTRY_DATA$date, "%Y-%m-%d"))
VACCINATIONS <- read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations.csv")
VACCINATIONS_US <- read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/us_state_vaccinations.csv")
intersection_col = intersect(colnames(VACCINATIONS), colnames(VACCINATIONS_US))
VACCINATIONS <- rbind(VACCINATIONS[,intersection_col], VACCINATIONS_US[,intersection_col] %>% filter(location != "United States"))



PERIOD_FOR_PREDICTING = as.numeric(as.Date(input$date_event)- as.Date(max(MAX_DATE, na.rm = TRUE), fmt="%Y-%m-%d"))
DECAY = max(0, (7.57+ 
                  1.41* (input$temperature-20.54)/10.66 +
                  0.0218 *(input$RH-45.24)/28.67 + 
                  7.55 *((input$UV*0.185)-50) / 50 +
                  (input$temperature-20.54)/10.66*(input$UV*0.185-50)/50*1.40) *60)  #https://www.dhs.gov/science-and-technology/sars-airborne-calculator

BREATHING_RATE = 0.012 * 60 
DEPOSITION = 0.24
MASK_INHALATION_EFFICIENCY = 0.5
PRESSURE = 0.95
######################################################################
######################################################################
##### Step 1: compute future prevalence
######################################################################
######################################################################


pred_incidence = c()

for (d in c("2021-03-15", "2021-01-15", "2020-08-15")){
  print(d)
  input$date_event = as.Date(d)
  data2fit = COUNTRY_DATA %>% 
    dplyr::filter(location=="United Kingdom", date >  input$date_event- PERIOD_FOR_FITTING, date <= input$date_event) 
  
  res_1 <- predict_prevalence(origin=input$date_event -28, country=country,
                              country_data = COUNTRY_DATA,nb_curves=NCURVES, 
                              distance=Difference_function, 
                              period4predicting=28, period4fitting =PERIOD_FOR_FITTING,
                              distance_type = "MSE")
  
  res_1$res["Date"] = input$date_event
  res_1$output["Date"] = input$date_event
  
  
  prevalence_df = res_1$res
  prevalence_df["region"] = country
  samples_prev = res_1$samples
  samples_prev["region"] = country
  
  
  
  
  
  
  bias<- compute_underascertainment_bias(input$date_event - (PERIOD_FOR_FITTING + 28), "United Kingdom", COUNTRY_DATA, date_max = input$date_event - 28, 
                                         Case_to_death_delay= 21, plot=FALSE)
  bias_corr  = mean(bias$value)
  bias_sd  = sd(bias$value)
  
  future_prevalence_df = prevalence_df %>%
    mutate(prevalence = 1/bias_corr * prevalence,
           sd_prevalence = 1/bias_corr * sd_prevalence,
           Observed = 1/bias_corr *  Observed,
           q975 = 1/bias_corr *  q975,
           q50 = 1/bias_corr *  q50,
           q25 = 1/bias_corr *  q25,
    )
  future_prevalence_df = future_prevalence_df %>% mutate(ymin = prevalence -2*sd_prevalence,
                                                         ymax = prevalence + 2*sd_prevalence)
  future_prevalence_df["ymin"] = sapply(future_prevalence_df["ymin"], function(x){ifelse(x<0,0,x)})
  future_prevalence_df["ymax"] = sapply(future_prevalence_df["ymax"], function(x){ifelse(x>1,1,x)})
  

  
  future_prevalence_df["Date_event"] = input$date_event
  pred_incidence = rbind(future_prevalence_df,pred_incidence)
}


ggplot(pred_incidence)+
  geom_line(data= pred_incidence %>% filter(time>=0), aes(x=Date_of_cases, y=1e6 * prevalence, colour="Expected"), size=1)+
  geom_ribbon(data= pred_incidence %>% filter(time>=0),aes(x=Date_of_cases, ymin=1e6 * q25, ymax=1e6*q975), colour="grey", alpha=0.5)+
  geom_line(aes(x=Date_of_cases, y=1e6 * Observed, colour="Observed"), size=1) + 
  theme_bw()+ xlab("") + ylab("Incidence (per Million)")+
  facet_wrap(.~Date_event, scales = "free") +
  scale_color_manual(values = c( "red", "black")) + 
  theme(legend.text=element_text(size=16))+labs(colour="Data Type", fill="95% Prediction Range", size=14) + 
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1)) 
