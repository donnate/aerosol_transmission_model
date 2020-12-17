list.of.packages  <-  c("ggplot2", "doBy", "tidyverse","slider")
#new.packages  <-  list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.package, lib="~/R_libs")

library(ggplot2)
#install.packages('doBy', repos='http://cran.us.r-project.org')
library(doBy)
library(tidyverse)
#install.packages('slider', lib="~/R_libs", repos='http://cran.us.r-project.org')
library(pbapply)
library('slider')

# set working directory (on Freddy's computer)
setwd("C:/Users/fb370/Documents/Non-PostDoc work/Codi/aerosol_transmission_model")


compute_prevalence <- function(event_date, country, nb_curves=20){

  EVENT_DATE = as.Date(event_date)
  MIN_DATE = as.Date("2020-03-01")
  MAX_DATE = as.Date("2021-03-01")
  COUNTRY = country
  PERIOD_FOR_FITTING = 14
  PERIOD_FOR_PREDICTING = max(23, min(EVENT_DATE - Sys.Date(), MAX_DATE - Sys.Date()))
  NB_OF_CASE_CURVES = nb_curves
  DAYS_TO_EVENT = 0
  TIME_TO_SYMPTOM_ONSET = 5
  TIME_FROM_SYMPT_ONSET_TO_TEST_RESULT = 4
  PROPORTION_CASES_DETECTED = 1
  PREDICTION_TYPE = ifelse(EVENT_DATE >= Sys.Date(), "future", "past")
    
    
  # Convert date to numeric
  COUNTRY_DATA <- read.csv(file="https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv", header=T, sep=",")
  COUNTRY_DATA$date <- as.numeric(as.Date(COUNTRY_DATA$date, "%Y-%m-%d"))
  # Select the smoothed new cases per million in your country of interest
  
  print(PREDICTION_TYPE)
  filename=paste0("prevalence_", COUNTRY, "_data.csv")
  if (PREDICTION_TYPE == "future"){
    chosen_location_data <- COUNTRY_DATA %>% 
      dplyr::filter((location==COUNTRY) & (date >=  Sys.Date() - PERIOD_FOR_FITTING ) &  (date < Sys.Date())) %>%
      dplyr::select(new_cases_smoothed_per_million)
    # Select all the other data from other countries which might fit your country of interest
    country_data_historic <-  COUNTRY_DATA%>% 
      dplyr::filter( (date > MIN_DATE) &  (date <= Sys.Date() - PERIOD_FOR_PREDICTING))
    
    # Convert the data to wide format by location name
    country_data_historic_wide <- country_data_historic %>% 
      dplyr::select(date, location, new_cases_smoothed_per_million) %>% 
      pivot_wider(.,names_from=location, values_from=new_cases_smoothed_per_million)
    
    # Define the difference function
    Difference_function <- function(data) {
      if(sum(is.na(data))>0){
        return(NA)
        ### might want to 
      }
      sum((data-chosen_location_data$new_cases_smoothed_per_million)^2)
  
    }
    # Calculate the differences matrix
    #pboptions(type = "txt", style = 1, char = "=")
    Differences_matrix <- apply(country_data_historic_wide[,2:ncol(country_data_historic_wide)], 
                                MARGIN=2, 
                                FUN=function(x){ 
                                  sapply((PERIOD_FOR_FITTING):length(x),
                                         function(d){Difference_function(x[(d-PERIOD_FOR_FITTING + 1):d])})
                                  }
                            )
    # Format, date and melt the differences matrix to a long dataframe
  
    Differences_matrix = data.frame(Differences_matrix)
    Differences_matrix$date = unlist(country_data_historic_wide[PERIOD_FOR_FITTING:nrow(country_data_historic_wide),1],
                                     use.names = FALSE)
    diff_vec = reshape2::melt(as_tibble(Differences_matrix), id.vars=c("date"))
                
    
    # Find and select the end time points of the closes case curves
    closest_case_curves <- which.minn(diff_vec$value,NB_OF_CASE_CURVES)  # decreasing = FALSE, index.return=TRUE)$ix[1:NB_OF_CASE_CURVES]
    diff_vec2 = diff_vec[closest_case_curves, ]
    
    # Find all the case data before and after the last date of the best fit curves
    full_closest_case_curves = sapply(1:NB_OF_CASE_CURVES, function(x){
      print(gsub(".", " ", toString(diff_vec2$variable[x]), fixed = TRUE))
      unlist(COUNTRY_DATA %>% 
               filter(location == gsub(".", " ", toString(diff_vec2$variable[x]), fixed=TRUE),
                      date >= diff_vec2$date[x] - PERIOD_FOR_FITTING,
                      date <= diff_vec2$date[x]  +  PERIOD_FOR_PREDICTING) %>%
               dplyr::select(new_cases_smoothed_per_million), use.names = FALSE)
      })
    
    # Convert to dataframe, add a time column and melt to long format
    full_closest_case_curves <- data.frame(full_closest_case_curves)
    full_closest_case_curves$time = seq(from=-PERIOD_FOR_FITTING,
                                        to=PERIOD_FOR_PREDICTING,
                                        by=1)
    melted_case_curves <- reshape2::melt(full_closest_case_curves, 
                                         id.vars="time")
    
    # Calculate the mean and standard deviation at each timepoint
    Summarised_case_predictions <- melted_case_curves %>% 
      group_by(time) %>% 
      dplyr::summarise(prevalence=mean(value/1e6), sd_prevalence=sd(value/1e6))
    
    # Parametrise the delay between infections and cases being reported
    Infection_to_test_result_delay <- TIME_TO_SYMPTOM_ONSET + TIME_FROM_SYMPT_ONSET_TO_TEST_RESULT
    
    #Calculate the date of the cases
    Summarised_case_predictions$Date_of_cases=as.Date(Summarised_case_predictions$time + 
                                                      EVENT_DATE,
                                                      origin = "1970-01-01")
    
    # Create a dataframe of infection
    Infections_df <- data.frame(Date_of_infection=Summarised_case_predictions$Date_of_cases - Infection_to_test_result_delay, 
                              Infection_prevalence=Summarised_case_predictions$prevalence/PROPORTION_CASES_DETECTED, 
                              sd_Infection_prevalence=Summarised_case_predictions$sd_prevalence/PROPORTION_CASES_DETECTED)
    
    # Reverse sort by date of infection
    sorted_case_predictions <- Infections_df[order(Infections_df$Date_of_infection, decreasing=T), ]
    # save as csv file
    #write.csv(x=sorted_case_predictions, file=filename)
    #print(sorted_case_predictions)
  }else{
    print("yata")
    chosen_location_data <- COUNTRY_DATA %>% 
      dplyr::filter((location==COUNTRY) & (date >=  EVENT_DATE - PERIOD_FOR_FITTING ) &  (date <= min(EVENT_DATE + PERIOD_FOR_PREDICTING, Sys.Date()))) %>%
      dplyr::select(new_cases_smoothed_per_million, date)
    names(chosen_location_data) <- c("Infection_prevalence","Date_of_cases")
    chosen_location_data$Infection_prevalence = chosen_location_data$Infection_prevalence/1e6
    chosen_location_data$Date_of_cases=as.Date(chosen_location_data$Date_of_cases,
                                                      origin = "1970-01-01")
    chosen_location_data["sd_Infection_prevalence"] = 0
    chosen_location_data$Date_of_infection = as.Date(chosen_location_data$Date_of_cases - TIME_FROM_SYMPT_ONSET_TO_TEST_RESULT - TIME_TO_SYMPTOM_ONSET,
                                                      origin = "1970-01-01")
    sorted_case_predictions = chosen_location_data
  }
  return(sorted_case_predictions)
}

