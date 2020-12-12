list.of.packages  <-  c("ggplot2", "doBy", "tidyverse","slider")
new.packages  <-  list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(ggplot2)
library(doBy)
library(tidyverse)
library(slider)


#Defining the variables
#Load the Our World In Data datas
#COUNTRY_DATA <- read.csv("owid-covid-data.csv", header=T)


predict_prevalence  <-  function( Selected_country="Estonia",
                                Period_for_fitting=14,
                                Period_for_predicting=23,
                                Number_of_case_curves=20,
                                Days_to_event =0,
                                Time_to_symptom_onset=5,
                                Time_from_symptom_to_test_result=4,
                                Cases_detected=1){
  
  
  filename=paste0("prevalence_",Selected_country, "_data.csv")
  # Convert date to numeric
  COUNTRY_DATA <- read.csv(file="https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv", header=T, sep=",")
  COUNTRY_DATA$date <- as.numeric(as.Date(COUNTRY_DATA$date, "%Y-%m-%d")) 
  # Select the smoothed new cases per million in your country of interest
  chosen_location_data <- COUNTRY_DATA %>% 
    dplyr::filter((location==Selected_country) & (date> max(date) - Period_for_fitting)) %>%
    dplyr::select(new_cases_smoothed_per_million)
  # Select all the other data from other countries which might fit your country of interest
  country_data_historic <- filter(COUNTRY_DATA, date< max(date)-Period_for_predicting)
  
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
  Differences_matrix <- apply(country_data_historic_wide[,2:ncol(country_data_historic_wide)], 
                              MARGIN=2, 
                              FUN=function(x){ 
                                sapply((Period_for_fitting):length(x),
                                       function(d){Difference_function(x[(d-Period_for_fitting + 1):d])})
                                }
                              )
  # Format, date and melt the differences matrix to a long dataframe
  Differences_matrix = data.frame(Differences_matrix)
  Differences_matrix$date = unlist(country_data_historic_wide[Period_for_fitting:nrow(country_data_historic_wide),1],
                                   use.names = FALSE)
  diff_vec = reshape2::melt(as_tibble(Differences_matrix), id.vars=c("date"))
  
  # Find and select the end time points of the closes case curves
  closest_case_curves <- which.minn(diff_vec$value,  Number_of_case_curves)
  diff_vec2 = diff_vec[closest_case_curves, ]
  
  # Find all the case data before and after the last date of the best fit curves
  full_closest_case_curves = sapply(1:Number_of_case_curves, function(x){
    print(gsub(".", " ", toString(diff_vec2$variable[x]), fixed = TRUE))
    unlist(COUNTRY_DATA %>% 
             filter(location == gsub(".", " ", toString(diff_vec2$variable[x]), fixed=TRUE),
                    date >= diff_vec2$date[x] - Period_for_fitting,
                    date <= diff_vec2$date[x]  +  Period_for_predicting) %>%
             dplyr::select(new_cases_smoothed_per_million), use.names = FALSE)
    })
  
  # Convert to dataframe, add a time column and melt to long format
  full_closest_case_curves <- data.frame(full_closest_case_curves)
  full_closest_case_curves$time = seq(from=-Period_for_fitting,
                                      to=Period_for_predicting,
                                      by=1)
  melted_case_curves <- reshape2::melt(full_closest_case_curves, 
                                       id.vars="time")
  
  # Calculate the mean and standard deviation at each timepoint
  Summarised_case_predictions <- melted_case_curves %>% 
    group_by(time) %>% 
    dplyr::summarise(prevalence=mean(value/1e6), sd_prevalence=sd(value/1e6))
  
  # Parametrise the delay between infections and cases being reported
  Infection_to_test_result_delay <- Time_to_symptom_onset + Time_from_symptom_to_test_result
  
  #Calculate the date of the cases
  Summarised_case_predictions$Date_of_cases=as.Date(Summarised_case_predictions$time + 
                                                      max(COUNTRY_DATA$date),
                                                    origin = "1970-01-01")
  
  # Create a dataframe of infection
  Infections_df <- data.frame(Date_of_infection=Summarised_case_predictions$Date_of_cases - Infection_to_test_result_delay, 
                            Infection_prevalence=Summarised_case_predictions$prevalence/Cases_detected, 
                            sd_Infection_prevalence=Summarised_case_predictions$sd_prevalence/Cases_detected)
  
  # Reverse sort by date of infection
  sorted_case_predictions <- Infections_df[order(Infections_df$Date_of_infection, decreasing=T), ]
  # save as csv file
  write.csv(x=sorted_case_predictions, file=filename)
  return(sorted_case_predictions)
}


