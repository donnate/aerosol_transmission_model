list.of.packages <- c("ggplot2", "doBy", "tidyverse","slider")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(ggplot2)
library(doBy)
library(tidyverse)
library(slider)


#Load the Our World In Data dataset

#country_data<-read.csv("owid-covid-data.csv", header=T)
country_data<-read.csv(file="https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv", header=T, sep=",")

# Convert date to numeric
country_data$date<-as.numeric(as.Date(country_data$date, "%Y-%m-%d")) 


#Defining the variables

Selected_country="Estonia"
Period_for_fitting=14
Period_for_predicting=23
Number_of_case_curves=20


# Select the smoothed new cases per million in your country of interest
chosen_location_data<-filter(country_data, location==Selected_country & date>max(date)-Period_for_fitting) %>% 
  dplyr::select(new_cases_smoothed_per_million)

# Select all the other data from other countries which might fit your country of interest
country_data_historic<-filter(country_data, date<max(date)-Period_for_predicting)

# Convert the data to wide format by location name
country_data_historic_wide<-country_data_historic %>% 
  dplyr::select(date, location, new_cases_smoothed_per_million) %>% 
  pivot_wider(.,names_from=location, values_from=new_cases_smoothed_per_million)

# Define the difference function
Difference_function <- function(data) {
  if(sum(is.na(data))>0){
    return(NA)
  }
  sum((data-chosen_location_data)^2)
}

# Calculate the differences matrix
Differences_matrix<-apply(country_data_historic_wide[,2:ncol(country_data_historic_wide)], MARGIN=2, FUN=function(x){ 
  sapply(Period_for_fitting:length(x), function(d){Difference_function(x[(d-Period_for_fitting):d])})
})

# Format, date and melt the differences matrix to a long dataframe
Differences_matrix = data.frame(Differences_matrix)
Differences_matrix$date = unlist(country_data_historic_wide[Period_for_fitting:nrow(country_data_historic_wide),1], use.names = FALSE)
diff_vec = reshape2::melt(as.tibble(Differences_matrix), id.vars=c("date"))

# Find and select the end time points of the closes case curves
closest_case_curves<-which.minn(diff_vec$value,Number_of_case_curves)
diff_vec2 = diff_vec[closest_case_curves,]

# Find all the case data before and after the last date of the best fit curves
full_closest_case_curves = sapply(1:Number_of_case_curves, function(x){
  print(gsub(".", " ", toString(diff_vec2$variable[x]), fixed = TRUE))
  unlist(country_data %>% 
           filter(location == gsub(".", " ", toString(diff_vec2$variable[x]), fixed=TRUE),
                  date >= diff_vec2$date[x] - Period_for_fitting,
                  date <= diff_vec2$date[x] + Period_for_predicting) %>%
           dplyr::select(new_cases_smoothed_per_million), use.names = FALSE)
})

# Convert to dataframe, add a time column and melt to long format
full_closest_case_curves<-data.frame(full_closest_case_curves)
full_closest_case_curves$time=seq(-Period_for_fitting,Period_for_predicting,by=1)
melted_case_curves<-reshape2::melt(full_closest_case_curves, id.vars="time")

# Calculate the mean and standard deviation at each timepoint
Summarised_case_predictions<-melted_case_curves %>% 
  group_by(time) %>% 
  dplyr::summarise(prevalence=mean(value/1000000), sd_prevalence=sd(value/1000000))

#Calculate the date of infection
Summarised_case_predictions$Date_of_infection=as.Date(Summarised_case_predictions$time+max(country_data$date),origin = "1970-01-01")

# Parametrise the delay between infections and cases being reported
Days_to_event<-0
Time_to_symptom_onset<-5
Time_from_symptom_to_test_result<-4
Infection_to_test_result_delay<-Time_to_symptom_onset+Time_from_symptom_to_test_result
Days_from_infection_to_event_delay<-Infection_to_test_result_delay+Days_to_event

# Calculate the day of infection on event day
Summarised_case_predictions$Days_since_infection=Days_from_infection_to_event_delay-Summarised_case_predictions$time

# Filter to days since infection to less than 14 then sort by days since infection
filtered_case_predictions<-filter(Summarised_case_predictions, Days_since_infection <=14)
sorted_case_predictions<-filtered_case_predictions[order(filtered_case_predictions$Days_since_infection),]
  
# save as csv file
write.csv(x=sorted_case_predictions,file = "chosen_prevalence_data.csv")

