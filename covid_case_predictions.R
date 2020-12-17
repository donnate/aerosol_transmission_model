list.of.packages <- c("ggplot2", "doBy", "tidyverse","slider")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(ggplot2)
library(doBy)
library(tidyverse)
library(slider)

# set working directory (on Freddy's computer)
setwd("C:/Users/fb370/Documents/Non-PostDoc work/Codi/aerosol_transmission_model")

#Defining the variables

Selected_country="Estonia"
Period_for_fitting=14
Period_for_predicting=23
Number_of_case_curves=20
Days_to_event<-0
Time_to_symptom_onset<-5
Time_from_symptom_to_test_result<-4
Cases_detected=1
filename="chosen_prevalence_data.csv"

#Load the Our World In Data dataset by first loading the previously saved dataset and then replacing it with a download if it is not up to date

country_data<-read.csv("owid-covid-data.csv", header=T)
if(max(as.numeric(as.Date(country_data$date, "%Y-%m-%d"))) < as.numeric(Sys.Date())-1){
  country_data=read.csv(file="https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv", header=T, sep=",")
} else {
  country_data=country_data
}
write.csv(x=country_data, file = "owid-covid-data.csv", row.names=FALSE)

# Convert date to numeric
country_data$date<-as.numeric(as.Date(country_data$date, "%Y-%m-%d"))

# Add a new column of log-transformed data
country_data$log_new_cases_smoothed_per_million=log(country_data$new_cases_smoothed_per_million)

# Select the smoothed new cases per million in your country of interest
chosen_location_data<-filter(country_data, location==Selected_country & date>max(date)-Period_for_fitting) %>% 
  dplyr::select(log_new_cases_smoothed_per_million)

# Select all the other data from other countries which might fit your country of interest
country_data_historic<-filter(country_data, date<max(date)-Period_for_predicting)

# Convert the data to wide format by location name
country_data_historic_wide<-country_data_historic %>% 
  dplyr::select(date, location, log_new_cases_smoothed_per_million) %>% 
  pivot_wider(.,names_from=location, values_from=log_new_cases_smoothed_per_million)

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
           dplyr::select(log_new_cases_smoothed_per_million), use.names = FALSE)
})

# Convert to dataframe, add a time column and melt to long format
full_closest_case_curves<-data.frame(full_closest_case_curves)
full_closest_case_curves$time=seq(-Period_for_fitting,Period_for_predicting,by=1)
melted_case_curves<-reshape2::melt(full_closest_case_curves, id.vars="time")

# Calculate the mean and standard deviation at each timepoint
Summarised_case_predictions<-melted_case_curves %>% 
  group_by(time) %>% 
  dplyr::summarise(prevalence=mean(exp(value)/1000000), sd_prevalence=sd(exp(value)/1000000))

# Parametrise the delay between infections and cases being reported
Infection_to_test_result_delay<-Time_to_symptom_onset+Time_from_symptom_to_test_result

#Calculate the date of the cases
Summarised_case_predictions$Date_of_cases=as.Date(Summarised_case_predictions$time+max(country_data$date),origin = "1970-01-01")

# Create a dataframe of infection
Infections_df<-data.frame(Date_of_infection=Summarised_case_predictions$Date_of_cases-Infection_to_test_result_delay, Infection_prevalence=Summarised_case_predictions$prevalence/Cases_detected, sd_Infection_prevalence=Summarised_case_predictions$sd_prevalence/Cases_detected)

# Reverse sort by date of infection
sorted_case_predictions<-Infections_df[order(Infections_df$Date_of_infection, decreasing=T),]
  
# save as csv file
write.csv(x=sorted_case_predictions,file = filename)

# plotting the individual curves
ggplot(melted_case_curves, aes(x=time, y=value, colour=variable))+
  geom_line()+
  theme_classic()

# plotting the summarised curves
ggplot(Summarised_case_predictions, aes(x=time, y=prevalence*1000000))+
  geom_line()+
  geom_errorbar(aes(ymin=10^6*(prevalence-sd_prevalence), ymax=10^6*(prevalence+sd_prevalence)))+
  theme_classic()+
  scale_y_continuous(limits=c(0,NA), labels = scales::comma)+
  labs(x="Time (Days)", y="COVID cases per million")

# read in the ifr data
ifr_data=read.csv("country_ifr_data.csv", header=T)

# merging the country case data with IFR data
cases_and_ifr=merge(country_data, ifr_data, by="location")

# calculating mean ifr
cases_and_ifr$mean_ifr<-rowMeans(cases_and_ifr[,c("ENE_COVID","COVID_US_CDC","COVID_Verity","COVID_Levin")], na.rm=TRUE)

# calculate the 'true cases' from ifr and death rate before time shifting
cases_and_ifr$estimated_cases_per_million<-cases_and_ifr$new_deaths_smoothed_per_million/cases_and_ifr$mean_ifr*100

# plot estimated vs actual cases
ggplot(cases_and_ifr, aes(x=new_cases_smoothed_per_million, y=estimated_cases_per_million))+
  geom_point()+
  theme_classic()+
  geom_smooth(method="lm")
