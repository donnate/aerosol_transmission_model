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
library(projections)


#Defining the variables
#Load the Our World In Data datas
#COUNTRY_DATA <- read.csv("owid-covid-data.csv", header=T)

# Define the difference function
Difference_function <- function(data,z, w=NULL) {
  if(sum(is.na(data))>0){
    return(NA)
    ### might want to get rid of these
  }
  if(is.null(w)){
    w=rep(1, length(z))/length(z)
  }
  #print(w)
  sum(w*(data-z)^2)
}

correlation_function <- function(data, z, w=NULL) {
  if(is.null(w)){
    w=rep(1, length(z))
  }
  if(sum(is.na(data))>0){
    return(NA)
    ### might want to get rid of these
  }else{
    if( (sum(z^2) ==0)){
      return(sum(1-cor(data,runif(length(z), 0, 0.1))))
    }else{
      if( (sum(data^2) ==0)){
        return((1-cor(w*z,runif(length(data), 0, 0.1))))
      }else{
        return((1-cor(data,w*z)))
      }
    }
  }
  
}

predict_prevalence <- function(origin, country="United Kingdom",
                               country_data = NULL,nb_curves=100, 
                               distance=Difference_function, 
                               period4predicting=NULL, period4fitting =NULL, 
                               div=2, distance_type="MSE", weights=NULL, agg=mean){
  ORIGIN = origin
  PERIOD_FOR_FITTING = ifelse(is.null(period4fitting), 28, period4fitting)
  PERIOD_FOR_PREDICTING = ifelse(is.null(period4predicting), 28,period4predicting )
  # Parametrise the delay between infections and cases being reported
  MIN_DATE = as.Date("2020-03-01")
  MAX_DATE = as.Date("2021-03-01")
  COUNTRY = country
  
  NB_OF_CASE_CURVES = nb_curves
  
  # Convert date to numeric
  if (is.null(country_data)){
    COUNTRY_DATA <- read.csv(file="https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv", header=T, sep=",")
    COUNTRY_DATA$date <- (as.Date(COUNTRY_DATA$date, "%Y-%m-%d"))
  }else{
    COUNTRY_DATA <- country_data
  }
  
  # Select the smoothed new cases per million in your country of interest
  
  #print(PREDICTION_TYPE)
  filename=paste0("prevalence_", COUNTRY, "_data.csv")
  #if (PREDICTION_TYPE == "future"){

  chosen_location_data <- COUNTRY_DATA %>% 
      dplyr::filter((location==COUNTRY) & (date >  ORIGIN - PERIOD_FOR_FITTING  ) &  (date <= ORIGIN)) %>%
      dplyr::select(new_cases_smoothed_per_million)
  print(nrow(chosen_location_data))
    if(nrow(chosen_location_data)<1){
      return(NULL)
    }else{
    future <- COUNTRY_DATA %>% 
      dplyr::filter((location==COUNTRY) & (date >  ORIGIN) &  (date <= ORIGIN + PERIOD_FOR_PREDICTING)) %>%
      dplyr::select(new_cases_smoothed_per_million)
    print(nrow(future))
    # Select all the other data from other countries which might fit your country of interest
    country_data_historic <-  COUNTRY_DATA%>% 
      dplyr::filter( (date > MIN_DATE + 1 ) &  (date <= ORIGIN - 1 - PERIOD_FOR_PREDICTING)) 
    
    # Convert the data to wide format by location name
    country_data_historic_wide <- country_data_historic %>% 
      dplyr::select(date, location, new_cases_smoothed_per_million) %>% 
      pivot_wider(.,names_from=location, values_from=new_cases_smoothed_per_million)
    
    
    # Calculate the differences matrix
    #pboptions(type = "txt", style = 1, char = "=")
    # Differences_matrix <- sapply(PERIOD_FOR_FITTING:nrow(country_data_historic_wide), 
    #                             FUN=function(i){ 
    #                               sapply(2:ncol(country_data_historic_wide),
    #                                     function(j){
    #                                       Difference_function(country_data_historic_wide[(i-PERIOD_FOR_FITTING + 1):i, j])})
    #                               }
    #                         )
    Differences_matrix <- apply(country_data_historic_wide[,2:ncol(country_data_historic_wide)],
                                MARGIN=2,
                                FUN=function(x){
                                  sapply(seq(PERIOD_FOR_FITTING,nrow(country_data_historic_wide),floor(PERIOD_FOR_FITTING/div)),
                                         function(d){
                                           #print(x[(d-PERIOD_FOR_FITTING + 1):d,1])
                                           distance(x[(d-PERIOD_FOR_FITTING + 1):d],
                                                                         chosen_location_data$new_cases_smoothed_per_million,
                                                    w=weights)})
                                }
    )
    # Format, date and melt the differences matrix to a long dataframe
    
    Differences_matrix = data.frame(Differences_matrix)
    print(nrow(Differences_matrix))
    colnames(Differences_matrix) <- colnames(country_data_historic_wide)[2:ncol(country_data_historic_wide)]
    Differences_matrix["date"] = unlist(country_data_historic_wide[seq(PERIOD_FOR_FITTING,nrow(country_data_historic_wide),floor(PERIOD_FOR_FITTING/div)),"date"],
                                        use.names = FALSE)
    
    diff_vec = reshape2::melt(as_tibble(Differences_matrix), id.vars=c("date"))
    
    
    # Find and select the end time points of the closes case curves
    closest_case_curves <- which.minn(diff_vec$value,NB_OF_CASE_CURVES)  # decreasing = FALSE, index.return=TRUE)$ix[1:NB_OF_CASE_CURVES]
    diff_vec2 = diff_vec[closest_case_curves, ]
    diff_vec2$date = as.Date(diff_vec2$date, origin = "1970-01-01")
    # Find all the case data before and after the last date of the best fit curves
    if (distance_type == "correlation"){
      full_closest_case_curves = sapply(1:NB_OF_CASE_CURVES, function(x){
        #print(c(gsub(".", " ", toString(diff_vec2$variable[x]), fixed = TRUE), diff_vec2$date[x], as.Date(as.numeric(diff_vec2$date[x]), origin = "1970-01-01")))
        t = unlist(COUNTRY_DATA %>% 
                     filter(location == gsub(".", " ", toString(diff_vec2$variable[x]), fixed=TRUE),
                            date >  diff_vec2$date[x] - PERIOD_FOR_FITTING,
                            date <= diff_vec2$date[x]  +  PERIOD_FOR_PREDICTING) %>%
                     dplyr::select(new_cases_smoothed_per_million), use.names = FALSE)
        t = mean(chosen_location_data$new_cases_smoothed_per_million) + max(0.01,sd(chosen_location_data$new_cases_smoothed_per_million))/max(sd(t), 0.001) * (t-mean(t))
        return(t)
      })
    }else{
      full_closest_case_curves = sapply(1:NB_OF_CASE_CURVES, function(x){
       # print(c(gsub(".", " ", toString(diff_vec2$variable[x]), fixed = TRUE), diff_vec2$date[x], as.Date(as.numeric(diff_vec2$date[x]), origin = "1970-01-01")))
        t = unlist(COUNTRY_DATA %>% 
                     filter(location == gsub(".", " ", toString(diff_vec2$variable[x]), fixed=TRUE),
                            date >  diff_vec2$date[x] - PERIOD_FOR_FITTING,
                            date <= diff_vec2$date[x]  +  PERIOD_FOR_PREDICTING) %>%
                     dplyr::select(new_cases_smoothed_per_million), use.names = FALSE)
        t = t + (chosen_location_data$new_cases_smoothed_per_million[PERIOD_FOR_FITTING] - t[PERIOD_FOR_FITTING])
        return(sapply(t, function(x){(max(0,x))}))
      })
    }
    # Convert to dataframe, add a time column and melt to long format
    full_closest_case_curves <- data.frame(full_closest_case_curves)
    full_closest_case_curves["time"] = seq(from=-PERIOD_FOR_FITTING + 1,
                                           to=PERIOD_FOR_PREDICTING,
                                           by=1)
    melted_case_curves <- reshape2::melt(full_closest_case_curves, 
                                         id.vars="time")
    
    #print(ggplot(melted_case_curves, aes(x=time, y=value, group = variable)) +
    #  geom_line() + geom_line())
    # Calculate the mean and standard deviation at each timepoint
    Summarised_case_predictions <- melted_case_curves %>% 
      group_by(time) %>% 
      dplyr::summarise(prevalence=agg(value/1e6), sd_prevalence=sd(value/1e6),
      q50 = quantile(value/1e6, 0.5),
      q975 = quantile(value/1e6, 0.98),
      q25 = quantile(value/1e6, 0.02))
    
    
    
    #Calculate the date of the cases
    Summarised_case_predictions$Date_of_cases=as.Date(Summarised_case_predictions$time + 
                                                        ORIGIN,
                                                      origin = "1970-01-01")
    if (nrow(future) == 28){
      Summarised_case_predictions["Observed"] = c(chosen_location_data$new_cases_smoothed_per_million/1e6, future$new_cases_smoothed_per_million/1e6)
    }else{
      Summarised_case_predictions["Observed"] = c(chosen_location_data$new_cases_smoothed_per_million/1e6, rep(NA,28))
    }
    
      
      # Create a dataframe of infection
      #Infections_df <- data.frame(Date_of_infection=Summarised_case_predictions$Date_of_cases - Infection_to_test_result_delay, 
      #                         Infection_prevalence=Summarised_case_predictions$prevalence/PROPORTION_CASES_DETECTED, 
      #                          sd_Infection_prevalence=Summarised_case_predictions$sd_prevalence/PROPORTION_CASES_DETECTED)
      
      # Reverse sort by date of infection
      #sorted_case_predictions <- Infections_df[order(Infections_df$Date_of_infection, decreasing=T), ]
      # save as csv file
      #write.csv(x=sorted_case_predictions, file=filename)
      #print(sorted_case_predictions)
    return(list(res = Summarised_case_predictions, output =  melted_case_curves))
  
    }
}


##### C

