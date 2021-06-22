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


compute_prevalence <- function(event_date, data2fit, 
                               country_data = NULL,nb_curves=100, 
                               distance=Difference_function, 
                               period4predicting=NULL, period4fitting =NULL, 
                               div=2, distance_type="MSE",plot=FALSE, weights=NULL){
  
  
  PERIOD_FOR_FITTING = ifelse(is.null(period4fitting), min(28, nrow( data2fit)), period4fitting)
  PERIOD_FOR_PREDICTING = ifelse(is.null(period4predicting), 28,period4predicting )
  ORIGIN = ifelse(Sys.Date() > as.Date(event_date), as.Date(event_date) - PERIOD_FOR_PREDICTING , Sys.Date() - 1 )
  # Parametrise the delay between infections and cases being reported
  MIN_DATE = as.Date("2020-03-01")
  MAX_DATE = ORIGIN - PERIOD_FOR_FITTING
  
  NB_OF_CASE_CURVES = nb_curves
  # Convert date to numeric

  if (is.null(country_data)){
    COUNTRY_DATA <- read.csv(file="https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv", header=T, sep=",")
    COUNTRY_DATA$date <- (as.Date(COUNTRY_DATA$date, "%Y-%m-%d"))
  }else{
    COUNTRY_DATA <- country_data
  }

  
  

  # Select all the other data from other countries which might fit your country of interest
  country_data_historic <-  COUNTRY_DATA%>% 
    dplyr::filter( (date > MIN_DATE + 1 ) &  (date <= ORIGIN - 1 - PERIOD_FOR_PREDICTING)) 
  
  # Convert the data to wide format by location name
  country_data_historic_wide <- country_data_historic %>% 
    dplyr::select(date, location, new_cases_smoothed_per_million) %>% 
    pivot_wider(.,names_from=location, values_from=new_cases_smoothed_per_million)
    
  Differences_matrix <- apply(country_data_historic_wide[,2:ncol(country_data_historic_wide)],
                                MARGIN=2,
                                FUN=function(x){
                                  sapply(seq(PERIOD_FOR_FITTING,nrow(country_data_historic_wide),floor(PERIOD_FOR_FITTING/div)),
                                         function(d){
                                           #print(x[(d-PERIOD_FOR_FITTING + 1):d,1])
                                           distance(x[(d-PERIOD_FOR_FITTING + 1):d],
                                                    data2fit$new_cases_smoothed_per_million[(nrow(data2fit) - PERIOD_FOR_FITTING+1):nrow(data2fit)],
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
    
    
    # Find and select the end time points of the closest case curves
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
        t = mean(data2fit$new_cases_smoothed_per_million) + max(0.01,sd(data2fit$new_cases_smoothed_per_million))/max(sd(t), 0.001) * (t-mean(t))
        return(t)
      })
    }else{
      full_closest_case_curves = sapply(1:NB_OF_CASE_CURVES, function(x){
        #print(c(gsub(".", " ", toString(diff_vec2$variable[x]), fixed = TRUE), diff_vec2$date[x], as.Date(as.numeric(diff_vec2$date[x]), origin = "1970-01-01")))
        t = unlist(COUNTRY_DATA %>% 
                     filter(location == gsub(".", " ", toString(diff_vec2$variable[x]), fixed=TRUE),
                            date >  diff_vec2$date[x] - PERIOD_FOR_FITTING,
                            date <= diff_vec2$date[x]  +  PERIOD_FOR_PREDICTING) %>%
                     dplyr::select(new_cases_smoothed_per_million), use.names = FALSE)
        t = t + (data2fit$new_cases_smoothed_per_million[nrow(data2fit)] - t[PERIOD_FOR_FITTING])
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
    #        geom_line() + geom_line())
    # Calculate the mean and standard deviation at each timepoint
    Summarised_case_predictions <- melted_case_curves %>% 
      group_by(time) %>% 
      dplyr::summarise(prevalence=mean(value/1e6), sd_prevalence=sd(value/1e6))
    
    
    
    #Calculate the date of the cases
    Summarised_case_predictions$Date_of_cases=as.Date(Summarised_case_predictions$time + 
                                                        ORIGIN,
                                                      origin = "1970-01-01")
    return(list(res = Summarised_case_predictions, output =  melted_case_curves))
    
}











compute_vaccinations <- function(origin,data2fit, 
                                 type_prediction = "fully vaccinated",
                                 country_data = NULL,nb_curves=100, 
                                 distance=Difference_function, 
                                 period4predicting=28, period4fitting =28, 
                                 div=2, distance_type="MSE"){
  ORIGIN = as.Date(origin)
  PERIOD_FOR_FITTING = period4fitting
  PERIOD_FOR_PREDICTING = period4predicting 
  # Parametrise the delay between infections and cases being reported
  MIN_DATE = as.Date("2020-03-01")
  MAX_DATE = as.Date("2021-03-01")
  
  NB_OF_CASE_CURVES = nb_curves
  
  # Convert date to numeric
  if (is.null(country_data)){
    VACCINATIONS <- read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations.csv")
    VACCINATIONS_US <- read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/us_state_vaccinations.csv") %>% 
      filter(location!="United States")
    intersection_col = intersect(colnames(VACCINATIONS), colnames(VACCINATIONS_US))
    VACCINATIONS <- rbind(VACCINATIONS[,intersection_col], VACCINATIONS_US[,intersection_col])
    VACCINATIONS$date <- (as.Date(VACCINATIONS$date, "%Y-%m-%d"))
  }else{
    VACCINATIONS <- country_data
  }
  
  # Select the smoothed new cases per hundred in your country of interest
  if (type_prediction == "fully vaccinated"){
    VACCINATIONS = VACCINATIONS %>% mutate(vaccinated = people_fully_vaccinated_per_hundred)
  }else{
    VACCINATIONS = VACCINATIONS %>% mutate(vaccinated = people_vaccinated_per_hundred - people_fully_vaccinated_per_hundred)
  }
  

    # Select all the other data from other countries which might fit your country of interest
    country_data_historic <-  VACCINATIONS%>% 
      dplyr::filter( (date > MIN_DATE + 1 ) &  (date <= ORIGIN - 1 - PERIOD_FOR_PREDICTING)) 

    
    country_data_historic_wide <- country_data_historic %>% 
        dplyr::select(date, location, vaccinated) %>%
        pivot_wider(.,names_from=location, values_from=vaccinated)
    
    
    
    Differences_matrix <- apply(country_data_historic_wide[,2:ncol(country_data_historic_wide)],
                                MARGIN=2,
                                FUN=function(x){
                                  sapply(seq(PERIOD_FOR_FITTING,nrow(country_data_historic_wide),floor(PERIOD_FOR_FITTING/div)),
                                         function(d){
                                           #print(x[(d-PERIOD_FOR_FITTING + 1):d,1])
                                           distance(x[(d-PERIOD_FOR_FITTING + 1):d],
                                                    data2fit)})
                                }
    )
    
    Differences_matrix = data.frame(Differences_matrix)
    colnames(Differences_matrix) <- colnames(country_data_historic_wide)[2:ncol(country_data_historic_wide)]
    Differences_matrix["date"] = unlist(country_data_historic_wide[seq(PERIOD_FOR_FITTING,nrow(country_data_historic_wide),floor(PERIOD_FOR_FITTING/div)),"date"],
                                        use.names = FALSE)
    
    diff_vec = reshape2::melt(as_tibble(Differences_matrix), id.vars=c("date"))
    
    
    # Find and select the end time points of the closes case curves
    closest_case_curves <- which.minn(diff_vec$value,NB_OF_CASE_CURVES)  # decreasing = FALSE, index.return=TRUE)$ix[1:NB_OF_CASE_CURVES]
    diff_vec2 = diff_vec[closest_case_curves, ]
    diff_vec2$date = as.Date(diff_vec2$date, origin = "1970-01-01")
    diff_vec2$variable = as.character(diff_vec2$variable )
    # Find all the case data before and after the last date of the best fit curves
    if (distance_type == "correlation"){
      full_closest_case_curves = sapply(1:NB_OF_CASE_CURVES, function(x){
        print(c(diff_vec2$variable[x], diff_vec2$date[x], as.Date(as.numeric(diff_vec2$date[x]), origin = "1970-01-01")))
        t = unlist(VACCINATIONS %>% 
                     filter(location == gsub(".", " ", toString(diff_vec2$variable[x]), fixed=TRUE),
                            date >  diff_vec2$date[x] - PERIOD_FOR_FITTING,
                            date <= diff_vec2$date[x]  +  PERIOD_FOR_PREDICTING) %>%
                     dplyr::select(vaccinated), use.names = FALSE)
        t = mean(data2fit) + max(0.01,sd(data2fit))/max(sd(t), 0.001) * (t-mean(t))
        return(t)
      })
    }else{
      full_closest_case_curves = sapply(1:NB_OF_CASE_CURVES, function(x){
        print(c(as.character(diff_vec2$variable[x]), diff_vec2$date[x], as.Date(as.numeric(diff_vec2$date[x]), origin = "1970-01-01")))
        t = unlist(VACCINATIONS %>% 
                     filter(location == gsub(".", " ", toString(diff_vec2$variable[x]), fixed=TRUE),
                            date >  diff_vec2$date[x] - PERIOD_FOR_FITTING,
                            date <= diff_vec2$date[x]  +  PERIOD_FOR_PREDICTING) %>%
                     dplyr::select(vaccinated), use.names = FALSE)
        
        print(x);print(t)
        #print(c(diff_vec2$variable[x], diff_vec2$date[x], as.Date(as.numeric(diff_vec2$date[x]), origin = "1970-01-01")))
        return(t)
      })
    }
    # Convert to dataframe, add a time column and melt to long format
    full_closest_case_curves <- data.frame(full_closest_case_curves)
    full_closest_case_curves["time"] = seq(from=-PERIOD_FOR_FITTING + 1,
                                           to=PERIOD_FOR_PREDICTING,
                                           by=1)
    melted_case_curves <- reshape2::melt(full_closest_case_curves, 
                                         id.vars="time")
    
    print(ggplot(melted_case_curves, aes(x=time, y=value, group = variable)) +
            geom_line() + geom_line())
    # Calculate the mean and standard deviation at each timepoint
    Summarised_vaccination_predictions <- melted_case_curves %>% 
      group_by(time) %>% 
      dplyr::summarise(vaccinations=mean(value/100), sd_vaccinations=sd(value/100))
    
    
    
    #Calculate the date of the cases
    Summarised_vaccination_predictions$date=as.Date(Summarised_vaccination_predictions$time + 
                                                        ORIGIN,
                                                      origin = "1970-01-01")
    
    return(Summarised_vaccination_predictions)
    
  
}



