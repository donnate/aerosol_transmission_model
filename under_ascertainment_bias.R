compute_underascertainment_bias <- function(origin,  country, country_data, date_max = NULL, Case_to_death_delay= 21, plot=FALSE){
    ifr_data=read.csv("country_ifr_data.csv", header=T)
    # merging the country case data with IFR data
    cases_and_ifr = merge(country_data, ifr_data, by="location")
    # calculating mean ifr
    cases_and_ifr$mean_ifr <- rowMeans(cases_and_ifr[,c("ENE_COVID","COVID_US_CDC","COVID_Verity","COVID_Levin")], na.rm=TRUE)
    # calculate the 'true cases' from ifr and death rate before time shifting
    cases_and_ifr$estimated_cases_per_million<-cases_and_ifr$new_deaths_smoothed_per_million/cases_and_ifr$mean_ifr * 100

        # shift the true cases by case_to_death_delay
    # first create new dataframe of estimated cases by shifted date
    estimated_cases=data.frame(location=cases_and_ifr$location, date=cases_and_ifr$date-Case_to_death_delay, estimated_cases_per_million=cases_and_ifr$estimated_cases_per_million)
    # remove estimated cases per million in cases_and_ifr dataset then merge with estimated_cases dataframe
    cases_and_ifr<-dplyr::select(cases_and_ifr,-c(estimated_cases_per_million)) %>%
      merge(.,estimated_cases, by=c("location", "date"))
    # Calculate the proportion of cases found
    cases_and_ifr$proportion_of_cases_detected=cases_and_ifr$new_cases_smoothed_per_million/cases_and_ifr$estimated_cases_per_million
    # Limit the maximum proportion of cases found to 1.0 and minimum to 0
    cases_and_ifr$proportion_of_cases_detected=ifelse(cases_and_ifr$proportion_of_cases_detected>1,1,ifelse(cases_and_ifr$proportion_of_cases_detected<0,0,cases_and_ifr$proportion_of_cases_detected))
    # specific countries example of ascertainment rate
    # select rows and columns of interest then pivot longer
    selected_countries_example = cases_and_ifr %>% filter(location %in% c(country)) %>%
      dplyr::select("location", "estimated_cases_per_million", "new_cases_smoothed_per_million", "new_deaths_smoothed_per_million", "date", "proportion_of_cases_detected") %>%
      pivot_longer(cols = c("estimated_cases_per_million", "new_cases_smoothed_per_million", "new_deaths_smoothed_per_million", "proportion_of_cases_detected"))
   
    if(plot){
     print(ggplot(selected_countries_example, aes(x=as.Date(date, origin = "1970-01-01"), y=value, colour=name))+
      geom_point(size=0.2)+
      geom_smooth(span=0.1)+
      scale_y_log10()+
      labs(y="value", x="date")+
      theme_bw()+
      theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)))
      
      print(ggplot(selected_countries_example %>% filter(name== "proportion_of_cases_detected"), aes(x=as.Date(date, origin = "1970-01-01"), y=value), colour="sierra")+
              geom_point(size=0.2)+
              geom_smooth(span=0.1)+
              labs(y="value", x="date")+
              theme_bw()+  ylab("Proportion of Cases Detected") + ylim(0,1) + 
              theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)))
    }
    if (is.null(date_max)){
      return(selected_countries_example %>% filter(date > origin, name == "proportion_of_cases_detected" ) %>% select(date, value) )
    }else{
      return(selected_countries_example %>% filter(date > origin, date<date_max, name == "proportion_of_cases_detected" ) %>% select(date, value) )
    }
    
}


compute_underascertainment_bias2 <- function(origin,  country, country_data, date_max = NULL, Case_to_death_delay= 21, plot=FALSE){
  ifr_data=read.csv("country_ifr_data.csv", header=T)
  # merging the country case data with IFR data
  cases_and_ifr = merge(country_data, ifr_data, by="location")
  # calculating mean ifr
  cases_and_ifr$mean_ifr <- rowMeans(cases_and_ifr[,c("ENE_COVID","COVID_US_CDC","COVID_Verity","COVID_Levin")], na.rm=TRUE)
  # calculate the 'true cases' from ifr and death rate before time shifting
  cases_and_ifr$estimated_cases_per_million<-cases_and_ifr$new_deaths_smoothed_per_million/cases_and_ifr$mean_ifr * 100
  
  # shift the true cases by case_to_death_delay
  # first create new dataframe of estimated cases by shifted date
  estimated_cases=data.frame(location=cases_and_ifr$location, date=cases_and_ifr$date-Case_to_death_delay, estimated_cases_per_million=cases_and_ifr$estimated_cases_per_million)
  # remove estimated cases per million in cases_and_ifr dataset then merge with estimated_cases dataframe
  cases_and_ifr<-dplyr::select(cases_and_ifr,-c(estimated_cases_per_million)) %>%
    merge(.,estimated_cases, by=c("location", "date"))
  # Calculate the proportion of cases found
  cases_and_ifr$proportion_of_cases_detected=cases_and_ifr$new_cases_smoothed_per_million/cases_and_ifr$estimated_cases_per_million
  # Limit the maximum proportion of cases found to 1.0 and minimum to 0
  cases_and_ifr$proportion_of_cases_detected=ifelse(cases_and_ifr$proportion_of_cases_detected>1,1,ifelse(cases_and_ifr$proportion_of_cases_detected<0,0,cases_and_ifr$proportion_of_cases_detected))
  # specific countries example of ascertainment rate
  # select rows and columns of interest then pivot longer
  selected_countries_example = cases_and_ifr %>% filter(location %in% c(country)) %>%
    dplyr::select("location", "estimated_cases_per_million", "new_cases_smoothed_per_million", "new_deaths_smoothed_per_million", "date", "proportion_of_cases_detected") %>%
    pivot_longer(cols = c("estimated_cases_per_million", "new_cases_smoothed_per_million", "new_deaths_smoothed_per_million", "proportion_of_cases_detected"))
  
  if(plot){
    print(ggplot(selected_countries_example, aes(x=as.Date(date, origin = "1970-01-01"), y=value, colour=name))+
            geom_point(size=0.2)+
            geom_smooth(span=0.1)+
            scale_y_log10()+
            labs(y="value", x="date")+
            theme_bw()+
            theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)))
    
    print(ggplot(selected_countries_example %>% filter(name== "proportion_of_cases_detected"), aes(x=as.Date(date, origin = "1970-01-01"), y=value), colour="sierra")+
            geom_point(size=0.2)+
            geom_smooth(span=0.1)+
            labs(y="value", x="date")+
            theme_bw()+  ylab("Proportion of Cases Detected") + ylim(0,1) + 
            theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)))
  }
  if (is.null(date_max)){
    return(selected_countries_example %>% filter(date > origin, name == "proportion_of_cases_detected" ) %>% select(date, value) )
  }else{
    return(selected_countries_example %>% filter(date > origin, date<date_max, name == "proportion_of_cases_detected" ) %>% select(date, value) )
  }
  
}

