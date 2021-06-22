
ODDS_RATIO_VACCINE_DOSE1 <- c(rep(1,6),rep(0.9,3),rep(0.75,4),rep(0.47,7), rep(0.3,14), rep(0.25,14), 0.07 )[1:31]
ODDS_RATIO_VACCINE_DOSE2 <- c(rep(0.5, 6), rep(0.15, 7), rep(0.1,28), rep(0.07, 8))[1:31]
ODDS_RATIO_VACCINE_DOSE1_BEST <- c(rep(0.8,6),rep(0.75,3),rep(0.6,4),rep(0.37,7), rep(0.25,14), rep(0.18,14) )[1:31]
ODDS_RATIO_VACCINE_DOSE2_BEST  <- rep(0.05, 31) 
ODDS_RATIO_VACCINE_DOSE1_WORST <- c(rep(1,6),rep(1,3),rep(0.97,4),rep(0.6,7), rep(0.47,14), rep(0.4,14), 0.07 )[1:31]
ODDS_RATIO_VACCINE_DOSE2_WORST <- c(rep(0.7, 6), rep(0.47, 7), rep(0.33,28), rep(0.07, 8))[1:31]
# SD_ODDS_RATIO_VACCINE_DOSE1 <- c(rep(0.25,6),rep(0.25,3),rep(0.2,4),rep(0.1,7), rep(0.1,14), rep(0.1,14), 0.14)[1:31]
# ALPHA_ODDS_RATIO_VACCINE_DOSE2 <- c(rep(0.5, 6), rep(0.15, 7), rep(0.1,28), rep(18, 8)) [1:31]
# BETA_ODDS_RATIO_VACCINE_DOSE2 =  c(rep(0.5, 6), rep(0.15, 7), rep(0.1,28), rep(2, 8)) [1:31]
VACCINATIONS <- read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations.csv")
VACCINATIONS_US <- read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/us_state_vaccinations.csv")
intersection_col = intersect(colnames(VACCINATIONS), colnames(VACCINATIONS_US))
VACCINATIONS <- rbind(VACCINATIONS[,intersection_col], VACCINATIONS_US[,intersection_col] %>% filter(location != "United States"))


immunity<-function(country, date_event, origin, vaccination_data, odds1 =ODDS_RATIO_VACCINE_DOSE1,
                   odds2 = ODDS_RATIO_VACCINE_DOSE2, plot=FALSE){
  
  vac_df <-  vaccination_data %>% filter(location == country, date > as.Date("2021-01-15"), date < origin)
  vac_df$people_fully_vaccinated_per_hundred[is.na(vac_df$people_fully_vaccinated_per_hundred)] =0
  vac_df$people_vaccinated_per_hundred[is.na(vac_df$people_vaccinated_per_hundred)] =0
  out1 <- lm(people_vaccinated_per_hundred ~ date,  data = vac_df)
  out2 <- lm(people_fully_vaccinated_per_hundred ~ date, data = vac_df)
  
  if (plot){
    vac_pred <- data.frame( date= seq(from=origin -PERIOD_FOR_FITTING, to = date_event, by="day"))
    vac_pred["people_vaccinated_per_hundred"] = predict(out1, newdata=vac_pred)
    vac_pred["people_fully_vaccinated_per_hundred"] = predict(out2, newdata=vac_pred)
    vac_pred = rbind(vac_pred, vac_df[c("date", "people_vaccinated_per_hundred", "people_fully_vaccinated_per_hundred")])
    temps = vac_df[c("date", "people_vaccinated_per_hundred", "people_fully_vaccinated_per_hundred")]
    colnames(temps) = c("date", "Observed_vaccinated_per_hundred", "Observed_fully_vaccinated_per_hundred")
    vac_pred = merge(vac_pred, temps, on="date",
                     all.x=TRUE)
    ggplot(vac_pred, aes(x=date)) +
         geom_line(aes(y=people_vaccinated_per_hundred, colour="Predicted",linetype="Vaccinated"),size=1.4)+
         geom_line(aes(y=Observed_vaccinated_per_hundred, colour="Observed", linetype="Vaccinated"),size=1.4)+
         geom_line(aes(y=people_fully_vaccinated_per_hundred, colour="Predicted", linetype="Fully Vaccinated"),size=1.4)+
         geom_line(aes(y=Observed_fully_vaccinated_per_hundred, colour="Observed", linetype="Fully Vaccinated"),size=1.4)+
         theme_bw()  + 
         scale_colour_manual(breaks = c("Predicted","Observed"),  values=c("black","red"))+ 
         theme(legend.text=element_text(size=16))+labs(colour="Line Type", linetype="Vaccination Status", fill="95% CI Range", size=14) + 
         theme(text = element_text(size=20),
               axis.text.x = element_text(angle=90, hjust=1)) + 
         ylab("Proportion of the Population Vaccinated (%)")
  }
  
  vac_pred <- data.frame( date= seq(from= date_event - length(odds1), to = date_event, by="day"))
  vac_pred["people_vaccinated_per_hundred"] = predict(out1, newdata=vac_pred)
  vac_pred["people_fully_vaccinated_per_hundred"] = predict(out2, newdata=vac_pred)
  
  vac_pred = vac_pred %>% mutate(one_dose = people_vaccinated_per_hundred - people_fully_vaccinated_per_hundred )
  vac_pred = vac_pred %>% mutate(delta_one_dose = people_vaccinated_per_hundred - lag(people_vaccinated_per_hundred),
                                   delta_second_dose = people_fully_vaccinated_per_hundred - lag(people_fully_vaccinated_per_hundred))
    #### People are assumed to receive their  28 days after the first
  vac_pred[1,"delta_one_dose"] = vac_pred[1,"one_dose"] 
  vac_pred[1,"delta_second_dose"] = vac_pred[1,"people_fully_vaccinated_per_hundred"]
    
  
  
  vac_pred$delta_one_dose[is.na(vac_pred$delta_one_dose)] = 0
  vac_pred$delta_second_dose[is.na(vac_pred$delta_second_dose)] = 0
  immunity = 0.01 * sum(vac_pred["delta_second_dose"] *rev(1-odds2), na.rm=TRUE) +
      0.01 * sum(vac_pred["delta_one_dose"] *rev(1-odds1), na.rm=TRUE) 
  return(immunity)
}
  
  


   