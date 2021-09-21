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
             "date_event" = as.Date("2021-08-20"),
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
B=50000
N_TOT = 5000
volume = 86650
surface = 20000

country = input$country
PERIOD_FOR_FITTING = 14 
NCURVES = 100
MAX_DATE = input$date_event - 28
COUNTRY_DATA <- read.csv(file="https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv", header=T, sep=",")
COUNTRY_DATA$date <- (as.Date(COUNTRY_DATA$date, "%Y-%m-%d"))



PERIOD_FOR_PREDICTING = 28
DECAY = max(0, (7.57+ 
                  1.41* (input$temperature-20.54)/10.66 +
                  0.0218 *(input$RH-45.24)/28.67 + 
                  7.55 *((input$UV*0.185)-50) / 50 +
                  (input$temperature-20.54)/10.66*(input$UV*0.185-50)/50*1.40) *60)  #https://www.dhs.gov/science-and-technology/sars-airborne-calculator

BREATHING_RATE = 0.012 * 60 
DEPOSITION = 0.24
MASK_INHALATION_EFFICIENCY = 0.5
PRESSURE = 0.95


data2fit = COUNTRY_DATA %>% 
  dplyr::filter(location=="United Kingdom", date >  input$date_event- PERIOD_FOR_FITTING, date <= input$date_event) 

res_1 <- predict_prevalence(origin=input$date_event -28, country=country,
                            country_data = COUNTRY_DATA,nb_curves=NCURVES, 
                            distance=Difference_function, 
                            period4predicting=28, period4fitting = PERIOD_FOR_FITTING,
                            distance_type = "MSE", agg=median)

res_1$res["Date"] = input$date_event
res_1$output["Date"] = input$date_event


prevalence_df = res_1$res
test = res_1$output
test["trajectory"] = sapply(1:nrow(test), FUN=function(x){x%/%(PERIOD_FOR_FITTING+PERIOD_FOR_PREDICTING)})
prevalence_df["region"] = country
samples_prev = res_1$samples
samples_prev["region"] = country

bias<- compute_underascertainment_bias(input$date_event - (PERIOD_FOR_FITTING+PERIOD_FOR_PREDICTING), "United Kingdom", 
                                       COUNTRY_DATA, date_max = input$date_event - 28, 
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


group_assignment = c(sapply(1:PERIOD_FOR_PREDICTING,
                            function(x){paste0(x)}))
df = pivot_wider(future_prevalence_df %>% dplyr::select(time, prevalence), names_from = c("time"), values_from = "prevalence")



nb_people_infected <- (N_TOT * apply(df,1,sum))
nb_people_infected_all <- (N_TOT * apply(df,1,sum))
infectiousness_all = compute_relative_infectiousness(input, PERIOD_FOR_PREDICTING,
                                                     plot=FALSE)
library(data.table)
df_sample = (sapply(1:B, FUN=function(x){ 
  return((test %>% dplyr::filter(trajectory==  sample(1:NCURVES, 1) , time>0))$value/1e6)}))
df_sample <- data.frame(matrix(unlist(df_sample), 
                               nrow=B, byrow=TRUE))
colnames(df_sample) = group_assignment







##########################################
##########################################
###### Step 2: TRANSMISSION DYNAMICS
##########################################
##########################################

##########################################
###### Step 2.a: compute room parameters for aerosolization
##########################################


res_tot = c()
for (N_TOT in c(10, 20, 50, 100, 200, 300, 500, 800, 1000, 2000, 3000, 4000, 5000, 7000, 10000 )){
  volume = 86650/ 5000 * N_TOT
  surface = 20000/5000 * N_TOT
  print(N_TOT)
  nb_people_infected <- (N_TOT * apply(df,1,sum))
  nb_people_infected_all <- (N_TOT * apply(df,1,sum))
  
  
  region="United Kingdom"
  
  #nb_people_infectious_at_the_event_all <- rpois(B, N_TOT* apply(df_sample[group_assignment ],1, sum))
  #nb_people_infectious_at_the_event <- mean(N_TOT* apply(df_sample[group_assignment ],1, sum))
  n_susc = rep(N_TOT- max(1, 0.1*N_TOT), B)
  nb_people_vulnerable_at_the_event <- mean(n_susc)
  ######### 
  
  occupant_density = N_TOT /surface 
  ventilation <- lookup_ventilation(input$ventilation, N_TOT, occupant_density, 
                                    surface, volume)
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
  
  
  ##### Maybe I should add a negative binomial
  nb_infective_people = rep( max(1, 0.1*N_TOT), B)
  
  library(EnvStats)
  alpha = log(5)/log(4)
  quanta_emission_rate0 <- compute_quanta_emission_rate(input$activity,
                                                        MASK_EFFICIENCY ,
                                                        input$prop_mask/100,
                                                        1)
  
  proba_infection_sim <- sapply(nb_infective_people, function(n){
    if (n>0){
      q_e = sum(rpareto(n, as.numeric(quanta_emission_rate0) /2^(1/alpha),   alpha))
      q_c <- as.numeric(compute_quanta_concentation(q_e, first_order_loss_rate = first_order_loss_rate, 
                                                    volume =volume  ,
                                                    duration =input$duration))
      
      quanta_inhaled_per_person <- compute_quanta_inhaled_per_person(quanta_concentration = q_c,  
                                                                     breathing_rate = BREATHING_RATE,
                                                                     duration=input$duration, 
                                                                     inhalation_mask_efficiency = MASK_INHALATION_EFFICIENCY,
                                                                     prop_mask = 0.01 * input$prop_mask)
      #print(c(q_e, q_c,quanta_inhaled_per_person))
      #print(paste0("quanta_inhaled_per_person: ",quanta_inhaled_per_person))
      ###### Now there is a lot of uncertainty around that parameter
      ###### We model potential super spreader events by using a pareto distirbution
      return(1-exp(-quanta_inhaled_per_person))
    }else{
      return(0)
    }
  })
  
  n_infections_sim <- sapply(1:B, function(b){
    sum(sapply(proba_infection_sim[b],
               FUN=function(x){rbinom(n_susc[b],1,x)}))
  })
  
  # ggplot(data = data.frame(x=as.factor(cut(n_tot[1:B], 
  #                                          breaks=c(1,20, 50, 100, 200, 1000))),
  #                          y= 100 * n_infections_sim[1:B]/ (n_susc[1:B]),
  #                          pop=n_tot[1:B],
  #                          density = n_density[1:B]), aes(x=x, y=y))+
  #   theme_bw()+
  #   geom_boxplot(outlier.shape=NA)+
  #   geom_jitter(aes(colour="Simulated"), shape=4,alpha=0.7,width = 0.25) + 
  #   geom_point(data=data.frame(x = cut(data_with_index$Total.Pop.at.Event, breaks=c(1,20, 50, 100, 200, 1000)), 
  #                              y = 100 * data_with_index$Secondary.Cases /(data_with_index$Total.Pop.at.Event - data_with_index$Index.Cases),
  #                              density = data_with_index$Setting),
  #              aes(x=x, y=y,colour="Observed"), shape=17, size=2)+
  #   geom_text(data=data.frame(x = cut(data_with_index$Total.Pop.at.Event, breaks=c(1,20, 50, 100, 200, 1000)), 
  #                              y = 100 * data_with_index$Secondary.Cases /(data_with_index$Total.Pop.at.Event - data_with_index$Index.Cases),
  #                              label = data_with_index$Setting1,
  #                             density = data_with_index$Setting),
  #              aes(x=x, y=y, label=label), size=3, hjust=0.5, vjust=1)+
  #   scale_y_log10() + 
  #   xlab("Event Size") +
  #   ylab("Percentage of Secondary Infections") + 
  #   scale_color_manual(values = c( "indianred", "grey")) + 
  #   labs(colour="Colour")+
  #   theme(text = element_text(size=20),
  #         axis.text.x = element_text(angle=90, hjust=1, vjust=0.3))  +
  #   facet_wrap(.~density)
  
  
  
  # L  = min(max(which(proba_dist_infectiousness > 1e-10)), N_TOT)
  # n_infections <- t(sapply(1:B, function(b){
  #   sapply(1:L, function(n){
  #     min(n_susc[b], rpareto(1, (n_susc[b] - n) * proba_infection[n] *0.16/1.16, 1.16))
  #   })
  # }))
  res_tot = rbind(res_tot,
                  data.frame(x=n_infections_sim, "n" = N_TOT, date=input$date_event))
}

library(ggplot2)

sumup = res_tot %>% group_by(n,date) %>% summarise(m=mean((x+ 0.1*n)), med=median((x+ 0.1*n)/n),
                                            q25=quantile( (x+ 0.1*n)/n, 0.025),
                                           q95=quantile((x+ 0.1*n)/n, 0.95),
                                           q99=quantile((x+ 0.1*n)/n, 0.99),
                                           q999=quantile((x+ 0.1*n)/n, 0.999))
ggplot(data=sumup,
       aes(x=n, y=med)) +
  #geom_jitter(data=res_tot, aes(x=n, y=x/n), alpha=0.2, shape=4, size=0.2)+
  geom_smooth(aes(y=q99, colour="99th Quantile")) +
  #geom_line(aes(y=q999, colour="99.9th Quantile")) +
  geom_smooth(aes(y=q999, colour="99.9th Quantile"))+
  geom_smooth(aes(y=med, colour="Median")) + theme_bw() +ylab("Percentage of Total Infections")


res = data.frame("Average Number of Transmissions" = mean(n_infections_sim ),
                 "Q97.5 Number of Transmissions" = quantile((n_infections_sim ), 0.975),
                 "Q2.5 Number of Transmissions" = quantile((n_infections_sim ), 0.025),
                 "Q99 Number of Transmissions" = quantile((n_infections_sim ), 0.99),
                 "Q1 Number of Transmissions" = quantile((n_infections_sim ), 0.01),
                 "Q50 Number of Transmissions" = quantile((n_infections_sim ), 0.5),
                 "date"=input$date_event)

print(res)
#### Tail behaviour analysis:
n_tails = n_infections_sim[n_infections_sim>=56]
###
hist(n_tails/N_TOT)
res_march = rbind(res_march, res)

print(mean(sapply(1:B, function(b){ qbinom( 0.99, n_susc[b], as.numeric(proba_null))})))

#n_simulated_infections <-c()
n_simulated_infections <- rbind(n_simulated_infections, 
                                data.frame("n" = t(sapply(1:B, function(b){
                                  ifelse(nb_infective_people[b] == 0, 0, min(n_susc[b], rpareto(1, (n_susc[b] - nb_infective_people[b]) * proba_infection[nb_infective_people[b]] *0.16/1.16, 1.16)))
                                })), type="No Masks", date=input$date_event)
)


library(data.table)
library(reshape2)
library(tibbletime)
prev_day_event =  future_prevalence_df %>% filter(Date_of_cases == input$date_event)
proba_baseline <- t(sapply(1:B, function(b){
  max(0,rnorm(n=1, mean=prev_day_event$prevalence,
              sd=prev_day_event$sd_prevalence))
}))

quantiles= c(0.99, 0.975, 0.5, 0.025, 0.01)
infections_baseline <- data.frame(rbindlist(lapply(1:B, function(b){
  a = data.frame(t(c(sapply(quantiles, function(x){qbinom(x, n_susc[b], proba_baseline[b])}), proba_baseline[b] *  n_susc[b])))
  colnames(a)= c("Q99.Number.of.Transmissions","Q97.5.Number.of.Transmissions", "Q50.Number.of.Transmissions", "Q2.5.Number.of.Transmissions" , "Q1.Number.of.Transmissions" , "Average.Number.of.Transmissions")
  return(a)
})))

infections_baseline <- sapply(1:B, function(b){
  a = rbinom(1, n_susc[b], proba_baseline[b])})
colnames(a)= c("Q99.Number.of.Transmissions","Q97.5.Number.of.Transmissions", "Q50.Number.of.Transmissions", "Q2.5.Number.of.Transmissions" , "Q1.Number.of.Transmissions" , "Average.Number.of.Transmissions")
res2 = data.frame("Average Number of Transmissions" = mean(infections_baseline),
                  "Q97.5 Number of Transmissions" = quantile((infections_baseline), 0.975),
                  "Q2.5 Number of Transmissions" = quantile((infections_baseline), 0.025),
                  "Q99 Number of Transmissions" = quantile((infections_baseline ), 0.99),
                  "Q1 Number of Transmissions" = quantile((infections_baseline), 0.01),
                  "Q50 Number of Transmissions" = quantile((infections_baseline ), 0.5),
                  "date"=input$date_event)


infections_baseline_res = apply(infections_baseline , 2, mean)
res = rbind(res, apply(infections_baseline , 2, mean))


STOP
n_simulated_infections = data.frame("n"= n_simulated_infections,
                                    "type"= c(rep("2020-08-03", B),rep("2021-01-18", B),rep("2021-03-08", B)))
n_simulated_infections ["type"]= c(rep("2020-08-03", B),rep("2021-01-18", B),rep("2021-03-08", B))

library(cowplot)
n_simulated_infections2 =pivot_longer(n_simulated_infections, cols = sapply(1:5000, function(x){paste0("n.", x)}))
n_simulated_infections2$type = factor(n_simulated_infections2$type , levels=c("No Masks", "Half Mask", "All Masks", "Half Capacity", "Half Duration"))
ggplot(n_simulated_infections2, aes(x = as.factor(type),y=value, fill=as.factor(type)))+ theme_bw() +
  stat_boxplot(geom ='errorbar', width = 0.3) +
  geom_boxplot(outlier.colour = NA) + ylim(0, 50) + scale_y_log10() + 
  xlab("Date of Event") +  scale_fill_brewer(palette="Dark2") + 
  theme(legend.text=element_text(size=16))+labs(colour="Date", fill="type", size=14) + 
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=0, hjust=0.5)) + 
  ylab("Number of infections")  + xlab("")




odds1 =ODDS_RATIO_VACCINE_DOSE1
odds2 = ODDS_RATIO_VACCINE_DOSE2
vac_df <-  VACCINATIONS %>% filter(location == region, date < input$date_event -28)
out1 <- lm(people_vaccinated_per_hundred ~ date, data = vac_df)
out2 <- lm(people_fully_vaccinated_per_hundred ~ date, data = vac_df)
vac_pred <- data.frame( date= seq(from= date_event -length(odds1), to = date_event, by="day"))
vac_pred["people_vaccinated_per_hundred"] = predict(out1, newdata=vac_pred)
vac_pred["people_fully_vaccinated_per_hundred"] = predict(out2, newdata=vac_pred)

vac_pred = vac_pred %>% mutate(one_dose = people_vaccinated_per_hundred - people_fully_vaccinated_per_hundred )
vac_pred = vac_pred %>% mutate(delta_one_dose = people_vaccinated_per_hundred - lag(people_vaccinated_per_hundred),
                               delta_second_dose = people_fully_vaccinated_per_hundred - lag(people_fully_vaccinated_per_hundred))
#### People are assumed to receive their  28 days after the first
vac_pred[1,"delta_one_dose"] = as.numeric(vac_df %>% filter(date == max(vac_df$date)) %>% mutate(m=people_vaccinated_per_hundred - people_fully_vaccinated_per_hundred ) %>% select(m))
vac_pred[1,"delta_second_dose"] = as.numeric(vac_df %>% filter(date == max(vac_df$date)) %>% mutate(m=people_fully_vaccinated_per_hundred ) %>% select(m))


p_immunity = 0.01 * sum(vac_pred["delta_second_dose"] *rev(1-odds2)) +
  0.01 * sum(vac_pred["delta_one_dose"] *rev(1-odds1)) 
#vac_pred[1,"delta_second_dose"] = 0


odds1 =ODDS_RATIO_VACCINE_DOSE1_WORST
odds2 = ODDS_RATIO_VACCINE_DOSE2_WORST
p_immunity_worst = 0.01 * sum(vac_pred["delta_second_dose"] *rev(1-odds2)) +
  0.01 * sum(vac_pred["delta_one_dose"] *rev(1-odds1)) 

odds1 =ODDS_RATIO_VACCINE_DOSE1_BEST
odds2 = ODDS_RATIO_VACCINE_DOSE2_BEST
p_immunity_best = 0.01 * sum(vac_pred["delta_second_dose"] * rev(1-odds2)) +
  0.01 * sum(vac_pred["delta_one_dose"] *rev(1-odds1)) 

