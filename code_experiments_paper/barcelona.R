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
input = list("country" =  "Spain",
             "date_event" = as.Date("2020-12-12"),
             "temperature" = 23,
             "UV"=0,
             "RH"=50, 
             "time2event"=1,
             "duration" = 3*60/60, #### time must be in hours
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
             "mask"=100.0,
             "prop_mask"=100,
             "mixing"=0,
             "p_symptoms" =0.4,
             "p_lie" = 0.5,
             file1=NULL 
)
B=50000
N_TOT = 472
surface =  1024
volume =5  *surface
country = input$country
PERIOD_FOR_FITTING = 14 
NCURVES = 100
MAX_DATE = input$date_event - 28
COUNTRY_DATA <- read.csv(file="https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv", header=T, sep=",")
COUNTRY_DATA$date <- (as.Date(COUNTRY_DATA$date, "%Y-%m-%d"))
VACCINATIONS <- read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations.csv")
VACCINATIONS_US <- read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/us_state_vaccinations.csv")
intersection_col = intersect(colnames(VACCINATIONS), colnames(VACCINATIONS_US))
VACCINATIONS <- rbind(VACCINATIONS[,intersection_col], VACCINATIONS_US[,intersection_col] %>% filter(location != "United States"))



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
######################################################################
######################################################################
##### Step 1: compute future prevalence
######################################################################
######################################################################

data2fit = COUNTRY_DATA %>% 
  dplyr::filter(location==input$country, date >  input$date_event- PERIOD_FOR_FITTING, date <= input$date_event) 

res_1 <- predict_prevalence(origin=input$date_event -PERIOD_FOR_PREDICTING, country=country,
                            country_data = COUNTRY_DATA,nb_curves=NCURVES, 
                            distance=Difference_function, 
                            period4predicting=PERIOD_FOR_PREDICTING, period4fitting = PERIOD_FOR_FITTING,
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

ggplot(future_prevalence_df)+
  geom_line(aes(x=Date_of_cases, y=1e6 * prevalence, colour="Predicted"))+
  geom_line(aes(x=Date_of_cases, y=1e6*Observed, colour="Observed"))+
  geom_ribbon(aes(x=Date_of_cases, ymin=1e6 * q25, ymax=1e6*q975), colour="grey", alpha=0.5)+
  theme_bw()+ xlab("Date") + ylab("Incidence (per Million)")+
  ggtitle("Predicted Incidence (per Million) Per Region")
future_prevalence_df["Date_event"] = input$date_event





##### If the event is less than 3 weeks into the future, some people mig

group_assignment = c(sapply(1:PERIOD_FOR_PREDICTING,
                            function(x){paste0(x)}))
df = pivot_wider(future_prevalence_df %>% dplyr::filter(time>=0) %>%
                   dplyr::select(time, prevalence), names_from = c("time"), values_from = "prevalence")
df = pivot_wider(future_prevalence_df %>% dplyr::filter(time>=0) %>%
                   dplyr::select(time, Observed), names_from = c("time"), values_from = "Observed")
proba_null <- future_prevalence_df[PERIOD_FOR_PREDICTING,"prevalence"]

#### Effect of the test + symptoms screening
nb_people_infected <- (N_TOT * apply(df,1,sum))
nb_people_infected_all <- rpois(B,N_TOT * apply(df,1,sum))

infectiousness_all = compute_relative_infectiousness(input, PERIOD_FOR_PREDICTING,
                                                     plot=FALSE)
library(data.table)
df_sample = (sapply(1:B, FUN=function(x){ 
  return((test %>% dplyr::filter(trajectory==  sample(1:NCURVES, 1) , time>=0))$value/1e6) * infectiousness_all$infectiousness}))
df_sample <- data.frame(matrix(rep( as.numeric(df), B), 
                               nrow=B, byrow=TRUE))

colnames(df_sample) = group_assignment
mult = sapply(1:29, function(i){ rnorm(n =B,infectiousness_all$infectiousness_event[i], infectiousness_all$infectiousness_event_sd[i])})
mult[mult<=0] =0
df_sample = df_sample * mult * 0.01

#### Sample  a tajectory
#for (i in 1:PERIOD_FOR_PREDICTING){
#  ii = PERIOD_FOR_PREDICTING - i
#  ind = which(infectiousness_all$Date.of.Infection == -(PERIOD_FOR_PREDICTING - i))
#  df_sample[as.character(i)] = rep(future_prevalence_df$prevalence[which(future_prevalence_df$time == i)], B)
#df_sample[as.character(i)]= sapply(1:B, function(b){
#       min(1,max(0,rnorm(1, future_prevalence_df$prevalence[which(future_prevalence_df$time == i)], 
#                                  future_prevalence_df$sd_prevalence[which(future_prevalence_df$time == i)]))) * max(min(1, 1e-2*rnorm(1, infectiousness_all$infectiousness_event[ind], infectiousness_all$infectiousness_event_sd[ind])),0)
#     })
#}


#for (i in 1:PERIOD_FOR_PREDICTING){
#  ii = PERIOD_FOR_PREDICTING - i
#  ind = which(infectiousness_all$Date.of.Infection == -(PERIOD_FOR_PREDICTING - i))
#  df_sample[as.character(i)] = rep(future_prevalence_df$prevalence[which(future_prevalence_df$time == i)], B)
#df_sample[as.character(i)]= sapply(1:B, function(b){
#       min(1,max(0,rnorm(1, future_prevalence_df$prevalence[which(future_prevalence_df$time == i)], 
#                                  future_prevalence_df$sd_prevalence[which(future_prevalence_df$time == i)]))) * max(min(1, 1e-2*rnorm(1, infectiousness_all$infectiousness_event[ind], infectiousness_all$infectiousness_event_sd[ind])),0)
#     })
#}

nb_people_infectious_at_the_event_all <- rpois(B, N_TOT* apply(df_sample[group_assignment ],1, sum))
nb_people_infectious_at_the_event <- mean(nb_people_infectious_at_the_event_all)
nb_people_detected <- nb_people_infected_all - nb_people_infectious_at_the_event_all 

#beta_immunity_params = beta.parms.from.quantiles(c(p_immunity_worst, p_immunity_best ))
# n_susc  <- (sapply(1:B, function(b){
#    #### ASSUME FIRST DOSE BECOMES FULLY VACCINATED AFTER FOUR WEEKS
#    rpois(1,sum(1-rbeta(N_TOT, beta_immunity_params$a, beta_immunity_params$b)))
#  }))



n_susc = N_TOT - nb_people_infectious_at_the_event_all
######### COMPUTE ASSOCIATED DISTRIBUTIONS
# dist_infectious<- bind_rows(lapply(1:B, function(b){
#   t = data.frame(sapply(REGIONS, function(region){sapply(1:100, function(y){dpois(y, sum(N[[region]] * df[region,group_assignment ] ))})}))
#   t["n"] = 1:100
#   return(t)
# }))
######### 






##########################################
##########################################
###### Step 2: TRANSMISSION DYNAMICS
##########################################
##########################################

##########################################
###### Step 2.a: compute room parameters for aerosolization
##########################################




occupant_density = N_TOT /surface
ventilation <- lookup_ventilation(input$ventilation, N_TOT, occupant_density, surface, volume)
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
nb_infective_people = apply(df_sample, 1, function(x){rpois(1,N_TOT * sum(x ))})

# proba_dist_infectiousness <- sapply(1:min(100, N_TOT), function(n){dpois(n,  
#                                                                          N_TOT * mean(apply(df_sample,1,sum)) )})
# 
# proba_infection <- sapply(1:min(100, N_TOT), function(n){
#   q_e = n * as.numeric(quanta_emission_rate0)
#   q_c <- as.numeric(compute_quanta_concentation(q_e, first_order_loss_rate = first_order_loss_rate, 
#                                                 volume =volume,
#                                                 duration =input$duration))
#   
#   quanta_inhaled_per_person <- compute_quanta_inhaled_per_person(quanta_concentration = q_c,  
#                                                                  breathing_rate = BREATHING_RATE,
#                                                                  duration=input$duration, 
#                                                                  inhalation_mask_efficiency = MASK_INHALATION_EFFICIENCY,
#                                                                  prop_mask = 0.01 * input$prop_mask)
#   #print(c(q_e, q_c,quanta_inhaled_per_person))
#   #print(paste0("quanta_inhaled_per_person: ",quanta_inhaled_per_person))
#   ###### Now there is a lot of uncertainty around that parameter
#   ###### We model potential super spreader events by using a pareto distirbution
#   return(quanta_inhaled_per_person)
# })
library(EnvStats)
alpha = log(5)/log(4)
proba_infection_sim <- sapply(nb_people_infectious_at_the_event_all, function(n){
  if (n>0){
    q_e = sum(rpareto(n, as.numeric(quanta_emission_rate0) /2^(1/alpha),   alpha))
    q_c <- as.numeric(compute_quanta_concentation(q_e, first_order_loss_rate = first_order_loss_rate, 
                                                  volume =volume,
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


# L  = min(max(which(proba_dist_infectiousness > 1e-10)), N_TOT)
# n_infections <- t(sapply(1:B, function(b){
#   sapply(1:L, function(n){
#     min(n_susc[b], rpareto(1, (n_susc[b] - n) * proba_infection[n] *0.16/1.16, 1.16))
#   })
# }))



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
proba_baseline <- df_sample$`28`

quantiles= c(0.99, 0.975, 0.5, 0.025, 0.01)
# infections_baseline <- data.frame(rbindlist(lapply(1:B, function(b){
#   a = data.frame(t(c(sapply(quantiles, function(x){qbinom(x, n_susc[b], 
#                                                           proba_baseline[b])}),
#                      proba_baseline[b] *  n_susc[b])))
#   colnames(a)= c("Q99.Number.of.Transmissions","Q97.5.Number.of.Transmissions", "Q50.Number.of.Transmissions", "Q2.5.Number.of.Transmissions" , "Q1.Number.of.Transmissions" , "Average.Number.of.Transmissions")
#   return(a)
# })))

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

print(res2)
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
vac_pred <- data.frame( date= seq(from= input$date_event -length(odds1), to = input$date_event, by="day"))
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

