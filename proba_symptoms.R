library(tidyverse)

##### Functions to compute the probability of having symptoms

### Data on incubation length: https://bmjopen.bmj.com/content/10/8/e039652
### Data on symptom length: https://www.nature.com/articles/s41467-020-20568-4
### The data associated to the latter is loaded through the following command:
data_sympt <- readxl::read_xlsx("41467_2020_20568_MOESM4_ESM.xlsx")
data_sympt <- data_sympt %>% filter(`virus culture result` == "POS")


# compute_proba_symptoms <- function(lag=2){
#   
#   PROBA_SYMPTOMS = apply(sapply(1:10000,function(b){
#     mu = rnorm(1,1.63, 0.12); 
#     s = rnorm(1,0.5, 0.05); 
#     ##### sample incubation rate
#     incubation = round(rlnorm(1,mu,s)) 
#     sympt = sample(data_sympt$`duration of symptoms in days`,1)
#     return(c(rep(0, incubation), rep(1, sympt), rep(0, 60- incubation - sympt)))
#   }),1,mean)
#   
#   
#   
#   PROBA_STILL_SYMPTOMS_IN_2_DAYS = apply(sapply(1:10000,function(b){
#     mu = rnorm(1,1.63, 0.12); 
#     s = rnorm(1,0.5, 0.05); 
#     ##### sample incubation rate
#     incubation = ceiling(rlnorm(1,mu,s)) 
#     sympt = sample(data_sympt$`duration of symptoms in days`,1)
#     x = c(rep(0, incubation), rep(1, sympt), rep(0, 60-incubation - sympt))
#     return(sapply(1:60, function(i){(x[i] * x[i+lag])}))
#   }),1,mean)
#   
#   PROBA_BUT_SYMPTOMS_IN_2_DAYS = apply(sapply(1:10000,function(b){
#     mu = rnorm(1,1.63, 0.12); 
#     s = rnorm(1,0.5, 0.05); 
#     ##### sample incubation rate
#     incubation = ceiling(rlnorm(1,mu,s)) 
#     sympt = sample(data_sympt$`duration of symptoms in days`,1)
#     x = c(rep(0, incubation), rep(1, sympt), rep(0, 60-incubation - sympt))
#     return(sapply(1:60, function(i){((1-x[i]) * x[i+lag])}))
#   }),1,mean)
#   
#   PROBA_AND_ASYMPTOM_IN_2_DAYS = apply(sapply(1:10000,function(b){
#     mu = rnorm(1,1.63, 0.12); 
#     s = rnorm(1,0.5, 0.05); 
#     ##### sample incubation rate
#     incubation = ceiling(rlnorm(1,mu,s)) 
#     sympt = sample(data_sympt$`duration of symptoms in days`,1)
#     x = c(rep(0, incubation), rep(1, sympt), rep(0, 60-incubation - sympt))
#     return(sapply(1:60, function(i){((1-x[i]) * (1-x[i+lag]))}))
#   }),1,mean)
#   
#   return(data.frame(rbind(PROBA_SYMPTOMS, PROBA_STILL_SYMPTOMS_IN_2_DAYS, PROBA_BUT_SYMPTOMS_IN_2_DAYS, PROBA_AND_ASYMPTOM_IN_2_DAYS), 
#                     row.names = c("proba_infection", "proba_sympt_and_still_in_d_days", "proba_asympt_but_sympt_in_d_days",
#                                   "proba_asympt_and_asympt_in_d_days")))
#   
# }

compute_proba_symptoms <- function(lag=2){
  uu = sapply(1:10000,function(b){
    mu = rnorm(1,1.63, 0.12); 
    s = rnorm(1,0.5, 0.05); 
    ##### sample incubation rate
    incubation = ceiling(rlnorm(1,mu,s)) 
    sympt = sample(data_sympt$`duration of symptoms in days`,1)
    return(c(rep(0, incubation), rep(1, sympt), rep(0, max(0,60- incubation - sympt))))
  })
  PROBA_SYMPTOMS = apply(uu ,1,mean)
  return(PROBA_SYMPTOMS)
}



