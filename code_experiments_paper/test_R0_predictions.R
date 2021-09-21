#### Experiments
#### Experiments
setwd("~/Dropbox/aerosol_transmission_model/")
library(outbreaks)
library(incidence)
library(ggplot2)
library(R0)
library(incidence)
library(earlyR)
library(projections)
set.seed(1)
library(distcrete)
library(epitrix)

source("code_experiments_paper/covid_projectons.R")
set.seed(12345)
B=100
mu <- 3.96
sigma <- 4.75
cv <- sigma / mu   ##sqrt(alpha)/beta / (alpha/beta) = 1/sqrt(alpha)=cv   beta = 
params <- gamma_mucv2shapescale(mu, cv)
si <- distcrete("gamma", shape = params$shape,
                scale =  params$scale,   ### scale = 1/beta 
                interval = 1, w = 0)
dates = c(as.Date("2020-06-01"), 
          as.Date("2020-09-01"),
          as.Date("2020-11-01"),
          as.Date("2021-01-01"), 
          as.Date("2021-03-01"),
          as.Date("2021-07-01"))
countries = unique(COUNTRY_DATA$location)
res = data.frame(matrix(0, 4* 4* length(dates) * length(countries), 127 ))
colnames(res)<-c("type", "date", "country", "period4fitting",
                 "ncurves",
                 "mean prev", "mean sd pred",
                 "RMSE","W-RMSE", "correlation",
                 "coverage", "bin bandwith", "max bandwidth",
                 "mean bandwdith", "sd bandwidth",
                 sapply(1:28, function(x){paste0("mean_predicted",x)}),
                 sapply(1:28, function(x){paste0("q25_predicted",x)}),
                 sapply(1:28, function(x){paste0("q975_predicted",x)}),
                 sapply(1:28, function(x){paste0("Observed",x)})
)
COUNTRY_DATA <- read.csv(file="https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv", header=T, sep=",")
COUNTRY_DATA$date <- (as.Date(COUNTRY_DATA$date, "%Y-%m-%d"))


names = c("correlation", "correlation_weighted", "MSE", "MSE_weighted", "R")
distance_type = c("correlation", "correlation", "MSE", "MSE")


i=1
for (it in 5:5){
  name = names[it]
    for (PERIOD_FOR_FITTING in c(7,14, 21, 28)){
      if (it%%2 == 1){
        w = NULL
      }else{
        w = (1:per)/sum(1:per)
      }
      print(c(PERIOD_FOR_FITTING))
      for (date in as.character(dates)){
        print(date)
        for (country in countries){
          ORIGIN = as.Date(date)
          ##### predict according to the different methods.
          tryCatch({
            chosen_location_data <- COUNTRY_DATA %>% 
              dplyr::filter((location==country) & (date >=  ORIGIN - PERIOD_FOR_FITTING  ) &  (date < ORIGIN)) %>%
              dplyr::select(new_cases_smoothed, date, population)
            ii <- unlist(sapply(1:PERIOD_FOR_FITTING, function(x){
              rep(chosen_location_data$date[x], max(0,round(chosen_location_data$new_cases_smoothed[x])))
            }))
            inc <- incidence(ii)
            
            obs <- COUNTRY_DATA %>% 
              dplyr::filter((location==country) & (date >  ORIGIN) &  (date <= ORIGIN + PERIOD_FOR_PREDICTING)) %>%
              dplyr::select(new_cases_smoothed_per_million)
            
            mGT<-generation.time("gamma", c(params$shape, params$scale))
            t <- estimate.R(chosen_location_data$new_cases_smoothed, mGT, begin=1, end=PERIOD_FOR_FITTING, methods=c("EG", "AR"),
                            pop.size=chosen_location_data$population[1], nsim=1000)
            
            
            pred_1 <- as.matrix(project(inc, R = t$estimates$EG$R, si = si, n_days = 28, n_sim = 1000))/ chosen_location_data$population[1]
            pred_2 <- as.matrix(project(inc, R = t$estimates$AR$R, si = si, n_days = 28, n_sim = 1000))/ chosen_location_data$population[1]

              q025 = apply(pred_1, 1, quantile, 0.025)
              q975 = apply(pred_1, 1, quantile, 0.975)
              res[i,] = c("EG", date, country, PERIOD_FOR_FITTING, NA,
                          1e6 * mean(apply(pred_1, 1, mean)), 
                          1e6 * mean(apply(pred_1, 1, sd)), 
                          1e6 * sqrt(mean((apply(pred_1, 1, mean) - obs$new_cases_smoothed_per_million/1e6 )^2)),
                          1e6 * sqrt(weighted.mean((apply(pred_1, 1, mean) - obs$new_cases_smoothed_per_million/1e6  )^2, 1:28/sum(1:28))),
                          cor(apply(pred_1, 1, mean), obs$new_cases_smoothed_per_million/1e6),
                          mean(sapply(1:28, function(ii){ (obs$new_cases_smoothed_per_million[ii]/1e6 > quantile(pred_1[ii,], 0.025)) * (obs$new_cases_smoothed_per_million[ii]/1e6 < quantile(pred_1[ii,], 0.975)) })),
                          1e6 * min(q975  - q025),
                          1e6 * max(q975  - q025),
                          1e6 * mean(q975  - q025),
                          1e6 * sd(q975  - q025),
                          apply(1e6 * pred_1, 1, mean),
                          apply(1e6 * pred_1, 1, quantile, 0.025),
                          apply(1e6 * pred_1, 1, quantile, 0.975),
                          obs$new_cases_smoothed_per_million)
              
              i = i+1
              
              q025 = apply(pred_2, 1, quantile, 0.025)
              q975 = apply(pred_2, 1, quantile, 0.975)
              res[i,] = c("AR", date, country, PERIOD_FOR_FITTING, NA,
                          1e6 *mean(apply(pred_2, 1, mean)), 
                          1e6 *mean(apply(pred_2, 1, sd)), 
                          1e6 *sqrt(mean((apply(pred_2, 1, mean) - obs$new_cases_smoothed_per_million/1e6 )^2)),
                          1e6 * sqrt(weighted.mean((apply(pred_2, 1, mean) - obs$new_cases_smoothed_per_million/1e6  )^2, 1:28/sum(1:28))),
                          cor(apply(pred_2, 1, mean), obs$new_cases_smoothed_per_million/1e6),
                          mean(sapply(1:28, function(ii){ (obs$new_cases_smoothed_per_million[ii]/1e6 > quantile(pred_2[ii,], 0.025)) * (obs$new_cases_smoothed_per_million[ii]/1e6 < quantile(pred_2[ii,], 0.975)) })),
                          1e6 *min(q975  - q025),
                          1e6 *max(q975  - q025),
                          1e6 *mean(q975  - q025),
                          1e6 *sd(q975  - q025),
                          apply(1e6 *pred_2, 1, mean),
                          apply(1e6 * pred_2, 1, quantile, 0.025),
                          apply(1e6 * pred_2, 1, quantile, 0.975),
                          obs$new_cases_smoothed_per_million)
              
              i = i+1
              
              print(i)
              write.csv(res,paste0("results_R0", ".csv"))
            
          }, error=function(e){cat("ERROR ",date, "\n")})
          
        }
        
      } 
      }
    
  }



