#### Experiments
setwd("~/Dropbox/aerosol_transmission_model/")
source("code_experiments_paper/covid_projectons.R")
set.seed(12345)
B=100

args <- commandArgs(trailingOnly = TRUE)
it = as.numeric(args[1])
COUNTRY_DATA <- read.csv(file="https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv", header=T, sep=",")
COUNTRY_DATA$date <- (as.Date(COUNTRY_DATA$date, "%Y-%m-%d"))
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


names = c("correlation", "correlation_weighted", "MSE", "MSE_weighted" )
distance_type = c("correlation", "correlation", "MSE", "MSE")
dist = c(correlation_function, 
         Difference_function)

i=1

  name = names[it]
  for (NCURVES in c(100,200, 300, 500)){
    for (per in c(7,14, 21, 28)){
      if (it%%2 == 1){
        w = NULL
      }else{
        w = (1:per)/sum(1:per)
      }
      
      for (date in as.character(dates)){
        for (country in countries){
          print(c(country, date))
          ##### predict according to the different methods.
          tryCatch({
            res_1 <- predict_prevalence(origin=as.Date(date), country=country,
                                        country_data = COUNTRY_DATA,nb_curves=NCURVES, 
                                        distance=dist[[(it-1)%/%2 + 1]], 
                                        period4predicting=28, period4fitting =per,
                                        distance_type = distance_type[it],
                                        weights= w, agg=median)
            if(is.null(res_1) == FALSE){
              temp_res = res_1$res%>% filter(time>0)
              res[i,] = c(names[it], date, country, per, NCURVES,
                          1e6 * mean(temp_res$prevalence), 
                          1e6 * mean(temp_res$sd_prevalence), 
                          1e6 * sqrt(mean((temp_res$prevalence - temp_res$Observed )^2)),
                          1e6 * sqrt(weighted.mean((temp_res$prevalence - temp_res$Observed )^2, 1:28/sum(1:28))),
                          cor(temp_res$prevalence, temp_res$Observed),
                          mean(sapply(1:nrow(temp_res), function(ii){ (temp_res$Observed[ii] > temp_res$q25[ii]) *(temp_res$Observed[ii] < temp_res$q975[ii]) })),
                          1e6 * min(temp_res$q975 - temp_res$q25),
                          1e6 * max(temp_res$q975 - temp_res$q25),
                          1e6 * mean(temp_res$q975 - temp_res$q25),
                          1e6 * sd(temp_res$q975 - temp_res$q25),
                          1e6 * temp_res$prevalence,
                          1e6 * temp_res$q25,
                          1e6 * temp_res$q975,
                          1e6 * temp_res$Observed)
              
              i = i+1
              print(i)
              write.csv(res,paste0("results_median", names[it], ".csv"), row.names = FALSE)
            }
          }, error=function(e){
            cat("ERROR ",date, "\n")
            res[i,1:5] = c(names[it],
                           date, country, per, NCURVES)
            i = i +1
            })
          
        }
      }
    }
  }



