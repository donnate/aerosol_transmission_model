#### Experiments
setwd("~/Dropbox/aerosol_transmission_model/")
source("covid_projectons.R")
set.seed(12345)
B=100
NCURVES = 100
dates = seq(from=as.Date("2020-06-01"),to=as.Date("2021-02-10"), by=7)
countries = unique(COUNTRY_DATA$location)
res = data.frame(matrix(0, length(countries) * length(dates), 3 + 56 +  4  ))
COUNTRY_DATA <- read.csv(file="https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv", header=T, sep=",")
COUNTRY_DATA$date <- (as.Date(COUNTRY_DATA$date, "%Y-%m-%d"))
i = 1
for (date in as.character(dates)){
  for (country in countries){
    print(c(country, date))
    ##### predict according to the different methods.
    tryCatch({
      res_1 <- predict_prevalence(origin=as.Date(date), country=country,
                                  country_data = COUNTRY_DATA,nb_curves=NCURVES, 
                                  distance=correlation_function, 
                                  period4predicting=28, period4fitting =28,
                                  distance_type = "correlation")
      if(is.null(res_1) == FALSE){
        res[i,] = c("correlation", date, country, res_1$prevalence[29:56], 
                    res_1$sd_prevalence[29:56], 
                    mean((res_1$prevalence[29:56] - res_1$Observed[29:56] )^2), cor(res_1$prevalence[29:56],res_1$Observed[29:56]),
                    coverage = mean(sapply(29:56, function(i){ (res_1$Observed[i] > res_1$prevalence[i] -2 * res_1$sd_prevalence[i]) *(res_1$Observed[i] < res_1$prevalence[i]  + 2 * res_1$sd_prevalence[i]) 
                    })), 4 * mean(res_1$sd_prevalence[29:56]))
        i = i+1
        print(i)
      }
    }, error=function(e){cat("ERROR ",date, "\n")})
    
    
  }
  
  
}


colnames(res)<-c("type", "date", "country", sapply(1:28, function(x){paste0("mean_predicted",x)}),
  sapply(1:28, function(x){paste0("sd_predicted",x)}), "MSE", "correlation",
  "coverage", "CI band width"
)
