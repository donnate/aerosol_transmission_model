#### Experiments
setwd("~/Dropbox/aerosol_transmission_model/")
source("code_experiments_paper/covid_projectons.R")
set.seed(12345)
B=100

dates = c(as.Date("2020-06-01"), 
          as.Date("2020-09-01"),
          as.Date("2020-11-01"),
          as.Date("2021-01-01"), 
          as.Date("2021-03-01"),
          as.Date("2021-07-01"))
countries = unique(COUNTRY_DATA$location)
res = data.frame(matrix(0,  length(countries) * length(dates), 42 ))
colnames(res)<-c("type", "date", "country", "period4fitting",
                 "ncurves",
                 "mean prev", "mean sd pred",
                 "MSE", "correlation",
                 "coverage", "bin bandwith", "max bandwidth",
                 "mean bandwdith", "sd bandwidth",
                 sapply(1:28, function(x){paste0("mean_predicted",x)})
)
COUNTRY_DATA <- read.csv(file="https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv", header=T, sep=",")
COUNTRY_DATA$date <- (as.Date(COUNTRY_DATA$date, "%Y-%m-%d"))
mu <- 3.96
sigma <- 4.75
cv <- sigma / mu   ##sqrt(alpha)/beta / (alpha/beta) = 1/sqrt(alpha)=cv   beta = 
params <- gamma_mucv2shapescale(mu, cv)
si <- distcrete("gamma", shape = params$shape,
                scale =  params$scale,   ### scale = 1/beta 
                interval = 1, w = 0)

i = 1
for (per in c(7,14, 21, 28)){
for (date in as.character(dates)){
  for (country in countries){
    print(c(country, date))
    ##### predict according to the different methods.
    tryCatch({
      res_1 <- predict_prevalence(origin=as.Date(date), country=country,
                                  country_data = COUNTRY_DATA,nb_curves=NCURVES, 
                                  distance=correlation_function, 
                                  period4predicting=28, period4fitting =per,
                                  distance_type = "correlation")
      if(is.null(res_1) == FALSE){
        temp_res = res_1$res%>% filter(time>0)
        res[i,] = c("correlation", date, country, 
                    mean(temp_res$prevalence), 
                    mean(temp_res$sd_prevalence), 
                    mean((temp_res$prevalence - temp_res$Observed )^2),
                    cor(temp_res$prevalence, temp_res$Observed),
                    mean(sapply(1:nrow(temp_res), function(i){ (temp_res$Observed[i] > temp_res$q25[i]) *(temp_res$Observed[i] < temp_res$q975[i]) })),
                    min(temp_res$q975 - temp_res$q25),
                    max(temp_res$q975 - temp_res$q25),
                    mean(temp_res$q975 - temp_res$q25),
                    sd(temp_res$q975 - temp_res$q25),
                    temp_res$prevalence)
                    
        i = i+1
        print(i)
      }
    }, error=function(e){cat("ERROR ",date, "\n")})
    
  }
  
  
}
}


colnames(res)<-c("type", "date", "country", sapply(1:28, function(x){paste0("mean_predicted",x)}),
  sapply(1:28, function(x){paste0("sd_predicted",x)}), "MSE", "correlation",
  "coverage", "CI band width"
)
