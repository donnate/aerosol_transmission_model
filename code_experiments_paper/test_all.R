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
res = data.frame(matrix(0, 4 * 4 * length(countries) * length(dates), 42 ))
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

countries1 = c("Mali", "France", "Germany", "United States", "Vatican",
               "Zimbabwe", "Sudan","South Africa","Poland" )
names = c("correlation", "correlation_weighted", "MSE", "MSE_weighted" )
distance_type = c("correlation", "correlation", "MSE", "MSE")


i=0
for (it in 1:1){
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
                                  distance=correlation_function, 
                                  period4predicting=28, period4fitting =per,
                                  distance_type = distance_type[it],
                                  weights= w)
      if(is.null(res_1) == FALSE){
        temp_res = res_1$res%>% filter(time>0)
        res[i,] = c(names[it], date, country, per, NCURVES,
                    mean(temp_res$prevalence), 
                    mean(temp_res$sd_prevalence), 
                    mean((temp_res$prevalence - temp_res$Observed )^2),
                    cor(temp_res$prevalence, temp_res$Observed),
                    mean(sapply(1:nrow(temp_res), function(ii){ (temp_res$Observed[ii] > temp_res$q25[ii]) *(temp_res$Observed[ii] < temp_res$q975[ii]) })),
                    min(temp_res$q975 - temp_res$q25),
                    max(temp_res$q975 - temp_res$q25),
                    mean(temp_res$q975 - temp_res$q25),
                    sd(temp_res$q975 - temp_res$q25),
                    temp_res$prevalence)
        
        i = i+1
        print(i)
        write.csv(res,paste0("results", names[it], ".csv"))
      }
    }, error=function(e){cat("ERROR ",date, "\n")})
    
  }
  
  
}
}
}
}


