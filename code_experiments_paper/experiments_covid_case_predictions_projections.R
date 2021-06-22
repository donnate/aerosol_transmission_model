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

source("covid_case_predictions.R")
set.seed(12345)
B=100
NCURVES = 100
dates = seq(from=as.Date("2020-06-01"),to=as.Date("2021-02-10"), by=28)

COUNTRY_DATA <- read.csv(file="https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv", header=T, sep=",")
COUNTRY_DATA$date <- (as.Date(COUNTRY_DATA$date, "%Y-%m-%d"))
countries = unique(COUNTRY_DATA$location)
PERIOD_FOR_FITTING  = 28
PERIOD_FOR_PREDICTING  = 28
mu <- 3.96
sigma <- 4.75
cv <- sigma / mu   ##sqrt(alpha)/beta / (alpha/beta) = 1/sqrt(alpha)=cv   beta = 
params <- gamma_mucv2shapescale(mu, cv)
si <- distcrete("gamma", shape = params$shape,
                scale =  params$scale,   ### scale = 1/beta 
                interval = 1, w = 0)

res = list()
samples_knn= list()
samples_R = list()
a = 1
for (date in as.character(dates)){
  for (country in countries){
    print(c(country, date))
    ORIGIN = as.Date(date)
    ##### predict according to the different methods.
    tryCatch({
      
      chosen_location_data <- COUNTRY_DATA %>% 
        dplyr::filter((location==country) & (date >=  ORIGIN - PERIOD_FOR_FITTING  ) &  (date < ORIGIN)) %>%
        dplyr::select(new_cases_smoothed, date, population)
      ii <- unlist(sapply(1:PERIOD_FOR_FITTING, function(x){
        rep(chosen_location_data$date[x], max(0,round(chosen_location_data$new_cases_smoothed[x])))
      }))
      i <- incidence(ii)
      
      obs <- COUNTRY_DATA %>% 
        dplyr::filter((location==country) & (date >  ORIGIN) &  (date <= ORIGIN + PERIOD_FOR_PREDICTING)) %>%
        dplyr::select(new_cases_smoothed_per_million)
      
      mGT<-generation.time("gamma", c(params$shape, params$scale))
      t <- estimate.R(chosen_location_data$new_cases_smoothed, mGT, begin=1, end=28, methods=c("EG", "AR"),
                      pop.size=chosen_location_data$population[1], nsim=1000)
      
      
      pred_1 <- as.matrix(project(i, R = t$estimates$EG$R, si = si, n_days = 28, n_sim = 1000))/ chosen_location_data$population[1]
      pred_2 <- as.matrix(project(i, R = t$estimates$AR$R, si = si, n_days = 28, n_sim = 1000))/ chosen_location_data$population[1]
      if(TRUE){
        res = rbind(res, data.frame(date= as.Date(date),
                                    time = as.Date(date) + 1:28, "pred" = 1e6 * apply(pred_1, 1, mean),"sd" =  1e6 *  apply(pred_1, 1,sd), "type"="EG", "obs"= obs$new_cases_smoothed_per_million))
        res = rbind(res, data.frame(date= as.Date(date),
                                    time = as.Date(date) + 1:28, "pred" = 1e6 * apply(pred_2, 1, mean),"sd" =  1e6 *  apply(pred_1, 1,sd), "type"="AR", "obs"= obs$new_cases_smoothed_per_million))
        samples_R <- rbind(samples_R, data.frame(date= as.Date(date),
                                                 time = as.Date(date) + 1:28,  pred_1, "type"="EG", "obs"= obs$new_cases_smoothed_per_million))
        samples_R <- rbind(samples_R, data.frame(date= as.Date(date),
                                         time = as.Date(date) + 1:28,  pred_2, "type"="AR", "obs"= obs$new_cases_smoothed_per_million))
        print(a)
      }
    }, error=function(e){cat("ERROR ",date, "\n")})
    
    tryCatch({
      chosen_location_data <- COUNTRY_DATA %>% 
        dplyr::filter((location==country) & (date >=  ORIGIN - PERIOD_FOR_FITTING  ) &  (date < ORIGIN)) %>%
        dplyr::select(new_cases_smoothed, date, population)
      res_1 <- predict_prevalence(origin=as.Date(date), country=country,
                                  country_data = COUNTRY_DATA,nb_curves=NCURVES, 
                                  distance=Difference_function, 
                                  period4predicting=28, period4fitting =28,
                                  distance_type = "MSE")
      res_3 <- predict_prevalence(origin=as.Date(date), country=country,
                                  country_data = COUNTRY_DATA,nb_curves=NCURVES, 
                                  distance=Difference_function, 
                                  period4predicting=28, period4fitting =28,
                                  distance_type = "MSE", weights=(1:PERIOD_FOR_FITTING)/sum(1:PERIOD_FOR_FITTING))
      res_2 <- predict_prevalence(origin=as.Date(date), country=country,
                                  country_data = COUNTRY_DATA,nb_curves=NCURVES, 
                                  distance=correlation_function, 
                                  period4predicting=28, period4fitting =28,
                                  distance_type = "correlation")
      res_4 <- predict_prevalence(origin=as.Date(date), country=country,
                                  country_data = COUNTRY_DATA,nb_curves=NCURVES, 
                                  distance=correlation_function, 
                                  period4predicting=28, period4fitting =28,
                                  distance_type = "correlation", weights=(1:PERIOD_FOR_FITTING)/sum(1:PERIOD_FOR_FITTING))
      if(is.null(res_1) == FALSE){
        res = rbind(res, data.frame(date= as.Date(date), time = as.Date(date) + 1:28, "pred" = res_1$res$prevalence ,"sd" =  res_1$res$sd_prevalence, "type"="Sum of Squares", "obs"= res_1$res$Observed))
        samples_knn = rbind( samples_knn, data.frame(date= as.Date(date), time_pred = as.Date(date) + 1:28, pivot_wider(res_1$output,names_from = variable,values_from = value),  type="Sum of Squares", "obs"= res_1$res$Observed))
      }
      if(is.null(res_2) == FALSE){
        res = rbind(res, data.frame(date= as.Date(date),
                                    time = as.Date(date) + 1:28, 
                                    "pred" = res_2$res$prevalence ,"sd" =  res_2$res$sd_prevalence,
                                    "type"="Correlations", "obs"=  res_2$res$Observed))
        samples_knn= rbind( samples_knn, data.frame(date= as.Date(date), time_pred = as.Date(date) + 1:28, 
                                                    pivot_wider(res_2$output,names_from = variable,values_from = value),  type="Correlation", 
                                                    "obs"= res_2$res$Observed))
        
      }
      if(is.null(res_3) == FALSE){
        res = rbind(res, data.frame(date= as.Date(date), time = as.Date(date) + 1:28, "pred" = res_3$res$prevalence ,
                                    "sd" =  res_3$res$sd_prevalence, "type"="Weighted Sum of Squares", "obs"= res_3$res$Observed))
        samples_knn = rbind( samples_knn, data.frame(date= as.Date(date), time_pred = as.Date(date) + 1:28, 
                                                     pivot_wider(res_3$output,names_from = variable,values_from = value),
                                                     type="Weighted Sum of Squares", "obs"= res_3$res$Observed))
      }
      if(is.null(res_4) == FALSE){
        res = rbind(res, data.frame(date= as.Date(date), time = as.Date(date) + 1:28,
                                    "pred" = res_4$res$prevalence ,"sd" =  res_4$res$sd_prevalence, "type"="Weighted Correlations", 
                                    "obs"= res_4$res$Observed))
        samples_knn = rbind( samples_knn, data.frame(date= as.Date(date), time_pred = as.Date(date) + 1:28,
                                                     pivot_wider(res_r$output,names_from = variable,values_from = value), 
                                                     type="Weighted Correlations", "obs"= res_4$res$Observed))
      }
    }, error=function(e){cat("ERROR ",date, "\n")})
    
    
  }
  
  
}




u = melt(samples_knn %>% 
           filter(type=="Sum of Squares", date %in% as.Date(c("2020-07-06","2020-12-21", "2021-02-08"))), id.vars =c("type", "obs","date", "time_pred", "time"))%>%
  group_by(time, date) %>%
  summarize(obs = mean(obs), mean = mean(value), q97 =quantile( value,0.975), q25 = quantile(value,0.025))
ggplot(u %>% mutate(event = date +28))+
  theme_bw() + 
  geom_ribbon(data = u, aes(x=time + date, ymin=q25, ymax=q97, fill="Fitting"),  alpha=0.5) + 
  geom_ribbon(data = u %>% filter(time>0), aes(x=time + date , ymin=q25, ymax=q97, fill="Prediction"), alpha=0.5) + 
  scale_fill_manual(breaks = c("Fitting", "Prediction"), 
                    values=c("gray61", "gray32"))+
  geom_line(aes(x=time + date, y=mean, colour="Predicted (average)"),size=1) + 
  geom_line(aes(x=time + date, y=1e6 * obs, colour="Observed"), size=1) + 
  scale_colour_manual(breaks = c("Predicted (average)", "Observed"), 
                    values=c("black", "red"))+
  theme_bw() + xlab("Date") + ylab("Cases Per Million") +
  facet_wrap(.~ (date+28), scales = "free")
  


colnames(res)<-c("type", "date", "country", sapply(1:28, function(x){paste0("mean_predicted",x)}),
                 sapply(1:28, function(x){paste0("sd_predicted",x)}), "MSE", "correlation",
                 "coverage", "CI band width"
)

res$pred[which(res$type %in% c("EG", "AR"))] = 1e-6 * res$pred[which(res$type %in% c("EG", "AR"))]
res$obs[which(res$type %in% c("EG", "AR"))] = 1e-6 * res$obs[which(res$type %in% c("EG", "AR"))]
res$sd[which(res$type %in% c("EG", "AR"))] = 1e-6 * res$sd[which(res$type %in% c("EG", "AR"))]
res = res %>% mutate(MSE = (obs- pred)^2)

#### cap
res$MSE[which(res$MSE>1e6)] =NA
ggplot(res %>% filter(type !="EG") %>% group_by(type, date) %>% summarise(MSE = 1e6* sqrt(mean(MSE))),
       aes(x=type, y=MSE, fill=type))+
  geom_boxplot() + theme_bw() +
  theme(legend.text=element_text(size=16))+labs( fill="Method", size=14) + 
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1, size=0))+
  scale_fill_brewer(palette="Dark2") + xlab("Method") + ylab("MSE")
  



res = res %>% mutate(coverage = (obs <  pred + 2 * sd)  *(obs  > pred - 2 * sd)  )
ggplot(res %>% filter(type !="EG") %>% group_by(type, date) %>% summarise(coverage = mean(coverage)),
       aes(x=type, y=coverage, fill=type))+
  geom_boxplot() + theme_bw() +
  theme(legend.text=element_text(size=16))+labs( fill="Method", size=14) + 
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1, size=0))+
  scale_fill_brewer(palette="Dark2") + xlab("Method") + ylab("Coverage")

ggplot(res %>% filter(type !="EG") %>% group_by(type, date) %>% summarise(coverage = mean(coverage)),
       aes(x=type, y=coverage, fill=type))+
  geom_boxplot() + theme_bw() +
  theme(legend.text=element_text(size=16))+labs( fill="Method", size=14) + 
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1, size=0))+
  scale_fill_brewer(palette="Dark2") + xlab("Method") + ylab("Coverage")


ggplot(res_tot %>% dplyr::filter(country %in% c("China", "United States", "United Kingdom" , "Russia",  "Italy",
                                         "Greece"  ,"France" ,  "European Union" , "Vatican","Brazil" ), 
                                 date%in% as.Date(c("2020-06-01", "2020-09-07", "2020-12-21", "2021-01-25"))))+
  geom_boxplot(aes(x=date, y=MSE, fill=country)) + 
  scale_y_log10() + 
  facet_wrap(.~ type) + 
  theme_bw()

res_tot$`CI band width` = as.numeric(res_tot$`CI band width`)
library(RColorBrewer)
mycolors <- colorRampPalette(brewer.pal(8, "Spectral"))(length(unique(res_tot$date)))
ggplot(data=res_tot %>% dplyr::filter(country %in% c("China", "United States", "United Kingdom" , "Russia",  "Italy",
                                                "Greece"  ,"France" ,  "European Union" , "Trinidad and Tobago",
                                                "Uganda" ,"Brazil", "Sri Lanka" )),
       aes(x=country, y=coverage))+
  geom_boxplot(fill="grey38", outlier.shape = NA) + 
  theme_bw() +geom_jitter(aes(colour=date), position=position_jitter(0.2), alpha=0.5) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_manual(values = mycolors) + 
  facet_wrap(.~type)



ggplot(res_tot %>% group_by(type, date) %>% summarise(w = mean(`CI band width`, na.rm=TRUE)))+
  geom_point(aes(x=date, y=w, colour=type)) + 
  ylim(0, 0.001) + theme_bw()+
  theme(axis.text.x = element_text(angle = 180, vjust = 0.5, hjust=1))

#### Let's see if we have good coverage across countries
res_tot$date = as.factor(res_tot$date)
ggplot(res_tot %>% filter(type%in% c('Sum of Squares', 'AR')))+
  geom_point(aes(x=country, y=`CI band width`, colour=date)) + ylim(0, 5000)+
  scale_y_log10()+theme_bw() + 
  facet_wrap(.~type)


ggplot(res_tot %>% filter(type%in% c('Sum of Squares', 'AR')))+
  +   geom_point(aes(x=country, y=`CI band width`, colour=date)) + ylim(0, 5000)+
  +   scale_y_log10()+theme_bw() + 
  +   facet_wrap(.~type) + theme(axis.text.x = element_text(angle = 90))



res_UK = melt(res_tot %>% 
                filter(country %in% c( "United Kingdom"  )) %>%
                select(c("type", colnames(res_tot)[1:58])), id.vars=c("date", "country", "type"))
library(tidyr)
extract_numeric(years)
res_UK$time = parse_number(as.character(res_UK$variable))
res_UK["pred"] = "mean"
res_UK$pred[grepl( "sd", as.character(res_UK$variable), fixed = TRUE)] ="sd"


data_UK = merge(data.frame(res_UK %>% filter(pred=="mean")%>% select(value, date, type, time)), data.frame(res_UK %>% filter(pred=="sd") %>% select(value, date, type, time)), by = c("date", "time", "type"))
ggplot(data_UK %>% filter( date == "2021-02-08") ,
       aes(x=time, y=1e6 * value.x)
       ) +
  geom_line(colour="red", size=1) +
  theme_bw() +
  geom_ribbon(aes(ymin=1e6 * value.x  - 2 * 1e6 * value.y, ymax=1e6 * value.x  + 2 * 1e6 * value.y),colour="black", fill="grey", alpha=0.5, size=0.5) +
  facet_wrap(.~type)  
  
  

  geom_line() + 
  geom_line(aes(y=ymin), color="red",size=2) +
  geom_line(aes(y=ymax), color="red",size=2) +
  geom_line(aes(y=ymax), color="blue",size=2) +
  theme_bw() + ylab("Cases Per Million")





