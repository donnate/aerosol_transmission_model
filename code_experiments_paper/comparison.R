PERRIOD_FOR_PREDICTING=28
setwd("~/Dropbox/aerosol_transmission_model/")
MASK_EFFICIENCY = 0.5  ### 50% is the recommended value
##### Set parameters
input = list("country" =  "United Kingdom",
             "date_event" = as.Date("2021-03-15"),
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
B=100000
N_TOT = 5000
volume = 86650
surface = 20000

country = input$country
PERIOD_FOR_FITTING = 14 
NCURVES = 200
MAX_DATE = input$date_event - 28
COUNTRY_DATA <- read.csv(file="https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv", header=T, sep=",")
COUNTRY_DATA$date <- (as.Date(COUNTRY_DATA$date, "%Y-%m-%d"))
VACCINATIONS <- read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations.csv")
VACCINATIONS_US <- read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/us_state_vaccinations.csv")
intersection_col = intersect(colnames(VACCINATIONS), colnames(VACCINATIONS_US))
VACCINATIONS <- rbind(VACCINATIONS[,intersection_col], VACCINATIONS_US[,intersection_col] %>% filter(location != "United States"))



PERIOD_FOR_PREDICTING = as.numeric(as.Date(input$date_event)- as.Date(max(MAX_DATE, na.rm = TRUE), fmt="%Y-%m-%d"))
DECAY = max(0, (7.57+ 
                  1.41* (input$temperature-20.54)/10.66 +
                  0.0218 *(input$RH-45.24)/28.67 + 
                  7.55 *((input$UV*0.185)-50) / 50 +
                  (input$temperature-20.54)/10.66*(input$UV*0.185-50)/50*1.40) *60)  #https://www.dhs.gov/science-and-technology/sars-airborne-calculator

BREATHING_RATE = 0.012 * 60 
DEPOSITION = 0.24
MASK_INHALATION_EFFICIENCY = 0.5
PRESSURE = 0.95
library(tidyverse)
data0 <- read_csv("results_R0.csv") %>% filter(type !=0)

ind = data0 %>% group_by(country, period4fitting, date) %>% tally()
data <- read_csv("results_R0.csv") %>% filter(type !=0)
data <- data[,2:128]
data <- rbind(data,
              read_csv("resultscorrelation.csv")%>% filter(type !=0))
data <- rbind(data,
              read_csv("resultscorrelation_weighted.csv")%>% filter(type !=0))
data <- rbind(data,
              read_csv("resultsMSE.csv")%>% filter(type !=0))
data <- rbind(data,
              read_csv("resultsMSE_weighted.csv")%>% filter(type !=0))
data2 =  read_csv("results_medianMSE.csv")%>% filter(type !=0)
data2["type"] = "Median-MSE"
data <- rbind(data, data2)

#### Enrich dataset with actual observation
#data[sapply(1:PERIOD_FOR_PREDICTING, FUN=function(x){paste0("Observed_",x)})]  = 0
data = as.data.frame(data)
dates = unique(data$date)

ind 



resAR = data2 %>% filter(type=="AR")
resMSE = (data2 %>% filter(type=="MSE"))[,1:43]

resAR2 = merge(resAR,resMSE, by=c("date", "country", "period4fitting") ) 

ggplot(resAR2, aes(x=as.factor(period4fitting), y= (1e6* mean_predicted1.y - Observed1)/Observed1)) + 
  geom_boxplot() + scale_y_log10() +
  theme_bw()+ xlab("") + ylab("Root Mean Square Error (RMSE)")+
  theme_bw() +  facet_wrap(.~ncurves.y, scales="free") +
  scale_color_manual(values = c( "red", "black")) + 
  theme(legend.text=element_text(size=16))+labs(colour="Data Type", fill="95% Prediction Range", size=14) + 
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1)) 


print(summary <- resAR2  %>% group_by(period4fitting, ncurves.y) %>%
        summarise(meanMSE=median(abs(1e6* mean_predicted1.y - 0.01 - Observed1)/(0.01+Observed1), na.rm=TRUE),
                  meanAR=median(abs(1e6* mean_predicted1.x - 0.01 - Observed1)/(0.01+Observed1), na.rm=TRUE)
        ))

print(summary <- resAR2  %>% group_by(period4fitting, ncurves.y) %>%
        summarise(meanMSE=median(abs(1e6* mean_predicted3.y - 0.01 - Observed3)/(0.01+Observed3), na.rm=TRUE),
                  meanAR=median(abs(1e6* mean_predicted3.x - 0.01 - Observed3)/(0.01+Observed3), na.rm=TRUE)
        ))

print(summary <- resAR2  %>% group_by(period4fitting, ncurves.y) %>%
        summarise(meanMSE=median(abs(1e6* mean_predicted14.y - 0.01 - Observed14)/(0.01+Observed14), na.rm=TRUE),
                  meanAR=median(abs(1e6* mean_predicted14.x - 0.01 - Observed14)/(0.01+Observed14), na.rm=TRUE)
        ))

print(summary <- resAR2  %>% group_by(period4fitting, ncurves.y) %>%
        summarise(meanMSE=median(abs(1e6* mean_predicted21.y - 0.01 - Observed21)/(0.01+Observed21), na.rm=TRUE),
                  meanAR=median(abs(1e6* mean_predicted21.x - 0.01 - Observed21)/(0.01+Observed21), na.rm=TRUE)
        ))

print(summary <- resAR2  %>% group_by(period4fitting, ncurves.y) %>%
        summarise(meanMSE=median(abs(1e6* mean_predicted7.y  - Observed7), na.rm=TRUE),
                  meanAR=median(abs(1e6* mean_predicted7.x - Observed7), na.rm=TRUE)
        ))


print(summary <- resAR2  %>% group_by(period4fitting, ncurves.y) %>%
        summarise(meanMSE=median(abs(1e6* mean_predicted21.y  - Observed21), na.rm=TRUE),
                  meanAR=median(abs(1e6* mean_predicted21.x - Observed21), na.rm=TRUE)
        ))


print(summary <- resAR2 %>% filter(country %in% c("United Kingdom", "France", "Italy", "United States", "Israel", 
                                                  "Spain")) %>% group_by(period4fitting, ncurves.y) %>%
  summarise(meanMSE=mean(abs(1e6* mean_predicted1.y - 0.01 - Observed1)/(0.01+Observed1), na.rm=TRUE),
meanAR=mean(abs(1e6* mean_predicted1.x - 0.01 - Observed1)/(0.01+Observed1), na.rm=TRUE)
))

print(summary <- data %>% group_by(period4fitting, ncurves.y) %>%
        summarise(meanMSE=median(abs(1e6* mean_predicted2.y - 0.01 - Observed2)/(0.01+Observed2), na.rm=TRUE),
                  meanAR=median(abs(1e6* mean_predicted1.x - 0.01 - Observed2)/(0.01+Observed2), na.rm=TRUE)
        ))


print(summary <- data %>% group_by(period4fitting, ncurves, type) %>%
        summarise_if(is.numeric,mean, na.rm=TRUE
        ))

data = merge(data, COUNTRY_DATA %>% dplyr::select(continent, location) %>% group_by(continent, location) %>% tally() , by.x="country",
             by.y="location", all.x=TRUE)

data2 = merge(data, ind, by=c("date", "country"))
print(summary2 <- data2 %>% group_by(period4fitting, ncurves, type) %>%
        summarise_if(is.numeric,median, na.rm=TRUE
        ))

print(summary <- data %>% group_by(period4fitting, ncurves, type) %>%
        summarise_if(is.numeric,quantile, 0.975, na.rm=TRUE
        ))
ggplot(summary %>% filter(type!="EG")) +
  geom_point(aes(x=as.factor(period4fitting), y=`RMSE`, 
                 colour=as.factor(ncurves)), size=2) + 
  facet_wrap(.~type) + theme_bw() +
  xlab("Number of training days") + labs(colour = "Number of k-NN\n curves")

summary %>% arrange(`W-RMSE` + `RMSE`)
hospital_names <- list(
  'AR'="AR",
  'correlation'="Correlation Distance  \n Mean-based prediction",
  'correlation_weighted'="Weighted Correlation Distance  \n Mean-based prediction",
  'MSE'="MSE Distance  \n Mean-based prediction",
  'MSE_weighted'="Weighted MSE Distance \n Mean-based prediction",
  'Median-MSE'="MSE Distance  \n Median-based prediction"
)
hospital_labeller <- function(variable,value){
  return(hospital_names[value])
}

ggplot(summary %>% filter(type!="EG")) +
  geom_point(aes(x=as.factor(period4fitting), y=`RMSE`, 
                 colour=as.factor(ncurves)), size=2) + 
  facet_wrap(.~type, labeller=hospital_labeller) + theme_bw() +
  xlab("Number of training days") +  theme(text = element_text(size=20),
                                           axis.text.x = element_text(angle=90, hjust=1))+
  labs(colour = "Number of\n k-NN curves") + theme(legend.text=element_text(size=16))


print(summary <- data %>% group_by(period4fitting, ncurves, type) %>%
        summarise_if(is.numeric,median, na.rm=TRUE
        ))

ggplot(summary %>% filter(type!="EG")) +
  geom_point(aes(x=as.factor(period4fitting),
                 y=`W-RMSE`,
                 colour=as.factor(ncurves)),
             size=3) + 
  facet_wrap(.~type) +  theme_bw() + 
  xlab("Number of training days") + labs(colour = "Number of k-NN\n curves")

ggplot(summary %>% filter(type!="EG")) +
  geom_point(aes(x=period4fitting, y=`RMSE`, colour=as.factor(ncurves))) + 
  facet_wrap(.~type) +  theme_bw()

ggplot(data %>% filter(ncurves %in% c(100, NA), period4fitting == 14), 
       aes(x=type, y=100 * `coverage`)) + 
  geom_boxplot(aes(fill = type), outlier.shape = NA) +  
  theme_bw()+ xlab("") + ylab("Coverage (%)")+
  theme_bw() +
  scale_color_manual(values = c( "red", "black")) + 
  theme(legend.text=element_text(size=16))+
  labs(fill="Data Type", size=14) + 
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=40, hjust=1)) 


ggplot(data2 %>% filter(ncurves %in% c(100, NA)), 
       aes(x=as.factor(period4fitting), y=abs(1e6*mean_predicted1 - Observed1)/Observed1))+ 
  geom_boxplot() + scale_y_log10() +
  theme_bw()+ xlab("") + ylab("Diff12")+
  theme_bw() +  facet_wrap(.~type, scales="free") +
  scale_color_manual(values = c( "red", "black")) + 
  theme(legend.text=element_text(size=16))+labs(colour="Data Type", fill="95% Prediction Range", size=14) + 
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1)) 

ggplot(summary, aes(x=as.factor(period4fitting), y=1e6 * sqrt(MSE))) + 
  geom_boxplot() + scale_y_log10() +
  theme_bw()+ xlab("") + ylab("Weighted RMSE)")+
  theme_bw() +  facet_wrap(.~type, scales="free") +
  scale_color_manual(values = c( "red", "black")) + 
  theme(legend.text=element_text(size=16))+labs(colour="Data Type", fill="95% Prediction Range", size=14) + 
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1))

ggplot(summary %>% filter(type!="EG"), aes(x=as.factor(period4fitting), y=1e6 * sqrt(MSE))) + 
  geom_violin() + 
  theme_bw()+ xlab("") + ylab(" RMSE)")+
  theme_bw() +  facet_wrap(.~type) +
  scale_color_manual(values = c( "red", "black")) + 
  theme(legend.text=element_text(size=16))+labs(colour="Data Type", fill="95% Prediction Range", size=14) + 
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1)) 

ggplot(data %>% filter(type!="EG"), aes(x=date, y=1e6 * sqrt(MSE), colour=country)) + 
  geom_point() + 
  theme_bw()+ xlab("") + ylab(" RMSE)")+
  theme_bw() +  facet_wrap(.~type) +
  theme(legend.text=element_text(size=16))+labs(colour="Data Type", fill="95% Prediction Range", size=14) + 
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1)) 


ggplot(summary, aes(x=as.factor(period4fitting), y=coverage)) + 
  geom_boxplot() + scale_y_log10() +
  theme_bw()+ xlab("") + ylab("Toot Mean Square Error (RMSE)")+
  theme_bw() +  facet_wrap(.~type, scales="free") +
  scale_color_manual(values = c( "red", "black")) + 
  theme(legend.text=element_text(size=16))+labs(colour="Data Type", fill="95% Prediction Range", size=14) + 
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1))


ggplot(summary, aes(x=type, y=1e6 * sqrt(MSE), colour=period4fitting)) + 
  geom_boxplot() + facet_wrap(.~period4fitting)+ scale_y_log10()


ggplot(summary, aes(x=type, y=correlation)) + 
  geom_boxplot() 

ggplot(summary, aes(x=type, y=1e6 * `max bandwidth`)) + 
  geom_boxplot() 

ggplot(summary, aes(x=1e6 * `mean prev`, y=1e6 * `mean bandwdith`, colour=type)) + 
  geom_point() 


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





