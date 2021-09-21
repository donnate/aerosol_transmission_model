library(tidyverse)

data <- read.csv("~/Downloads/superspreader_events.csv")
data_all  = data %>% filter(is.na(Total.Pop.at.Event) == FALSE,
                is.na(Cases.at.and.or.Before.Event) == FALSE)

data_with_index  = data %>% filter(is.na(Total.Pop.at.Event) == FALSE,
                            is.na(Index.Cases) == FALSE,(Index.Cases) != "",(Index.Cases) != "?",
                            is.na(Secondary.Cases) == FALSE)
data_with_index$Total.Cases = as.numeric(data_with_index$Total.Cases)
data_with_index$Index.Cases = as.numeric(data_with_index$Index.Cases)
data_with_index$Cases.at.and.or.Before.Event = as.numeric(data_with_index$Cases.at.and.or.Before.Event)
data_with_index["Setting"] = c("Close", "Medium", "Close", "Medium", "Medium",
                               "Medium", "Close", "Close" ,"Close","Close",
                               "Close", "Medium", "Close", "Medium", "Close",
                               "Close", "Close", "Close", "Close", "Close",
                               "Medium", "Close", "Medium", "Medium", "Close", "Close")
ggplot(data_with_index,aes(x = Total.Pop.at.Event, 
                           y = Secondary.Cases /(Total.Pop.at.Event - Index.Cases))) + 
  geom_point() + 
  scale_x_log10() +theme_bw()+
  geom_smooth()+
  theme(legend.text=element_text(size=16)) + 
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1)) + 
  xlab("Event Size") + geom_text(aes(label=Setting1), hjust=-0.1,size=3,fontface = "bold",position=position_jitter(width=0.05,height=0.05))+
  ylab("Percentage of Secondary Infections") + facet_wrap(.~ Setting)

ggplot(data,aes(x = Total.Pop.at.Event, 
                y = 100*as.numeric(Cases.at.and.or.Before.Event)/(Total.Pop.at.Event))) + 
  geom_point() + 
  scale_x_log10() +theme_bw() + 
  geom_smooth()+
  theme(legend.text=element_text(size=16)) + 
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1)) + 
  xlab("Event Size") + 
  ylab("Percentage of Participants Infected")

data_agg = data_with_index
data_agg[is.na(data_agg$Index.Cases), "Index.Cases"] <- data_agg[is.na(data_agg$Index.Cases), "Cases.at.and.or.Before.Event"] - data_agg[is.na(data_agg$Index.Cases), "Secondary.Cases"] 
data_agg = data_with_index%>% 
  dplyr::mutate(prop_cases =  (Cases.at.and.or.Before.Event - Index.Cases) /(Total.Pop.at.Event - Index.Cases))
data_agg["F"] = sapply(1:data_agg$prop_cases, function(x){mean(data_agg$prop_cases <=x)})
ggplot(data_with_index) + 
  geom_density(aes(x= (Cases.at.and.or.Before.Event - Index.Cases) /(Total.Pop.at.Event - Index.Cases))) + 
  geom_point(aes(x= , )) + 
theme_bw()

