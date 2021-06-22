source("screening_efficiency.R")
source("relative_infectiousness.R")



test = list("Never lies"= compute_sensitivity(0.0, input = input),
            "Sometimes Lies"=compute_sensitivity(0.3, input = input),
            "50/50" =  compute_sensitivity(0.5, input = input),
            "Often Lies"=compute_sensitivity(0.8, input = input),
            "Always Lies" =  compute_sensitivity(1, input = input))
infectivity = cbind(read_csv("infectiousness_mean.txt",col_names = FALSE),
                    read_csv("infectiousness_low.txt",col_names = FALSE)[,2],
                    read_csv("infectiousness_up.txt",col_names = FALSE)[,2])
######## Then
colnames(infectivity) <- c("Days", "mean", "lower", "upper" )

p_incub = table(factor(sapply(1:5000, function(b){
  mu = rnorm(1,1.63, 0.12); 
  s = rnorm(1,0.5, 0.05); 
  ##### sample incubation rate
  incubation = ceiling(rlnorm(1,mu,s)) 
} ), levels=1:31))/5000

p_incub = p_incub[1:10]/sum(p_incub[1:10])


infectiousness = data.frame("Days"=0:30,
                            "mean"=c(0, apply((sapply(1:10, function(i){
                              as.numeric(p_incub[i]) * c(rep(0,max(i-4+1,0)), unlist(infectivity %>% filter(Days> -4) %>% select(mean), use.names = FALSE), rep(0,100))[1:30]
                            })),1,sum)),
                            "lower"=c(0, apply((sapply(1:10, function(i){
                              as.numeric(p_incub[i]) * c(rep(0,max(i-4+1,0)), unlist(infectivity %>% filter(Days> -4) %>% select(lower), use.names = FALSE), rep(0,100))[1:30]
                            })),1,sum)),
                            "upper"=c(0, apply((sapply(1:10, function(i){
                              as.numeric(p_incub[i]) * c(rep(0,max(i-4+1,0)),unlist(infectivity %>% filter(Days> -4) %>% select(upper), use.names = FALSE), rep(0,100))[1:30]
                            })),1,sum))
)
infectiousness = infectiousness %>% mutate(tentative = (upper + lower)/2)

ggplot(infectiousness) +
  geom_line(aes(x=Days, y=mean, colour="Infectiousness"), size=1) +
  geom_ribbon(aes(x=Days, ymin = lower, ymax=upper), alpha=0.5) + 
  theme_bw() + ylab("Probability (%)") + xlab("Days since infection") + 
  #geom_line(data=infectivity %>% filter(Days>-6), aes(x=Days + 5, y=mean, colour="Percentage Culture Positives"), size=1)+
  #geom_ribbon(data=infectivity %>% filter(Days>-6), aes(x=Days + 5, ymin = lower, ymax=upper), alpha=0.2) + 
  #geom_line(data=incubation, aes(x=Days, y=100*mean, colour="Density of \n Incubation Length"), size=1)+
#  geom_errorbar(data=incubation, aes(x=Days, ymin=100*lower, ymax=100*upper, colour="Density of \n Incubation Length")) +
  scale_colour_manual(labels  = c("Probability Density of \nIncubation Length","Infectiousness"), values = c("red","black"), drop=TRUE)+
  theme(legend.text=element_text(size=16))+labs(colour="Curve Type", size=14) + 
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1))


u = c()
proba = list("Never lies" =0,     "Sometimes Lies"=0.3, "50/50"  =0.5,
             "Often Lies" =0.8,    "Always Lies" =1.0 )
for (pp in names(test)){
  u = rbind(u,
            data.frame("Date of Infection"=-rev(0:(PERIOD_FOR_PREDICTING-1)), 
                       "p"=  test[[pp]]$sens, 
                       "infectiousness" = rev(as.vector(infectiousness$mean[1:PERIOD_FOR_PREDICTING])),
                       "infectiousness_q025" = rev(as.vector(infectiousness$lower[1:PERIOD_FOR_PREDICTING])),
                       "infectiousness_q975" = rev(as.vector(infectiousness$upper[1:PERIOD_FOR_PREDICTING])),
                       "p97" = test[[pp]]$q975,
                       "p2" = test[[pp]]$q025,
                       "p_sd" = test[[pp]]$sd,
                       "type" = proba[[pp]],
                       "proba_lie" = pp))
}



u = u %>% mutate(infectiousness_event = infectiousness * p )
xx = unlist(u %>%filter(proba_lie == "Never lies") %>% select(infectiousness_event)
            )
u ["relative_infectiousness_event"] = u["infectiousness_event"] /rep(xx,3)
u ["diff_infectiousness_event"] = u["infectiousness_event"] -rep(xx,3)
u = u %>% mutate(infectiousness_sd = 1/4 * (infectiousness_q975 - infectiousness_q025))
u = u %>% mutate(infectiousness_event_sd = infectiousness_sd  * p_sd + infectiousness_sd*p + p_sd * infectiousness)
ggplot(u, aes(x=`Date.of.Infection`)) +
  geom_line(aes(y=diff_infectiousness_event, colour=proba_lie), size=2)+
  #geom_line(aes(y=100-100*relative_infectiousness_event, colour=proba_lie), size=1)+
  theme_bw()+
  xlab("Date of Infection\n (with respect to Event Date)") + 
  theme(legend.text=element_text(size=16))+labs(colour="Probability of Lying", fill="95% CI Range", size=14) + 
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1)) + 
  ylab("Difference in Probability (%)") 


u$proba_lie = factor(u$proba_lie, levels=names(test))
ggplot(u, aes(x=`Date.of.Infection`)) +
  geom_ribbon(aes(x = `Date.of.Infection`, ymin=100 * `p2`, ymax=100 * p97, fill="Escapes Screening"),alpha=0.3)+
  geom_ribbon(aes(x = `Date.of.Infection`, ymin=`infectiousness_q025`, ymax=infectiousness_q975, fill="Infectious"),alpha=0.3)+
  geom_ribbon(aes(x = `Date.of.Infection`, ymin=infectiousness * p - 2 *infectiousness_event_sd , ymax=infectiousness * p + 2 *infectiousness_event_sd, fill="Infectious at the event"),alpha=0.3)+
  scale_fill_manual(breaks = c("Escapes Screening", "Infectious", "Infectious at the event"), values=c( "gray16", "gray61", "red")) +
  geom_line(aes(x = `Date.of.Infection`, y=100 * `p` , colour="Escapes Screening"), size=1.2,se = FALSE )+theme_bw() +
  geom_line(aes(x = `Date.of.Infection`, y=infectiousness, colour = "Infectious"), size=1.2)+
  geom_line(aes(x = `Date.of.Infection`, y=infectiousness * p, colour = "Infectious at the event"), size=1.2)+
  scale_colour_manual(breaks = c("Escapes Screening","Infectious", "Infectious at the event"),  values=c("black", "gray","red"))+
  theme(legend.text=element_text(size=16))+labs(colour="Average Probability", fill="95% CI Range", size=14) + 
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1)) + 
  xlab("Date of Infection\n (with respect to Event Date)") + 
  ylab("Probability (%)") + facet_wrap(proba_lie~., nrow=1)




##### Figure Appendix

comp = rbind(data.frame(test), data.frame(test00),data.frame(test05),data.frame(test10))
comp["type"] = c(rep(0.3,28), rep(0,28),rep(0.5,28),rep(1,28) ) 
comp["Probability Symptom Check Fail"] = c(rep("Sometimes Lies",28), rep("Never Lies",28),rep("50/50",28),rep("Always Lies",28) )
comp["Day"] = rep(1:28,4 ) 

ggplot(comp, aes(x=Day ))+ theme_bw() +
  geom_errorbar(aes( x = Day + 0.2*type,ymin=q025, ymax=q975, colour=`Probability Symptom Check Fail` ),size=1)+
  scale_colour_manual(breaks= c("Sometimes Lies", "Never Lies","50/50", "Always Lies"),  values=c("black", "gray", "purple", "red"))+
  theme(legend.text=element_text(size=16))+labs(colour="Average Probability", fill="95% CI Range", size=14) + 
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1)) + 
  xlab("Date of Infection\n (with respect to Event Date)") + 
  ylab("Probability") 
