##### Computes the relative infectiousness


compute_relative_infectiousness <- function(input,period4predicting, plot=FALSE){
  test=  compute_sensitivity(input$p_lie, input, period4predicting = period4predicting +1)
  test = test[(nrow(test) - period4predicting-1): nrow(test),]
  PERIOD_FOR_PREDICTING = period4predicting
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
  } ), levels=0:31))/5000
  
  p_incub = p_incub[1:15]/sum(p_incub[1:15])
  
  
  infectiousness = data.frame("Days"=0:period4predicting,
                              "mean"=c(0, apply((sapply(1:10, function(i){
                                as.numeric(p_incub[i]) * c(rep(0,max(i-4+1,0)), unlist(infectivity %>% filter(Days> -4) %>% select(mean), use.names = FALSE), rep(0,100))[1:period4predicting]
                              })),1,sum)),
                              "lower"=c(0, apply((sapply(1:10, function(i){
                                as.numeric(p_incub[i]) * c(rep(0,max(i-4+1,0)), unlist(infectivity %>% filter(Days> -4) %>% select(lower), use.names = FALSE), rep(0,100))[1:period4predicting]
                              })),1,sum)),
                              "upper"=c(0, apply((sapply(1:10, function(i){
                                as.numeric(p_incub[i]) * c(rep(0,max(i-4+1,0)),unlist(infectivity %>% filter(Days> -4) %>% select(upper), use.names = FALSE), rep(0,100))[1:period4predicting]
                              })),1,sum))
                              )

  
  if (plot){
    ggplot(infectiousness) +
      geom_line(aes(x=Days, y=mean, colour="Infectiousness"), size=1) +
      geom_ribbon(aes(x=Days, ymin = lower, ymax=upper), alpha=0.5) + 
      theme_bw() + ylab("Probability (%)") + xlab("Days since infection") + 
      #geom_line(data=infectivity %>% filter(Days>-6), aes(x=Days + 5, y=mean, colour="Percentage Culture Positives"), size=1)+
      #geom_ribbon(data=infectivity %>% filter(Days>-6), aes(x=Days + 5, ymin = lower, ymax=upper), alpha=0.2) + 
      #geom_line(data=incubation, aes(x=Days, y=100*mean, colour="Density of \n Incubation Length"), size=1)+
      #geom_errorbar(data=incubation, aes(x=Days, ymin=100*lower, ymax=100*upper, colour="Density of \n Incubation Length")) +
      scale_colour_manual(labels  = c("Probability Density of \nIncubation Length","Infectiousness"), values = c("red","black"), drop=TRUE)+
      theme(legend.text=element_text(size=16))+labs(colour="Curve Type", size=14) + 
      theme(text = element_text(size=20),
            axis.text.x = element_text(angle=90, hjust=1))
  }

  infectiousness_all = data.frame("Date of Infection"=-rev(0:(period4predicting)), 
                                  "p"=  test$sens, 
                                  "infectiousness" = rev(as.vector(infectiousness$mean[1:(period4predicting+1)])),
                                  "infectiousness_q025" = rev(as.vector(infectiousness$lower[1:(period4predicting+1)])),
                                  "infectiousness_q975" = rev(as.vector(infectiousness$upper[1:(period4predicting+1)])),
                                  "p97" = c(test$q975),
                                  "p2" = c(test$q025),
                                  "p_sd" = c(test$sd),
                                  "type" = input$p_lie,
                                  "proba_lie" = input$p_lie) 
  
  infectiousness_all  = infectiousness_all  %>% mutate(infectiousness_event = infectiousness * p )
  infectiousness_all = infectiousness_all %>% mutate(infectiousness_sd = 1/4 * (infectiousness_q975 - infectiousness_q025))
  infectiousness_all = infectiousness_all %>% mutate(infectiousness_event_sd = infectiousness_sd  * p_sd + infectiousness_sd*p + p_sd * infectiousness)
  infectiousness_all = infectiousness_all %>% mutate(infectiousness_event_low = infectiousness * p - 2 *infectiousness_event_sd,
                                                     infectiousness_event_high = infectiousness * p + 2 *infectiousness_event_sd)
  
  infectiousness_all$infectiousness_event_low = sapply(infectiousness_all$infectiousness_event_low, function(x){max(x,0)})
  infectiousness_all$infectiousness_event_high = sapply(infectiousness_all$infectiousness_event_high, function(x){min(x,100)})
  
  if (plot){
    ggplot(infectiousness_all, aes(x=`Date.of.Infection`)) +
      geom_ribbon(aes(x = `Date.of.Infection`, ymin=100 * `p2`, ymax=100 * p97, fill="Escapes Screening"),alpha=0.3)+
      geom_ribbon(aes(x = `Date.of.Infection`, ymin=`infectiousness_q025`, ymax=infectiousness_q975, fill="Infectious"),alpha=0.3)+
      geom_ribbon(aes(x = `Date.of.Infection`, ymin=infectiousness_event_low , ymax=infectiousness_event_high, fill="Infectious at the event"),alpha=0.3)+
      scale_fill_manual(breaks = c("Escapes Screening", "Infectious", "Infectious at the event"), values=c( "gray16", "gray61", "red")) +
      geom_line(aes(x = `Date.of.Infection`, y=100 * `p` , colour="Escapes Screening"), size=1.2)+theme_bw() +
      geom_line(aes(x = `Date.of.Infection`, y=infectiousness, colour = "Infectious"), size=1.2)+
      geom_line(aes(x = `Date.of.Infection`, y=infectiousness * p, colour = "Infectious at the event"), size=1.2)+
      scale_colour_manual(breaks = c("Escapes Screening","Infectious", "Infectious at the event"),  values=c("black", "gray60","red"))+
      theme(legend.text=element_text(size=16))+labs(colour="Average Probability", fill="95% CI Range", size=14) + 
      theme(text = element_text(size=20),
            axis.text.x = element_text(angle=90, hjust=1)) + 
      xlab("Date of Infection\n (with respect to Event Date)") + 
      ylab("Probability (%)") 
  }
  return(infectiousness_all)
}


# infectivity <- 100 - data.frame(rbind(
# c(97.7, 58.8, 99.9),
# c(71.0, 29.6, 94.1),
# c(38.7, 18.4, 64.6),
# c(24.8, 14.0, 39.9),
# c(20.1, 12.5, 31.0),
# c(19.1, 12.0, 29.1),
# c(20.0, 12.8, 30.2),
# c(22.1, 14.5, 32.7),
# c(25.0, 16.7, 36.0),
# c(28.6, 19.5, 40.1),
# c(32.5, 22.6, 44.4),
# c(36.8, 26.2, 49.0),
# c(41.2, 30.0, 53.6),
# c(45.5, 33.8, 58.0),
# c(49.6, 37.6, 62.1),
# c(53.5, 41.2, 65.7),
# c(57.0, 44.6, 68.9),
# c(60.2, 47.8, 71.6),
# c(63.0, 50.8, 74.2)))


# infectivity = cbind(read_csv("infectiousness_mean.txt",col_names = FALSE),
#                     read_csv("infectiousness_low.txt",col_names = FALSE)[,2],
#                     read_csv("infectiousness_up.txt",col_names = FALSE)[,2])
# ######## Then
# colnames(infectivity) <- c("Days", "mean", "lower", "upper" )
# 
# 
# 
# # k t  = 5
# # (k - 1) t  = 4  => t = 1, k=4
# 
# p_incub = table(factor(sapply(1:5000, function(b){
#    mu = rnorm(1,1.63, 0.12); 
#    s = rnorm(1,0.5, 0.05); 
#    ##### sample incubation rate
#    incubation = ceiling(rlnorm(1,mu,s)) 
#  } ), levels=1:31))/5000
#  
# p_incub = p_incub[1:10]/sum(p_incub[1:10])
# 
# 
# infectiousness = data.frame("Days"=1:30,
#                 "mean"=apply((sapply(1:10, function(i){
#   as.numeric(p_incub[i]) * c( unlist(infectivity %>% filter(Days> -i) %>% select(mean), use.names = FALSE), rep(0,100))[1:30]
# })),1,sum),
# "lower"=apply((sapply(1:10, function(i){
#   as.numeric(p_incub[i]) * c( unlist(infectivity %>% filter(Days> -i) %>% select(lower), use.names = FALSE), rep(0,100))[1:30]
# })),1,sum),
# "upper"=apply((sapply(1:10, function(i){
#   as.numeric(p_incub[i]) * c( unlist(infectivity %>% filter(Days> -i) %>% select(upper), use.names = FALSE), rep(0,100))[1:30]
# })),1,sum))
# 
# 
# uncertainty_incub = sapply(1:5000, function(b){
#   print(b)
#   table(factor(sapply(1:5000, function(b){
#     mu = rnorm(1,1.63, 0.12); 
#     s = rnorm(1,0.5, 0.05); 
#     ##### sample incubation rate
#     incubation = ceiling(rlnorm(1,mu,s)) 
#   } ), levels=1:31))/5000
# })
# #######  Compute uncertainty estimates from the simulations (per bin)
# uncertainty_incub = uncertainty_incub / apply(uncertainty_incub, 2, sum)
# upper = apply(uncertainty_incub, 1, quantile, 0.975)
# lower = apply(uncertainty_incub, 1, quantile, 0.025)
# mean  = apply(uncertainty_incub, 1, mean)
# 
# incubation <- data.frame("Days" = 1:31,
#                          "mean" =mean ,
#                          "lower"= lower,
#                          "upper"= upper)
# infectivity= rep(0,60)
# infectivity_q975= rep(0,60)
# infectivity_q025 = rep(0,60)
# infectivity_q50 = rep(0,60)
# for(i in 1:31){
#   print(i)
#  uu = table(factor(sapply(1:5000, function(b){ceiling(rgamma(1, i, scale= runif(1,0.5,1.5)))}), levels=1:60))/5000
#  temp <- apply(uu, 1, mean)
#  infectivity = infectivity + temp/max(temp) * mean[i]
#  infectivity_q975 = infectivity_q975 + temp/max(temp) * upper[i]
#  infectivity_q025 = infectivity_q025 + temp/max(temp) * lower[i]
# }
# 
# write_csv(infectiousness, file =  "infectiousness.csv")


# colnames(infectivity) <- c("average", "q025", "q975")
# infectivity <- infectivity %>% mutate(average = average/max(average),
#                                       q025 = q025/max(q025),
#                                       q975 = q975/max(q975))
# 
# infectivity_matrix <- matrix(c(infectivity$average, rep(0,30)), 1)
# infectivity_matrix_025 <- matrix(c(infectivity$q025, rep(0,30)), 1)
# infectivity_matrix_975 <- matrix(c(infectivity$q975, rep(0,30)), 1)
# for(i in 1:30){
#   infectivity_matrix<- rbind(infectivity_matrix,
#                              c(rep(0,i), infectivity$average, rep(0, 30-i) ))
#   infectivity_matrix_025<- rbind(infectivity_matrix_025,
#                              c(rep(0,i), infectivity$q025, rep(0, 30-i) ))
#   infectivity_matrix_975<- rbind(infectivity_matrix_975,
#                              c(rep(0,i), infectivity$q975, rep(0, 30-i) ))
# }
# 
# p_incub = table(factor(sapply(1:5000, function(b){
#   mu = rnorm(1,1.63, 0.12); 
#   s = rnorm(1,0.5, 0.05); 
#   ##### sample incubation rate
#   incubation = ceiling(rlnorm(1,mu,s)) 
# } ), levels=1:31))/5000
# 
# p_incub = p_incub/sum(p_incub)
# 
# p_infectiousness = t(as.matrix(p_incub, nrow=1)) %*% infectivity_matrix
# p_infectiousness = as.vector(p_infectiousness)
# 
# p_infectiousness_025 = t(as.matrix(p_incub, nrow=1)) %*% infectivity_matrix_025
# p_infectiousness_025 = as.vector(p_infectiousness_025)
# p_infectiousness_975 = t(as.matrix(p_incub, nrow=1)) %*% infectivity_matrix_975
# p_infectiousness_975 = as.vector(p_infectiousness_975)

