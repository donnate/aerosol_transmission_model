#### PARAMETER CHECK

decay = list()
for (UV in 0:10){
  for (rh in 20:70){
    for (temp in 10:30){
      decay = rbind(decay, data.frame("UV" = UV, "Relative Humidity" = rh, "Temperature"= temp, "Decay"=  (7.57 + 
                        1.41 * (temp-20.54)/10.66 +
                        0.0218 *(rh-45.24)/28.67 + 
                        7.55 *(UV *0.185-50) / 50 +
                        (temp-20.54)/10.66*(UV*0.185-50)/50*1.39 ) *60))  #https://www.dhs.gov/science-and-technology/sars-airborne-calculator
      
    }
  }
}


ggplot(decay %>% filter(Temperature ==20 )) +
  geom_tile(aes(x=`Relative.Humidity`, y=UV, fill=Decay)) +
  theme_classic() +
  scale_fill_gradientn(colours=hcl.colors(21, palette = "RdYlBu", rev=TRUE),na.value = "transparent",
                       breaks=0:20,
                       limits=c(0,20))+
  ylab("UV index") + xlab("Relative Humidity")

ggplot(decay %>% filter(Relative.Humidity ==40 )) +
  geom_tile(aes(x=Temperature, y=UV, fill=Decay)) +
  theme_classic() +
  scale_fill_gradientn(colours=hcl.colors(21, palette = "RdYlBu", rev=TRUE),na.value = "transparent",
                                        breaks=0:20,
                                        limits=c(0,20))+
  ylab("UV index") + xlab("Temperature (Celsius)")
  
