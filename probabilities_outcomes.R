# Based on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7273257/#APP1

infectiousness_probability <- function(){
  
}


hospitalization_probability <- function(age, Pregnant, Chronic_Renal_Insufficiency,
                            Diabetes, Immunosuppression, COPD, Obesity,
                            Hypertension, Tobacco, Cardiovascular_Disease,
                            Asthma, Gender){
  a = age * 1.000 
  a = a + Pregnant * 0.172
  a = a +  Chronic_Renal_Insufficiency * 	0.167
  a = a + Diabetes * 0.165
  a = a +  Immunosuppression *	0.139
  a = a +  COPD	* 0.094
  a = a + Obesity * 0.083
  a = a +  Hypertension	 * 0.039
  a = a + Tobacco 	* 0.007
  a = a + Cardiovascular_Disease *	- 0.005
  a = a + Asthma * 	- 0.065
  a = a + Gender * -0.121
  return(1/(1+exp(-a)))
}

death_probability <- function(age, Pregnant, Chronic_Renal_Insufficiency,
                                        Diabetes, Immunosuppression, COPD, Obesity,
                                        Hypertension, Tobacco, Cardiovascular_Disease,
                                        Asthma, Gender){
  a = age * 1.000 
  a = a + Immunosuppression * 0.175
  a = a + Asthma * 0.114
  a = a + Chronic_Renal_Insufficiency	 * 0.110
  a = a + Obesity	 *  0.104
  a = a + Hypertension * 0.087
  a = a + Pregnant * 0.087
  a = a + Diabetes * 0.077
  a = a + COPD * 	0.042
  a = a + Cardiovascular_Disease * 	0.027
  a = a + Tobacco_Use * -0.055
  a = a + Gender * -0.116
  return(1/(1+exp(-a)))
}
