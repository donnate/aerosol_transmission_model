## Based on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7273257/#APP1

BETA_PROFESSION <- c(0,0,0)
BETA_CONTACT <- 0
EFFICIENCY_MASK <- 0.7
compute_infectiousness_probability <- function(sensitivity, prevalence, profession=1,
                                               high_risk_contact = 1,
                                               nb_people_hh=1){
  ##### This function computes the probability for an individual to be infectious,
  ##### based on his/her personal information
  ##### --------------------------------------------------------
  ##### sensitivity       : array of size 15, that characterizes the sensitivity of 
  #####                     the antigen test as a function of time since infection
  ##### prevalence        : array of size 15 that characterizes the prevalence of the disease in
  #####                     the participant's region in the days before the test
  ##### profession        : category of the participant's occupational sector
  ##### high_risk_contact : has the participant had any high risk contact?
  ##### mask_wearing      : how likely has the participant been to wear a mask
  ##### nb_people_hh      : nb of people living in the same household as the participant 
  prev = prevalence 
  #### We model the log odds ratio
  l = log(prev/(1-prev))
  l = l + BETA_PROFESSION[profession + 1] + BETA_CONTACT * high_risk_contact
  p = 1/(1 + exp(-l))
  #print(p)
  return( p * (1 - sensitivity))
}


compute_hospitalization_probability <- function(age, race, ethnicity, smoking, 
                                                bmi,gender,
                                                diabetes){
  ### original 2.5 is assuming median income
  a =  2.5 + (race == "White") * 9.3 + (race == "Black") * 13.8 + (race == "Asian") * 26.6
  a = a + (ethnicity == "Hispanic") * 11.1 + (ethnicity == "Non Hispanic") * 34.3 
  a = a +  (smoking == "Current Smoker") * 5.6 + (smoking == "Former Smoker") * 11.2 + (smoking == "Non Smoker") * 13.3 
  a = a + ifelse(age>80, 80, age)/80 * 3.18 * 100/(7.43 -3)
  a = a + (bmi>=30) * (bmi - 30)/30  * 3.12 * 100/(7.43 -3) +  (bmi<30) * (30 - bmi)/30 * 2.94 * 100/(7.43 -3)
  a = a + (gender == "Male") * 7.7
  a = a + (diabetes == 1) * 15.8
  points = 250 + a
  return(points)
}

compute_death_probability <- function(age, Pregnant, Chronic_Renal_Insufficiency,
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
  a = a + Tobacco * (-0.055)
  a = a + Gender * (-0.116)
  return(1/(1+exp(-a)))
}




compute_hospitalization_probability_bis <- function(age, Pregnant, Chronic_Renal_Insufficiency,
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
  a = a + Tobacco * 0.007
  a = a + Cardiovascular_Disease *	(- 0.005)
  a = a + Asthma * 	(- 0.065)
  a = a + Gender * (- 0.121)
  return(1/(1+exp(-a)))
}

compute_death_probability_bis <- function(age, Pregnant, Chronic_Renal_Insufficiency,
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
  a = a + Tobacco * (-0.055)
  a = a + Gender * (-0.116)
  return(1/(1+exp(-a)))
}


