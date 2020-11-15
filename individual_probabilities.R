## Based on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7273257/#APP1

BETA_PROFESSION <- c(0,0,0)
BETA_CONTACT <- 0

compute_infectiousness_probability <- function(sensitivity, prevalence, profession=1,
                                               high_risk_contact = 1,
                                               mask_wearing = 0,
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
  prev = sum(prevalence * (1 - sensitivity))
  #### We model the log odds ratio
  l = log(prev/(1-prev) * (1-mask_wearing))
  l = l + BETA_PROFESSION[profession + 1] + BETA_CONTACT * high_risk_contact
  p = 1/(1 + exp(-l))
  return( p )
}


compute_hospitalization_probability <- function(age, Pregnant, Chronic_Renal_Insufficiency,
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