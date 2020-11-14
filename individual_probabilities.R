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
  l = log(prev/(1-prev) * (1-mask_wearing))  + BETA_PROFESSION[profession] + BETA_CONTACT * high_risk_contact
  p = 1/(1 + exp(-l))
  return( p )
}


compute_hospitalization_probability <- function(age, sex, race, diabetes){
  ##### This function computes the probability for an individual to have an adverse covid case
  ##### based on his/her personal information and comorbities
  ##### --------------------------------------------------------
  ##### sensitivity       : array of size 15, that characterizes the sensitivity of 
  #####                     the antigen test as a function of time since infection
  ##### prevalence        : array of size 15 that characterizes the prevalence of the disease in
  #####                     the participant's region in the days before the test
  ##### profession        : category of the participant's occupational sector
  ##### high_risk_contact : has the participant had any high risk contact?
  ##### mask_wearing      : how likely has the participant been to wear a mask
  return( 0.2 )
}


compute_death_probability <- function(age, sex, race, diabetes){
  ##### This function computes the probability for an individual to have an adverse covid case
  ##### leading to death based on his/her personal information and comorbities
  ##### --------------------------------------------------------
  ##### sensitivity       : array of size 15, that characterizes the sensitivity of 
  #####                     the antigen test as a function of time since infection
  ##### prevalence        : array of size 15 that characterizes the prevalence of the disease in
  #####                     the participant's region in the days before the test
  ##### profession        : category of the participant's occupational sector
  ##### high_risk_contact : has the participant had any high risk contact?
  ##### mask_wearing      : how likely has the participant been to wear a mask
  return( 0.2 )
}
