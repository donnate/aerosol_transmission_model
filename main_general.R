#### Main file where we compute the ouputs of interest


compute_quanta_emission_rate <- function(cv, ci, inhalation, activity){
  ### This fucnction computes the quanta emission rate (per hour) given room parameters
  D <- c(0.8 * 10^(-6), 1.8 * 10^(-6), 3.5 * 10^(-6), 5.5 * 10^(-6))
  Vd <- 4/3 * pi * (D / 2)^3 
  ### 'activity' is a vector of 4 droplets' distribution of 4 different size during a certain activity
  return(cv * ci * inhalation * sum(activity * Vd))
}

compute_quanta_inhaled_per_person <- function(emission, removal, infectious_people, inhalation, volume, duration){
  ### This fucnction computes the number of quanta inhaled by each participant 
  ### over the course of the event
  n <- function(t) {emission * infectious_people / (volume * removal) * (1 - exp(- removal * t))}   
  ### 'inhalation' is a vector of all susceptible people's inhalation rate
  return(integrate(n, 0, duration) * inhalation)
}

compute_distribution_infections <- function(inhaled, susceptible_people){
  #### This function coomputes the distribution (not just the average) of people that 
  #### will be infected over the course of the event
  ### 'inhaled' is a vector of all susceptible people's inhaled quanta
  p <- 1 - exp(- inhaled)
  r <- runif(1000 * susceptible_people)           ### random number pool should be larger than number of people
  N = 500
  c <- numeric(length(N))
  for(i in 1:N){                                ## N random samples to show a probability distribution
    d <- sample(r, susceptible_people)
    for(j in length(susceptible_people)){
      if(d[j] <= p[j]){
        c[i] <- c[i] + 1
      }
    }
  }
  P <- as.data.frame(table(c))
  P[,2] <- P[,2] / N
  return(P)
}

compute_distribution_hospitalizations <- function(){
  #### This function coomputes the distribution of people that 
  #### will be hospitalized as a result of the event
  ### TO DO
}

compute_distribution_deaths<- function(){
  #### This function coomputes the distribution of people that 
  #### will die  as a result  of the event
  ### TO DO
}