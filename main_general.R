#### Main file where we compute the ouputs of interest


#compute_quanta_emission_rate <- function(cv, ci, inhalation, activity){
  ### This fucnction computes the quanta emission rate (per hour) given room parameters
#  D <- c(0.8 * 10^(-6), 1.8 * 10^(-6), 3.5 * 10^(-6), 5.5 * 10^(-6))
#  Vd <- 4/3 * pi * (D / 2)^3 
  ## cv = 10^(3)~10^(11) copies/mL, average 10^(7) copies/mL 
  ## ci = 3.66*10^(-5)
  ### 'activity' is a vector of 4 droplets' distribution of 4 different size during a certain activity
#  return(cv * ci * inhalation * sum(activity * Vd))
#}
## A more direct result from paper is that for resting, standing, light exercise, moderate exercise, and heavy exercise, the average 
## inhalation rate is 0.49, 0.54, 1.38, 2.35, and 3.30 m3 h−1

compute_quanta_emission_rate <- function(inhalation, Vocal) {
  ## alternative, Vd is the droplet volume concentration(breathing, speaking, loud speaking)
  ## cv average 10^(7), sd 10^(7*0.71)
  # Use normal distribution for Vd, cv, crna, cpfu
  ci <- 1 / (crna * cpfu)
  N <- 500
  cv <- rnorm(N, 10^ (7), 10^ (7 * 0.71))
  crna <- rnorm(N, 1.3 * 10^ (2), 1.3 * 10)
  cpfu <- rnorm(N, 2.1 * 10^ (2), 2.1 * 10)
  VD <- c(rnorm(N, 2 * 10^ (-3), 5 * 10^ (-4)), rnorm(N, 9 * 10^ (-3), 2 * 10^ (-3)), rnorm (N, 6 * 10^ (-2), 2 * 10^ (-2)))
  Vd <- numeric(length(N))
  V <- VD * Vocal
  for(i in 1:3) {
    Vd <- Vd + V[i]
  }
  q <- sum(cv * ci * Vd) / N
  return(q * inhalation)  ### Got a single number for average emission rate
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
  N <- 500
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
  P[,2] <- P[,2] / susceptible_people
  return(P)
}

compute_distribution_hospitalizations <- function(conditions, age, susceptible_people){
  #### This function coomputes the distribution of people that 
  #### will be hospitalized as a result of the event
  r <- runif(1000 * number_of_conditions)    # from csv data
  # need to import characteristic data for the probability of medical conditions
  N <- 500   #steps
  c <- numeric(length(N))
  for(i in 1:N) {
    for(j in length(susceptible_people)) {
      for(k in length(number_of_conditions)) {
        p <- sample(r, 1)
        if(p < P[k]) {
          c[i] <- c[i] + 1
          break
        }        ## P[k] is the imported csv data, the percentage of hospitalized with a certain medical condition and age group
      }
    }
  }
  H <- as.data.frame(table(c))
  H[,2] <- H[,2] / susceptible_people
}

compute_distribution_deaths<- function(){
  #### This function coomputes the distribution of people that 
  #### will die  as a result  of the event
  ### TO DO
}