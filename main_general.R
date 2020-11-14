#### Main file where we compute the ouputs of interest


compute_quanta_emission_rate <- function(quanta_exhalation_rate, mask_efficiency , prop_mask,nb_infective_people){
  ### This function computes the quanta emission rate (per hour) given room parameters
  return(quanta_exhalation_rate * (1 - mask_efficiency * prop_mask) * nb_infective_people)
}

compute_quanta_concentation <- function(quanta_emission_rate, first_order_loss_rate, volume,   duration){
  ### This function computes the quanta emission rate (per hour) given room parameters
  return(quanta_emission_rate/first_order_loss_rate/volume * (1 - 1/first_order_loss_rate / duration)* (1 -exp(- first_order_loss * duration)) * nb_infective_people)
}

compute_quanta_inhaled_per_person <- function(quanta_concentration,  breathing_rate, duration, inhalation_mask_efficiency, prop_mask){
  ### This function computes the quanta emission rate (per hour) given room parameters
  return(quanta_concentration * breathing_rate * duration * (1 - inhalation_mask_efficiency * prop_mask) )
}
## A more direct result from paper is that for resting, standing, light exercise, moderate exercise, and heavy exercise, the average 
## inhalation rate is 0.49, 0.54, 1.38, 2.35, and 3.30 m3 h−1

# compute_quanta_emission_rate <- function(inhalation, Vocal) {
#   ## alternative, Vd is the droplet volume concentration(breathing, speaking, loud speaking)
#   ## cv average 10^(7), sd 10^(7*0.71)
#   # Use normal distribution for Vd, cv, crna, cpfu
#   ci <- 1 / (crna * cpfu)
#   N <- 500
#   cv <- rnorm(N, 10^ (7), 10^ (7 * 0.71))
#   crna <- rnorm(N, 1.3 * 10^ (2), 1.3 * 10)
#   cpfu <- rnorm(N, 2.1 * 10^ (2), 2.1 * 10)
#   VD <- c(rnorm(N, 2 * 10^ (-3), 5 * 10^ (-4)), rnorm(N, 9 * 10^ (-3), 2 * 10^ (-3)), rnorm (N, 6 * 10^ (-2), 2 * 10^ (-2)))
#   Vd <- numeric(length(N))
#   V <- VD * Vocal
#   for(i in 1:3) {
#     Vd <- Vd + V[i]
#   }
#   q <- sum(cv * ci * Vd) / N
#   return(q * inhalation)  ### Got a single number for average emission rate
# }



compute_quanta_inhaled_per_person <- function(emission, removal, infectious_people, inhalation, volume, duration){
  ### This fucnction computes the number of quanta inhaled by each participant 
  ### over the course of the event
  nb_quanta_inhaled <- function(t) {emission * infectious_people / (volume * removal) * (1 - exp(- removal * t))}   
  return(integrate(nb_quanta_inhaled, 0, duration) * inhalation)
}

compute_distribution_infections <- function(quanta_inhaled_per_person){
  #### This function coomputes the distribution (not just the average) of people that 
  #### will be infected over the course of the event

  p <- 1 - exp(-quanta_inhaled_per_person)        ### 'inhaled' is a vector of all susceptible people's inhaled quanta
  return(p)
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
  H[,2] <- H[,2] / N
}

compute_distribution_deaths<- function(){
  #### This function coomputes the distribution of people that 
  #### will die  as a result  of the event
  ### TO DO
}