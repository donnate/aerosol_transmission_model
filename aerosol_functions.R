#### Main file where we compute the ouputs of interest


compute_quanta_emission_rate <- function(activity, mask_efficiency ,
                                         prop_mask,nb_infective_people){
  ### This function computes the quanta emission rate (per hour) given room parameters
  ### Step 1. read in the file with the different categories
  quanta_rates = read_csv("quanta_emission_rates.csv")
  quanta_exhalation_rate = quanta_rates %>% dplyr::filter(Activity == activity) %>% dplyr::select(`Quanta/h`)
  return(quanta_exhalation_rate * (1 - mask_efficiency * prop_mask) * nb_infective_people)
}

compute_quanta_concentation <- function(quanta_emission_rate,
                                        first_order_loss_rate, volume,
                                        duration, nb_infective_people){
  ### This function computes the quanta emission rate (per hour) given room parameters
  return(quanta_emission_rate/first_order_loss_rate/volume * 
         (1 - 1/first_order_loss_rate / duration)* (1 -exp(- first_order_loss_rate * duration)) * 
           nb_infective_people)
}

compute_quanta_inhaled_per_person <- function(quanta_concentration,  breathing_rate,
                                              duration, inhalation_mask_efficiency, prop_mask){
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

