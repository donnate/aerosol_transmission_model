#### Main file where we compute the ouputs of interest


compute_quanta_emission_rate <- function(activity, mask_efficiency ,
                                         prop_mask,nb_infective_people){
  ### This function computes the quanta emission rate (per hour) given room parameters
  ### Step 1. read in the file with the different categories
  quanta_rates = read_csv("quanta_emission_rates.csv", col_types = cols())
  quanta_exhalation_rate = quanta_rates %>% dplyr::filter(Activity == activity) %>% dplyr::select(`Quanta/h`)
  return(quanta_exhalation_rate * (1 - mask_efficiency * prop_mask) * nb_infective_people)
}

compute_quanta_concentation <- function(quanta_emission_rate,
                                        first_order_loss_rate, volume,                            
                                        duration, nb_infective_people, distance, vr){                 #### Added (distance, vr) to the original function
  ### This function computes the quanta emission rate (per hour) given room parameters
  ## return(quanta_emission_rate/first_order_loss_rate/volume * 
  ##        (1 - 1/first_order_loss_rate / duration)* (1 -exp(- first_order_loss_rate * duration)) * 
  ##          nb_infective_people)
  if(distance<2){   #### 2 metres/6 feet is tipping point for close contact, within close range we calculate the escalated concentration
      return(quanta_emission_rate/first_order_loss_rate/(pi * distance^2 * 1.5) *               #### Volume=pi*distance^2*H, H is set to be 1.5 m.
        (1 - exp(- first_order_loss_rate * distance/vr))* (1 -exp(- first_order_loss_rate * duration)))
  }
  if(distance>=2){
      return(quanta_emission_rate/first_order_loss_rate/volume * 
          (1 - 1/first_order_loss_rate / duration)* (1 -exp(- first_order_loss_rate * duration)) * 
            nb_infective_people)
  }
}

compute_quanta_inhaled_per_person <- function(quanta_concentration,  breathing_rate,
                                              duration, inhalation_mask_efficiency, prop_mask){
  ### This function computes the quanta emission rate (per hour) given room parameters
  return(quanta_concentration * breathing_rate * duration * (1 - inhalation_mask_efficiency * prop_mask) )
}