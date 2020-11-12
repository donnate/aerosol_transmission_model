#### Set of Helper functions to compute all relevant parameters

extract_volume <- function(length, width, height){
  return(length * width * height)
}

extract_surface <- function(length, width){
  return(length * width)
}

extract_first_order <- function(Ventilation, control){
  ventilation <- c(0.5, 3, 10)          # (natural, mechanical, mechanical)
  decay <- 0.63
  deposition <- 0.24
  return((ventilation[Ventilation] + decay + deposition) * control)
}

extract_ventilation_rate_per_person <- function(activity, age){
  # find inhalation rate using 'inhalation rate.csv'
  

  return()
}
