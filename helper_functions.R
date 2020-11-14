#### Set of Helper functions to compute all relevant parameters

extract_volume <- function(length, width, height){
  return(length * width * height)
}

extract_surface <- function(length, width){
  return(length * width)
}

extract_first_order <- function(ventilation, control, decay=0.63, deposition=0.24){
  return( ventilation + decay + deposition + control)
}

extract_ventilation_rate_per_person <- function(volume, ventilation, control, n){
  # find inhalation
  return(volume * (ventilation + control) * 1000/36000/n)
}
