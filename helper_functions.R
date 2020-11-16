#### Set of Helper functions to compute all relevant parameters

extract_volume <- function(length, width, height, metric="m"){
  if (metric== "ft"){
    return(length * width * height * (0.3048)^3)
  }
  return(length * width * height)
}

extract_surface <- function(length, width, metric="m"){
  if (metric== "ft"){
    return(length * width  * (0.3048)^2)
  }
  return(length * width)
}

extract_first_order <- function(ventilation, control, decay=0.63, deposition=0.24){
  return( ventilation + decay + deposition + control)
}

extract_ventilation_rate_per_person <- function(volume, ventilation, control, n){
  # find inhalation
  return(volume * (ventilation + control) * 1000/36000/n)
}

lookup_ventilation <- function(setting, n_occupants, occupant_density, surface, volume){
  vent_rates = read_csv("ventilation_parameters.csv")
  print(c(setting))
  if (length(setting) == 1){
    t = vent_rates %>% dplyr::filter(category == setting)
    rp= t$rp
    ra = t$ra
  }else{
    t = vent_rates %>% dplyr::filter(category %in% setting)
    rp= max(t$rp)
    ra = max(t$ra)
  }
  

  vent_rate = n_occupants * rp + surface * ra
  return(vent_rate  * 3600 * 0.01 /volume)
}


