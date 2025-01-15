init.gen <- function() {
  # Nomi dei posti
  yini.names <- c("N", "biomass_e")
  
  # Valori iniziali
  y_ini <- c(1,1)
  
  # Assegna i nomi ai valori iniziali
  names(y_ini) <- yini.names
  y_ini = y_ini[yini.names]
  
  # Ritorna i dati iniziali
  return(y_ini)
}

