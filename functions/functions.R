
init.gen <- function() {
  # places' names
  yini.names <- c("E_coli", "E_coli_biomass_e")
  
  # initial marking
  y_ini <- c(1, 1) # places: "E_coli", "E_coli_biomass_e"
  
  names(y_ini) <- yini.names
  y_ini = y_ini[yini.names]
  
  return(y_ini)
}

psensitivty = function(n, min, max, dup, death) {
  starv = runif(n, min, max)
  return(matrix(c(starv, dup, death), ncol = 3))
}
