plot_analysis <- function(wd = getwd(),
                          trace_file = NULL,
                          flux_file = NULL,
                          reactions_of_interest = c(
                            "EX_biomass_e_f", "EX_biomass_e_r", 
                            "EX_glc_D_e_r", "EX_glc_D_e_f",
                            "EX_lcts_e_r", "EX_lcts_e_f"),
                          custom_colors = c("black", "#098962", "#6e2076", "#583725", "#FF7700"),
                          ncol_places = 1,
                          ncol_reactions = 2,
                          base_size = 12,
                          base_family = "Arial") {
  
  # Function to list all .log files in the current directory
  get_log_files <- function() {
    log_files <- list.files(pattern = "\\.log$")
    return(log_files)
  }
  
  filename = get_log_files()
  
  process_log_file <- function(filename) {
    # Read the log file
    log_data <- readLines(filename)
    
    # Extract starvation rates and times
    starv_results <- sapply(log_data, function(line) {
      if (grepl("Starvation rate", line)) {
        # Extract starvation rate
        rate_match <- regexpr("---> *[0-9.]+(?=\\s*\\(pgDW/h\\))", line, perl=TRUE)
        # Extract time
        time_match <- regexpr("time ---> *[0-9.e+-]+", line, perl=TRUE)
        
        if (rate_match != -1 && time_match != -1) {
          rate_str <- regmatches(line, rate_match)
          time_str <- regmatches(line, time_match)
          
          # Clean up the extracted values
          rate <- as.numeric(gsub("---> *", "", rate_str))
          time <- as.numeric(gsub("time ---> *", "", time_str))
          
          return(c(rate=rate, time=time))
        }
      }
      return(c(rate=NA, time=NA))
    })
    
    # Extract duplication rates and times
    dup_results <- sapply(log_data, function(line) {
      if (grepl("Duplication rate", line)) {
        # Extract duplication rate
        rate_match <- regexpr("---> *[0-9.]+(?=\\s*\\(cell/h\\))", line, perl=TRUE)
        # Extract time
        time_match <- regexpr("time ---> *[0-9.e+-]+", line, perl=TRUE)
        
        if (rate_match != -1 && time_match != -1) {
          rate_str <- regmatches(line, rate_match)
          time_str <- regmatches(line, time_match)
          
          # Clean up the extracted values
          rate <- as.numeric(gsub("---> *", "", rate_str))
          time <- as.numeric(gsub("time ---> *", "", time_str))
          
          return(c(rate=rate, time=time))
        }
      }
      return(c(rate=NA, time=NA))
    })
    
    # Extract duplication rates and times
    death_results <- sapply(log_data, function(line) {
      if (grepl("Death rate", line)) {
        # Extract duplication rate
        rate_match <- regexpr("---> *[0-9.]+(?=\\s*\\(cell/h\\))", line, perl=TRUE)
        # Extract time
        time_match <- regexpr("time ---> *[0-9.e+-]+", line, perl=TRUE)
        
        if (rate_match != -1 && time_match != -1) {
          rate_str <- regmatches(line, rate_match)
          time_str <- regmatches(line, time_match)
          
          # Clean up the extracted values
          rate <- as.numeric(gsub("---> *", "", rate_str))
          time <- as.numeric(gsub("time ---> *", "", time_str))
          
          return(c(rate=rate, time=time))
        }
      }
      return(c(rate=NA, time=NA))
    })
    
    # Convert results to matrices for easier handling
    starv_matrix <- matrix(unlist(starv_results), ncol=2, byrow=TRUE)
    dup_matrix <- matrix(unlist(dup_results), ncol=2, byrow=TRUE)
    death_matrix <- matrix(unlist(death_results), ncol=2, byrow=TRUE)
    
    # Print diagnostic information
    cat("\nNumber of lines read:", length(log_data))
    cat("\nStarvation rates extracted:", sum(!is.na(starv_matrix[,1])))
    cat("\nDuplication rates extracted:", sum(!is.na(dup_matrix[,1])))
    cat("\nDeath rates extracted:", sum(!is.na(death_matrix[,1])))
    
    # Create data frames for both measurements
    # Starvation rate dataframe
    valid_starv <- which(!is.na(starv_matrix[,1]))
    time_points_starv <- seq_along(valid_starv)
    
    df_starv <- data.frame(
      Time = starv_matrix[valid_starv,2],
      TimePoint = time_points_starv,
      StarvationRate = starv_matrix[valid_starv,1],
      File = rep(filename, length(time_points_starv))
    )
    
    # Duplication rate dataframe
    valid_dup <- which(!is.na(dup_matrix[,1]))
    time_points_dup <- seq_along(valid_dup)
    
    df_dup <- data.frame(
      Time = dup_matrix[valid_dup,2],
      TimePoint = time_points_dup,
      DuplicationRate = dup_matrix[valid_dup,1],
      File = rep(filename, length(time_points_dup))
    )
    
    # Duplication rate dataframe
    valid_death <- which(!is.na(death_matrix[,1]))
    time_points_death <- seq_along(valid_death)
    
    df_death <- data.frame(
      Time = death_matrix[valid_death,2],
      TimePoint = time_points_death,
      DeathRate = death_matrix[valid_death,1],
      File = rep(filename, length(time_points_death))
    )
    
    # Return both dataframes in a list
    return(list(starvation = df_starv, duplication = df_dup, death = df_death))
  }
  
  # FILE = basename(flux_file)
  extract_organism_name <- function(FILE) {
    # Pattern to match Genus_species format
    # Looks for word_word within the filename
    pattern <- "([A-Z][a-z]+_[a-z]+)"
    matches <- regexpr(pattern, FILE, perl = TRUE)
    if (matches > 0) {
      # Extract the matched portion
      organism <- regmatches(FILE, matches)[[1]]
      return(organism)
    }
    return(NULL)
  }
  
  # Create a custom publication theme
  custom_theme <- function(base_size = 12, base_family = "Arial") {
    theme_minimal(base_size = base_size, base_family = base_family) %+replace%
      theme(
        # Text elements
        text = element_text(family = base_family, color = "black"),
        plot.title = element_text(
          size = rel(1.3),
          face = "bold",
          margin = margin(b = 15),
          hjust = 0
        ),
        plot.subtitle = element_text(
          size = rel(1.1),
          margin = margin(b = 10),
          hjust = 0
        ),
        axis.title = element_text(size = rel(1.1)),
        axis.text = element_text(size = rel(0.9)),
        
        # Legend formatting
        legend.position = "bottom",
        legend.title = element_text(size = rel(0.9)),
        legend.text = element_text(size = rel(0.8)),
        legend.box.spacing = unit(0.5, "lines"),
        
        # Grid lines and panel
        panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey90", fill = NA, linewidth = 0.5),
        
        # Facet formatting
        strip.text = element_text(
          size = rel(0.9),
          face = "bold",
          margin = margin(b = 5, t = 5)
        ),
        strip.background = element_rect(fill = "grey95", color = NA),
        
        # Plot margins
        plot.margin = margin(15, 15, 15, 15),
        
        complete = TRUE
      )
  }
  
  col_place = custom_colors[1]
  col_flux = custom_colors[2]
  col_starv = custom_colors[3]
  col_dup = custom_colors[4]
  col_death = custom_colors[5]
  
  # Input validation (same as before)
  if (is.null(trace_file) || is.null(flux_file)) {
    trace_files <- list.files(path = file.path(wd, "BacPN_analysis"), 
                              pattern = "\\.trace$", 
                              full.names = TRUE)
    flux_files <- list.files(path = file.path(wd, "BacPN_analysis"), 
                             pattern = "\\.flux$", 
                             full.names = TRUE)
    
    if (length(trace_files) == 0 || length(flux_files) == 0) {
      stop("No .trace or .flux files found and none provided")
    }
    
    trace_file <- trace_files[1]
    flux_file <- flux_files[1]
  }
  
  # Extract organism name
  organism_name <- extract_organism_name(basename(flux_file))
  
  if (is.null(organism_name)) {
    warning("Could not extract organism name from filename, using default")
    organism_name <- "Unknown_organism"
  }
  
  # Process starvation rate data
  all_data <- process_log_file(filename)
  
  # Enhanced starvation rate plot
  p_starv <- ggplot(all_data[["starvation"]], aes(x = Time, y = StarvationRate, color = File)) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1, alpha = 0.5) +
    labs(title = "Starvation Rate",
         x = "Time (h)",
         y = expression("Starvation Rate"~(pgDW/h))) +
    scale_colour_manual(values = col_starv) +
    scale_y_continuous() +
    scale_x_continuous() +
    custom_theme(base_size = base_size, base_family = base_family)
  
  p_dup <- ggplot(all_data[["duplication"]], aes(x = Time, y = DuplicationRate, color = File)) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1, alpha = 0.5) +
    labs(title = "Duplication Rate",
         x = "Time (h)",
         y = expression("Duplication Rate"~(cell/h))) +
    scale_colour_manual(values = col_dup) +
    scale_y_continuous() +
    scale_x_continuous() +
    custom_theme(base_size = base_size, base_family = base_family)
  
  p_death <- ggplot(all_data[["death"]], aes(x = Time, y = DeathRate, color = File)) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1, alpha = 0.5) +
    labs(title = "Death Rate",
         x = "Time (h)",
         y = expression("Death Rate"~(cell/h))) +
    scale_colour_manual(values = col_death) +
    scale_y_continuous() +
    scale_x_continuous() +
    custom_theme(base_size = base_size, base_family = base_family)
  
  # Process and enhance trace plot
  trace <- utils::read.table(trace_file, header = TRUE) %>%
    tidyr::gather(key = "Places", value = "Marking", -Time)
  
  p_pl <- ggplot(trace, aes(Time, Marking)) + 
    geom_line(linewidth = 0.75, color = col_place) +
    labs(title = "Place Markings",
         x = "Time (h)",
         y = "Marking") +
    facet_wrap(~Places, ncol = ncol_places, scales = "free_y") +
    custom_theme(base_size = base_size, base_family = base_family) +
    theme(strip.text = element_text(face = "bold"),
          panel.spacing = unit(1.5, "lines"))
  
  # Process and enhance flux plot
  flux <- utils::read.table(flux_file, header = TRUE) %>%
    tidyr::gather(key = "Reaction", value = "Flux", -Time) %>%
    mutate(Organism = organism_name)
  
  flux_filtered <- flux %>% filter(Reaction %in% reactions_of_interest)
  
  p_ex <- ggplot(flux_filtered, aes(x = Time, y = Flux, color = Organism)) +
    geom_line(linewidth = 1) +
    labs(title = "Metabolic Fluxes", x = "Time (h)", y = "Flux (mmol/gDW*h)") +
    scale_colour_manual(values = col_flux) +
    facet_wrap(~Reaction, ncol = ncol_reactions, scales = "free_y") +
    custom_theme(base_size = base_size, base_family = base_family) +
    theme(strip.text = element_text(face = "bold"),
          panel.spacing = unit(1.5, "lines"))
  
  # Return enhanced plots
  return(list(
    places = p_pl,
    exchange = p_ex,
    starvation = p_starv,
    duplication = p_dup,
    death = p_death
  ))
}
