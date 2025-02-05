plot_analysis <- function(wd = getwd(),
                          trace_file = NULL,
                          flux_file = NULL,
                          reactions_of_interest,
                          custom_colors = c("black", "#098962", "#6e2076", "#583725", "#FF7700"),
                          ncol_places,
                          ncol_reactions,
                          place2plot,
                          base_size = 12) {
  
  # Function to list all .log files in the current directory
  get_log_files <- function() {
    log_files <- list.files(pattern = "\\.log$")
    return(log_files)
  }
  
  filename = get_log_files()
  
  process_log_file <- function(filename) {
    
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
  custom_theme <- function(base_size = 12) {
    theme_minimal(base_size = base_size) %+replace%
      theme(
        # Text elements
        text = element_text(color = "black"),
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
  if ( is.null(trace_file) || is.null(flux_file) ) {
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
  p_starv <- ggplot(all_data[["starvation"]], 
                    aes(x = as.double(Time), y = StarvationRate, color = File)) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1, alpha = 0.5) +
    labs(title = "Starvation Rate",
         x = "Time (h)",
         y = expression("Starvation Rate"~(pgDW/h))) +
    scale_colour_manual(values = col_starv) +
    scale_y_continuous() +
    scale_x_continuous() +
    custom_theme(base_size = base_size)
  
  p_dup <- ggplot(all_data[["duplication"]], aes(x = Time, y = DuplicationRate, color = File)) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1, alpha = 0.5) +
    labs(title = "Duplication Rate",
         x = "Time (h)",
         y = expression("Duplication Rate"~(cell/h))) +
    scale_colour_manual(values = col_dup) +
    scale_y_continuous() +
    scale_x_continuous() +
    custom_theme(base_size = base_size)
  
  p_death <- ggplot(all_data[["death"]], aes(x = Time, y = DeathRate, color = File)) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1, alpha = 0.5) +
    labs(title = "Death Rate",
         x = "Time (h)",
         y = expression("Death Rate"~(cell/h))) +
    scale_colour_manual(values = col_death) +
    scale_y_continuous() +
    scale_x_continuous() +
    custom_theme(base_size = base_size)
  
  trace <- utils::read.table(trace_file, header = TRUE) %>%
    tidyr::gather(key = "Places", value = "Marking", -Time) %>%
    # Convert Places to factor to maintain order
    mutate(Places = factor(Places, levels = unique(Places)))
  
  # Create a professional theme
  professional_theme <- function(base_size = 12) {
    theme_minimal(base_size = base_size) +
      theme(
        # Text elements
        text = element_text(color = "black"),
        plot.title = element_text(size = rel(1.2), face = "bold", hjust = 0.5),
        axis.title = element_text(size = rel(1.1), face = "bold"),
        axis.text = element_text(size = rel(1)),
        
        # Panel elements
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        panel.grid.minor = element_line(color = "gray95", linewidth = 0.2),
        panel.border = element_rect(fill = NA, color = "gray80", linewidth = 0.5),
        
        # Facet elements
        strip.text = element_text(face = "bold", size = rel(1)),
        strip.background = element_rect(fill = "gray95", color = NA),
        panel.spacing = unit(2, "lines"),
        
        # Legend elements
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        
        # Margins
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
      )
  }
  
  # Assuming place2plot is your vector of place names to plot
  variable_names <- setNames(
    list(
      paste0(levels(trace$Places)[1], " (cell)"),
      paste0(levels(trace$Places)[2], " (pgDW)")
    ),
    place2plot[1:2]  # Using the names from place2plot vector
  )
  
  variable_labeller <- function(variable,value){
    return(variable_names[value])
  }
  
  p_pl <- ggplot(trace, aes(Time, Marking)) + 
    geom_line(linewidth = 0.75, color = "darkgrey") +
    geom_point(size = 1, color = "#908432", alpha = 0.25) +  # Added points for better readability
    labs(title = "Markings Over Time",
         x = "Time (h)",
         y = "Marking") +
    # Use facet_grid instead of facet_wrap
    facet_grid(Places ~ ., scales = "free_y", labeller = variable_labeller) +
    professional_theme(base_size = 12) +
    # Add subtle background shading for visual interest
    geom_ribbon(aes(ymin = min(Marking), ymax = Marking), 
                alpha = 0.1, fill = "#908432")
  
  # Process and enhance flux plot
  flux <- utils::read.table(flux_file, header = TRUE) %>%
    tidyr::gather(key = "Reaction", value = "Flux", -Time) %>%
    mutate(Organism = organism_name)
  
  # Process the data as before
  flux_filtered <- flux %>%
    filter(Reaction %in% reactions_of_interest) %>%
    # Split reactions into forward and reverse
    mutate(
      Type = ifelse(grepl("_f$", Reaction), "forward", "reverse"),
      Base_Reaction = sub("_[fr]$", "", Reaction)
    )
  
  # Function to create individual plot
  create_subplot <- function(data, reaction_base, y_max) {
    f_data <- data %>% filter(Base_Reaction == reaction_base)
    
    ggplot(f_data, aes(x = Time, y = Flux, color = Organism)) +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = col_flux) +
      coord_cartesian(ylim = c(0, y_max)) +
      labs(x = "Time (h)", y = "Flux (mmol/gDW*h)") +
      facet_wrap(~Reaction, ncol = 1) +
      custom_theme(base_size = base_size) +
      theme(
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(1.5, "lines"),
        legend.position = "bottom"
      )
  }
  
  # Get unique base reactions
  base_reactions <- unique(sub("_[fr]$", "", reactions_of_interest))
  
  # Create list to store plots
  plot_list <- list()
  
  # Create plots for each base reaction
  for(reaction in base_reactions) {
    # Get maximum y value for this reaction pair
    y_max <- flux_filtered %>%
      filter(Base_Reaction == reaction) %>%
      pull(Flux) %>%
      max(na.rm = TRUE)
    
    # Create subplot and add to list
    plot_list[[reaction]] <- create_subplot(flux_filtered, reaction, y_max)
  }
  
  # Combine plots using patchwork
  final_plot <- wrap_plots(plot_list, ncol = ncol_reactions) +
    plot_annotation(
      title = "Metabolic Fluxes",
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
      )
    )
  
  # Add shared legend at the bottom
  p_ex <- final_plot & theme(legend.position = "bottom")
  
  # Return enhanced plots
  return(list(
    places = p_pl,
    exchange = p_ex,
    starvation = p_starv,
    duplication = p_dup,
    death = p_death
  ))
}
