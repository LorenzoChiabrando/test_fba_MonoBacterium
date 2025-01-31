
analyze_flux_sensitivity <- function(wd, 
                                     model_name, 
                                     mn_name,
                                     reactions_of_interest) {
  
  # Helper Functions
  extract_file_number <- function(filename, pattern) {
    as.numeric(gsub(pattern = pattern, replacement = "", x = basename(filename)))
  }
  
  validate_config <- function(n, config) {
    if(!exists("config") || is.null(config[[2]][[n]][[3]][1])) {
      stop(paste("Invalid config structure for n =", n))
    }
    return(config[[2]][[n]][[3]][1])
  }
  
  create_flux_plot <- function(data, reaction_name, parameter_name = "param_sen") {
    ggplot(data[data$Reaction == reaction_name,], 
           aes(x = as.numeric(Time), y = as.numeric(Flux), 
               group = config, color = param_sen)) +
      geom_line(alpha = 0.7, linewidth = 0.5) +
      scale_color_gradient(low = "#1d8", high = "#126", name = "Starv") +
      labs(title = reaction_name, x = "Time", y = "Flux") +
      theme_minimal() +
      theme(legend.position = "right",
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "gray90"),
            plot.title = element_text(hjust = 0.5))
    }
  
  process_flux_files <- function(wd, config) {
    flux_files <- sample(list.files(pattern = "\\.flux$", 
                                    path = paste0(wd, "/", model_name, "_analysis"), 
                                    full.names = TRUE))
    
    lapply(flux_files, function(file) {
      tryCatch({
        # Extract file number
        n <- extract_file_number(file, 
                                 pattern = paste("(", model_name, "-analysis-)|(-1-", 
                                                 mn_name, ".flux)", sep = ""))
        
        # Read and process flux data
        flux <- read.table(file, header = FALSE)
        reactions_names <- gsub("\\(|\\)", "", gsub("\\(", "_", flux[1, -c(1, 2)]))
        colnames(flux) <- flux[1, ]
        flux <- flux[-1, ]
        
        # Filter and transform flux data
        filtered_flux <- flux %>% 
          select(Time, any_of(reactions_of_interest))
        
        temp_flux <- filtered_flux %>%
          tidyr::gather(key = "Reaction", value = "Flux", -Time)
        
        # Add configuration information
        param_sen_value <- validate_config(n, config)
        
        cbind(temp_flux,
              config = rep(as.character(file), nrow(temp_flux)),
              param_sen = rep(param_sen_value, nrow(temp_flux)))
        
      }, error = function(e) {
        warning(paste("Error processing file", file, ":", e$message))
        return(NULL)
      })
    })
  }
  
  # Main Execution
  tryCatch({
    # Load configuration
    load(paste0(wd, "/", model_name, "_analysis/", model_name, "-analysis.RData"))
    
    # Process flux files
    sensitivity_flux <- process_flux_files(wd, config)
    sensitivity_flux <- Filter(Negate(is.null), sensitivity_flux)
    
    # Save results
    saveRDS(sensitivity_flux, 
            file = paste0(wd, "/", model_name, "_analysis/sensitivity_flux.rds"))
    
    # Combine results and create plots
    combined_flux <- do.call(rbind, sensitivity_flux)
    unique_react <- unique(combined_flux$Reaction)
    
    # Generate plots for each reaction
    plots <- lapply(unique_react, function(reaction) {
      create_flux_plot(combined_flux, reaction)
    })
    
    # Return results
    return(list(
      sensitivity_data = sensitivity_flux,
      combined_data = combined_flux,
      plots = plots
    ))
    
  }, error = function(e) {
    stop(paste("Error in flux sensitivity analysis:", e$message))
  })
}