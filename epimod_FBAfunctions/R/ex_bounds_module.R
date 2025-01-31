# @Author: Mandringo02

# ============================================================
#  File: ex_bounds_module.R
#  Description: A module in R that extracts EX_ reactions from
#               text files and processes them to generate 
#               CSV files with dynamically computed bounds.
# ============================================================

# ------------------------------------------------------------------
# 1) EXTRACT EX_ REACTIONS
#    Reads the first line from each file in 'files' (a character 
#    vector of file paths), splits by space, filters reactions 
#    starting with "EX_", and writes them into 'output_file' 
#    with "_r" and "_f" appended.
# ------------------------------------------------------------------
extract_ex_reactions <- function(
  files,
  output_file = file.path(getwd(), "all_ex_reactions.txt")
) {
  # Initialize an empty vector to store all EX_ reactions
  all_ex_reactions <- character()
  
  # Iterate over each file in 'files'
  for (f in files) {
    if (!file.exists(f)) {
      warning(paste("File not found:", f, "- skipping."))
      next
    }
    # Read only the first line
    first_line <- readLines(f, n = 1, warn = FALSE)
    if (length(first_line) == 0) {
      # If the file is empty or has no first line, skip
      next
    }
    
    # Split the first line by spaces
    reactions <- unlist(strsplit(first_line, " "))
    
    # Filter only those that start with "EX_"
    ex_filtered <- reactions[grepl("^EX_", reactions)]
    
    # Accumulate them into the main vector
    all_ex_reactions <- c(all_ex_reactions, ex_filtered)
  }
  
  # Get unique entries and sort them
  all_ex_reactions <- unique(all_ex_reactions)
  all_ex_reactions <- sort(all_ex_reactions)
  
  # Create or overwrite the output file
  con <- file(output_file, open = "w")
  
  # For each EX_ reaction, write "_r" and "_f" variants
  for (rxn in all_ex_reactions) {
    writeLines(paste0(rxn, "_r"), con = con)
    writeLines(paste0(rxn, "_f"), con = con)
  }
  close(con)
  
  message("Extracted unique EX_ reactions written to: ", output_file)
}


# ------------------------------------------------------------------
# 2) PROCESS EX_ REACTIONS TO CREATE CSV FILES
#    This function reads a file ('input_file') containing lines 
#    like "EX_something_r" or "EX_something_f". 
#    Two output CSVs are generated:
#      - FBA CSV:   includes only reactions whose *base name* 
#                   is in 'fba_reactions'. We only keep the "_r" 
#                   version, using 'fba_upper_bound'.
#      - nonFBA CSV: includes other reactions (not in 'fba_reactions'),
#                    with bounds = 'not_shared_base_bound / population'.
#    Reactions containing "EX_biomass_e" are skipped entirely.
# ------------------------------------------------------------------
process_ex_reactions <- function(
  reaction_file,
  bacteria_files,
  output_dir           = getwd(),
  fba_reactions        = character(),
  bacteria_counts      = c(1),
  not_shared_base_bound   = 1000,
  fba_upper_bound      = 0.015,
  reaction_versions    = c("both", "r", "f") # 1) nuovo parametro
) {
  # 2) Forza uno dei 3 valori
  reaction_versions <- match.arg(reaction_versions)

  # (resto invariato)  
  n_bact <- length(bacteria_files)
  if (n_bact < 1) {
    stop("You must provide at least one bacterial model in 'bacteria_files'.")
  }
  if (!file.exists(reaction_file)) {
    stop("Reaction file does not exist: ", reaction_file)
  }
  if (length(bacteria_counts) != n_bact) {
    stop("Length of 'bacteria_counts' must match length of 'bacteria_files'.")
  }
  if (length(not_shared_base_bound) != 1 && length(not_shared_base_bound) != n_bact) {
    stop("'not_shared_base_bound' must be either a single numeric or a vector of length = number of bacteria.")
  }
  if (length(fba_upper_bound) != 1 && length(fba_upper_bound) != n_bact) {
    stop("'fba_upper_bound' must be either a single numeric or a vector of length = number of bacteria.")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  output_ub_fba_path    <- file.path(output_dir, "EX_upper_bounds_FBA.csv")
  output_ub_nonfba_path <- file.path(output_dir, "EX_upper_bounds_nonFBA.csv")
  
  f_fba    <- file(output_ub_fba_path, "w")
  f_nonfba <- file(output_ub_nonfba_path, "w")
  
  # Validazione e "replicazione" del base bound per le non-FBA
  if (length(not_shared_base_bound) == 1) {
    base_values <- rep(not_shared_base_bound, n_bact)
  } else {
    base_values <- not_shared_base_bound
  }

  first_line <- paste(c("base_upper_bounds", base_values), collapse = ",")
  writeLines(first_line, con = f_nonfba)
  
  get_base_name <- function(rxn) {
    sub("(_r|_f)$", "", rxn)
  }

  # Leggi tutte le reazioni
  reactions <- readLines(reaction_file, warn = FALSE)
  
  # Itera
  for (reaction in reactions) {
    # 1) Skip in base a reaction_versions
    if (reaction_versions == "r" && grepl("_f$", reaction)) next
    if (reaction_versions == "f" && grepl("_r$", reaction)) next

    # 2) Skip eventuale se è "EX_biomass_e"
    if (grepl("EX_biomass_e", reaction)) {
      next
    }
    
    # 3) Determina la base (senza _r/_f)
    base_rxn <- get_base_name(reaction)
    
    # 4) Se la base è in fba_reactions
    if (base_rxn %in% fba_reactions) {
      # ... e finisce in _r, allora scrivi su FBA
      if (grepl("_r$", reaction)) {
        ub_values <- numeric(n_bact)
        for (i in seq_len(n_bact)) {
          ub_values[i] <- if (length(fba_upper_bound) == 1) {
            fba_upper_bound
          } else {
            fba_upper_bound[i]
          }
        }
        line_to_write <- paste(c(reaction, ub_values), collapse = ",")
        writeLines(line_to_write, con = f_fba)
      }
    } else {
      # NON-FBA
      ub_values <- numeric(n_bact)
      for (i in seq_len(n_bact)) {
        ub_values[i] <- if (length(not_shared_base_bound) == 1) {
          not_shared_base_bound / bacteria_counts[i]
        } else {
          not_shared_base_bound[i] / bacteria_counts[i]
        }
      }
      line_to_write <- paste(c(reaction, ub_values), collapse = ",")
      writeLines(line_to_write, con = f_nonfba)
    }
  }
  
  # Chiudi i file
  close(f_fba)
  close(f_nonfba)
  
  message("Processing completed.\n  FBA file:    ", output_ub_fba_path,
          "\n  nonFBA file: ", output_ub_nonfba_path)
}





# ------------------------------------------------------------------
# 3) COMBINED FUNCTION
#    A single call that executes both:
#      - extraction of EX_ reactions from given text files
#      - processing of those reactions to create CSV files
# ------------------------------------------------------------------
run_full_ex_bounds <- function(
  extraction_output  = "all_ex_reactions.txt",
  bacteria_files,
  output_dir         = getwd(),
  fba_reactions      = character(),
  bacteria_counts    = c(1),
  not_shared_base_bound = 1000,
  fba_upper_bound    = 0.015,
  reaction_versions  = c("both", "r", "f")  
) {
  extract_ex_reactions(
    files       = bacteria_files,
    output_file = file.path(output_dir, extraction_output)
  )
  
  process_ex_reactions(
    reaction_file      = file.path(output_dir, extraction_output),
    bacteria_files     = bacteria_files,
    output_dir         = output_dir,
    fba_reactions      = fba_reactions,
    bacteria_counts    = bacteria_counts,
    not_shared_base_bound = not_shared_base_bound,
    fba_upper_bound    = fba_upper_bound,
    reaction_versions  = reaction_versions 
  )
}


