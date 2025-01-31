
library(devtools)
library(dplyr)
library(R.matlab)
library(ggplot2)
library(dplyr)
library(tidyr)
library(epimod)
library(patchwork)
library(ggplot2)
library(scales)
# remove.packages("epimod")
# install_github("https://github.com/qBioTurin/epimod", ref="epimod_pFBA")
# downloadContainers()

wd = getwd()

# Load custom R scripts for Flux Balance Analysis (FBA) functions and utilities
source(paste0(wd, "/epimod_FBAfunctions/R/FBAgreatmodeClass.R"))    # Core FBA class functions
source(paste0(wd, "/epimod_FBAfunctions/R/class_generation.R"))     # Functions for generating FBA models
source(paste0(wd, "/epimod_FBAfunctions/R/readMat.R"))              # MATLAB file reading utilities
source(paste0(wd, "/epimod_FBAfunctions/R/ex_bounds_module.R"))
source(paste0(wd, "/functions/compiling_mn.R"))

model_name = "BacPN"
mn_name = "Escherichia_coli_str_K_12_substr_MG1655"
org_name = "Escherichia_coli"
abbr = "E_coli"

source(paste0(wd, "/functions/", model_name, "_plot_analysis.R"))
source(paste0(wd, "/functions/", model_name, "_plot_sensitivity_trace.R"))
source(paste0(wd, "/functions/", model_name, "_plot_sensitivity_flux.R"))
source(paste0(wd, "/functions/save_sensitivity_results.R"))

# Define the models to be processed
models <- list(list(name = mn_name, type = org_name))

# Create a results directory if it does not already exist
dir.create(paste0(wd, "/results", showWarnings = FALSE))

# Define the output directory for results
dest_dir <- paste0(wd, "/results/")

# Loop through each model for FBA processing
for (model_info in models) {
  fba_name <- model_info$name    # Extract the model name
  model_type <- model_info$type    # Extract the model type
  
  # Construct file paths for FBA input
  fba_fname <- paste0(fba_name, ".txt")
  matfile <- paste0(wd, '/input/', model_type, "/", fba_name, ".mat")
  cat("METABOLIC NETWORK: ", fba_name, " ( ", model_type, " )\n")
  
  # Check if the .mat file exists
  if (!file.exists(matfile)) {
    stop(paste("File not found:", matfile))  # Stop execution if file is missing
  }
  
  # Generate the model for the specified .mat file
  cat("Generating model for:", fba_name, "\n")
  model <- FBA4Greatmod.generation(fba_mat = matfile)
  
  # Set biomass parameters for the FBA model
  model <- setBiomassParameters(model, bioMax = 1.172, bioMean = 0.489, bioMin = 0.083)
  
  # Prepare the output path for .rds files
  output_rds_path <- paste0(wd, "/input/", model_type)
  files_to_move <- list.files(pattern = "\\.rds$", full.names = TRUE)
  
  # Move .rds files to the appropriate directory
  if (length(files_to_move) > 0) {
    sapply(files_to_move, function(f) {
      file.rename(from = f, to = file.path(output_rds_path, basename(f)))
    })
  } else {
    cat("No .rds files found for", fba_name, "\n")
  }

}

### Setting Upper Bounds for FBA reactions and Not FBA reactions ###

# To compile metabolic model
#
# FBA Model Processing Functions
# 
# This script contains functions for processing FBA  models.
# It includes functionality for model generation, parameter setting, and file handling.

FBA_PARAMETERS <- list(biomass = list(max = 1.172, mean = 0.489, min = 0.083),
                       model_tags = list(name = mn_name, type = org_name))

process_fba_models(working_dir = NULL, list(FBA_PARAMETERS$model_tags))

# This function extracts all the specified EX_ reactions from the provided models (bacteria_files)
# and calculates their lower bound as (bacteria_counts*biomass)/MW.

# NotFBA exchange reactions LOWER BOUNDS are denoting metabolite concentration units. 
# Essential unit conversions include:
#
# - Biomass in grams of dry weight (gDW)
# - Environmental volume in milliliters (mL)
# - Metabolite concentration in millimolar (mM = mmol/L = 0.001 mmol/mL)
# 
# The environment is defined by volume (V) bacterial cell counts, and medium composition. 
# Conversion between molar concentrations and absolute quantities follows:

# Define the FBA models and associated counts
bacteria_models <- paste0(wd, "/results/", mn_name, ".txt")
bacteria_counts <- 1 # (cells)

molar = 0.001 # mmol/mL (1 mM)
V = 0.001 # mL (1 microL)
C = molar*V # mmol
biomass = 1e-12 # gWD

new_bounds = C

run_full_ex_bounds(
  extraction_output  = "extracted_ex_reactions.txt", # Output file for extracted reactions
  bacteria_files     = bacteria_models,              # List of FBA models to process
  output_dir         = "results_ex_reactions",       # Directory to store results
  bacteria_counts    = (bacteria_counts*biomass),
  non_fba_base_bound = C,
  reaction_version   = "r" # "r", "f" and "both" (default)
)

Bacteria_Parameters <- read.csv(paste0(wd, "/input/Bacteria_Parameters.csv"), head = F)
Bacteria_Parameters[1] <- 0.21 # starv
Bacteria_Parameters[2] <- 1 # dup
Bacteria_Parameters[3] <- 0.018 # death

write.table(Bacteria_Parameters,
            paste0(wd, "/input/Bacteria_Parameters.csv"),
            sep = ",", row.names = FALSE, col.names = FALSE, quote = F)

# Adjust General_functions.cpp

delta = 1e+10 # density

# Read the C++ file
cpp_file_path <- paste0(wd, "/input/General_functions.cpp")
cpp_content <- readLines(cpp_file_path)

# Define the new variables to insert
new_variables <- c(
  paste0("double V = ", V, "; // (mL) (1 microL)"),
  paste0("long long int delta = ", delta, "; // (cell/mL)")
)

# Find position after includes and namespace (looking for first empty line)
insert_position <- which(cpp_content == "")[1]

# Insert the new variables after the includes
cpp_content <- c(
  cpp_content[1:17],
  new_variables,
  cpp_content[(20):length(cpp_content)]
)

# Write the modified content back to the file
writeLines(cpp_content, cpp_file_path)

model.generation(net_fname = paste0(wd, "/net/", model_name, ".PNPRO"),
                 transitions_fname = paste0(wd, "/input/General_functions.cpp"),
                 fba_fname = paste0(wd, "/results/", mn_name, ".txt"))

system(paste0("mv ", 
              model_name, ".def ", 
              model_name, ".fbainfo ", 
              model_name, ".net ", 
              model_name, ".PlaceTransition ", 
              model_name, ".solver ", wd, "/net/"))

model.analysis(solver_fname = paste0(wd, "/net/", model_name, ".solver"),
               i_time = 0, f_time = 48, s_time = 0.5,
               parameters_fname = paste0(wd, "/input/initData.csv"),
               functions_fname = paste0(wd, "/functions/functions.R"),
               debug = T,
               fba_fname = paste0(wd, "/results/", mn_name, ".txt"),
               user_files = c(paste0(wd, "/input/Bacteria_Parameters.csv"),
                              paste0(wd, "/net/", model_name, ".fbainfo"),
                              paste0(wd, "/results_ex_reactions/EX_upper_bounds_FBA.csv"),
                              paste0(wd, "/results_ex_reactions/EX_upper_bounds_nonFBA.csv")))

p = plot_analysis()
# First ensure your combined plot is stored in a variable
p[[1]] + p[[2]] + (p[[3]] / p[[4]] / p[[5]])

Bacteria_Parameters <- read.csv(paste0(wd, "/input/Bacteria_Parameters.csv"), head = F)

starv_int = c(0.15, 0.5)
dup = Bacteria_Parameters[2]
death = Bacteria_Parameters[3]
  
init_data <- read.csv(paste0(wd, "/input/initData.csv"), head = F)
init_data[2, ] <- paste0("g; Bacteria_Parameters.csv; psensitivty; n = 1; min = ", 
                         starv_int[1], "; max = ", starv_int[2], 
                         "; dup = ", dup, "; death = ", death)
write.table(init_data,
            paste0(wd, "/input/initDataSensitivity.csv"),
            sep = ",", row.names = FALSE, col.names = FALSE, quote = F)

model.sensitivity(solver_fname = paste0(wd, "/net/", model_name, ".solver"),
                  fba_fname = paste0(wd, "/results/", mn_name, ".txt"),
                  i_time = 0, f_time = 72, s_time = 0.5, 
                  n_config = 250, debug = T,
                  target_value = c(abbr, paste0(abbr, "_biomass_e")),
                  parameters_fname = paste0(wd, "/input/initDataSensitivity.csv"),
                  parallel_processors = parallel::detectCores(),
                  functions_fname = paste0(wd, "/functions/functions.R"),
                  user_files = c(paste0(wd, "/net/", model_name, ".fbainfo"),
                                 paste0(wd, "/results_ex_reactions/EX_upper_bounds_FBA.csv"),
                                 paste0(wd, "/results_ex_reactions/EX_upper_bounds_nonFBA.csv")))

# results <- process_sensitivity_analysis(wd = wd, save_type = "both")
# Save only trace files
trace_results <- process_sensitivity_analysis(wd = wd, save_type = "trace")

# Get the plots
plots <- analyze_trace_sensitivity(wd = wd, model_name = model_name)

# Combine plots manually as needed
plots[[1]] + plots[[2]]

# Run the analysis
results <- analyze_flux_sensitivity(wd = wd,
                                    model_name = model_name,
                                    mn_name = mn_name,
                                    reactions_of_interest = c("EX_biomass_e_f", "EX_biomass_e_r",
                                                              "EX_glc_D_e_r", "EX_glc_D_e_f",
                                                              "EX_lcts_e_r", "EX_lcts_e_f"))

# Access results
# sensitivity_data <- results$sensitivity_data
# combined_data <- results$combined_data
plots <- results$plots

(plots[[1]] + plots[[2]]) / (plots[[3]] + plots[[4]]) / (plots[[5]] + plots[[6]]) 

# file.remove(list.files(path = wd, pattern = "\\.log$", full.names = TRUE))
# file.remove(list.files(path = wd, pattern = "\\ID$", full.names = TRUE))
# file.remove(list.files(path = wd, pattern = "\\StatusFile$", full.names = TRUE))
# gc()
