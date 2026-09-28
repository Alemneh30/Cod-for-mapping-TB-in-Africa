# =============================================================================
# MASTER SCRIPT: Mapping TB in Africa
# =============================================================================
# This script defines the main analysis parameters and runs all components
# of the replication workflow.
#
# For consistent results across analyses, user-defined parameters should be 
# changed here rather than separately in the individual scripts.
# =============================================================================


# 1. USER PARAMETERS ----------------------------------------------------------

# Number of posterior samples
nn <- 10000              # use 10000 for final results

# Spatial aggregation factor to generate grid-cell map outputs
mainaggfactor <- 2

# logistic regression model parameter
mu <- 0.025

# Prevalence scaling (prevalence per 1,000)
prevunit <- 1000

# Population-density threshold for pop. filtering
popt <- 5

# Run all sensitivity configurations?
allrun <- TRUE


# 2. MODEL CONFIGURATIONS -----------------------------------------------------

# popfilter = apply population-density filter
# mozout    = exclude Mozambique

if (allrun) {
  combinations <- expand.grid(popfilter = c(FALSE, TRUE),mozout = c(FALSE, TRUE))
} else {
  # Main specification only
  combinations <- data.frame(popfilter = FALSE,mozout = FALSE)
}


# 3. MAIN ANALYSIS ------------------------------------------------------------
start_time <- Sys.time()
message("Analysis started at: ", format(start_time, "%Y-%m-%d %H:%M:%S"))

message("1/5 Running main TB mapping analysis")

source("Mapping_TB_Africa_code.R")

# 4. ASSESSING THE EFFECTS OF EXCLUDING ETHIOPIAN DATA------------------------

message("2/5 Running prior robustness analysis")

source("Mapping_TB_Africa_robustness_Ethiopia.R")

# 5. PRIOR ROBUSTNESS ANALYSIS ------------------------------------------------

message("3/5 Running prior robustness analysis")

source("Mapping_TB_Africa_robustness.R")


# 6. REPORTED VERSUS RE-ANALYSIS -------------------------------------------------

message("4/5 Comparing reported vs updated values")

source("Analysis_reported_vs_updated.R")


# 7. FIGURES ------------------------------------------------------------------

message("5/5 Generating maps (main and SI)")

source("Plots_TB_Africa.R")


# 7. COMPLETE -----------------------------------------------------------------

end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

hours <- floor(elapsed / 3600)
minutes <- floor((elapsed %% 3600) / 60)
seconds <- round(elapsed %% 60)

message("All analyses completed successfully.")
message("Analysis finished at: ", format(end_time, "%Y-%m-%d %H:%M:%S"))
message(sprintf("Total running time: %02d h %02d min %02d sec",hours, minutes, seconds))

#End
