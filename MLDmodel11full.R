## MLDmodel11full.R - Main Control Script
## Refactored for Robustness and Portability

# --- 1. Configuration Section ---
# Adjust these parameters as needed
config <- list(
  model_number = "BD2005sps",
  n_chains = 3,
  n_iter = 120000,
  burn_in = 20000,
  thining_rate = 1,
  nz = 100, # Number of zeros (augmented individuals)

  # Input files
  file_statespace = "new_mask_file.csv",
  file_traps = "Traps.csv",
  file_captures = "Capture.csv",
  file_sex = "Sex.csv",

  # Model settings
  theta = 1,
  Msigma = 1,
  Mb = 0,
  Msex = 1,
  Msexsigma = 1,

  # Spatial settings
  coord_scale = 1000,
  area_per_pixel = 0.336,
  thin_statespace = 1,
  max_nn = 40,
  dumprate = 1000
)

# --- 2. Environment Setup ---

# Remove hardcoded setwd()
# Instead, verify we are in a directory that contains the necessary files.
# If running interactively, the user should setwd() to the project root.

required_files <- c(
  config$file_statespace,
  config$file_traps,
  config$file_captures,
  config$file_sex,
  "e2dist.R",
  "SCRi.fn.par1-cheetah_sex.R",
  "scrDataWOeffort.R"
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(paste0(
    "Error: The following required files are missing in the current directory:\n",
    paste(missing_files, collapse = "\n"),
    "\n\nTip: If you need dummy data for testing, run 'simulate_data.R'."
  ))
}

cat("All required files found. Proceeding...\n")

# --- 3. Load Data ---

cat("Loading data...\n")
statespace <- read.csv(config$file_statespace)
traps <- read.csv(config$file_traps)
captures <- read.csv(config$file_captures)
sex <- read.csv(config$file_sex)

#effort <- read.csv("maralionEffort1.csv") # Commented out in original

## Extract specific information from the input files ##
Xsex <- sex[,2]
#Xeffort <- effort[,4:ncol(effort)]
Xeffort <- NULL # Explicitly set to NULL since it's commented out

## Get required SCR functions from the directory ###
source("e2dist.R")
source("SCRi.fn.par1-cheetah_sex.R")
source("scrDataWOeffort.R")

# --- 4. Prepare Data Object ---

cat("Formatting data...\n")
scrMaraLionData <- scrData(
  traps = traps,
  captures = captures,
  statespace = statespace,
  Xsex = Xsex,
  Xeff = Xeffort
)

# --- 5. Run Model ---

cat(paste0("Running model: ", config$model_number, " with ", config$n_chains, " chains...\n"))

scrMaraLionAnal <- SCRi.fn.par1(
  scrMaraLionData,
  modelno = config$model_number,
  nc = config$n_chains,
  ni = config$n_iter,
  burn = config$burn_in,
  skip = config$thining_rate,
  nz = config$nz,
  theta = config$theta,
  Msigma = config$Msigma,
  Mb = config$Mb,
  Msex = config$Msex,
  Msexsigma = config$Msexsigma,
  Xsex = Xsex,
  ss.prob = NULL,
  coord.scale = config$coord_scale,
  area.per.pixel = config$area_per_pixel,
  thinstatespace = config$thin_statespace,
  maxNN = config$max_nn,
  dumprate = config$dumprate
)

cat("Analysis complete.\n")
