## test_pipeline.R
## This is a lightweight version of MLDmodel11full.R designed for CI/CD testing.
## It runs a minimal number of iterations to verify that the code executes without errors.

# --- 1. Configuration (Test Mode) ---
config <- list(
  model_number = "TEST_RUN",
  n_chains = 1,          # Single chain for speed
  n_iter = 50,           # Very few iterations
  burn_in = 0,           # No burn-in
  thining_rate = 1,
  nz = 20,               # Reduced augmented individuals

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
  max_nn = 10,           # Reduced neighbors for speed
  dumprate = 1000
)

# --- 2. Environment Setup ---

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
    "Error: The following required files are missing:\n",
    paste(missing_files, collapse = "\n"),
    "\n\nMake sure to run 'simulate_data.R' before running this test."
  ))
}

cat("Test Environment: All files found.\n")

# --- 3. Load Data ---

cat("Loading test data...\n")
statespace <- read.csv(config$file_statespace)
traps <- read.csv(config$file_traps)
captures <- read.csv(config$file_captures)
sex <- read.csv(config$file_sex)

Xsex <- sex[,2]
Xeffort <- NULL

## Get required SCR functions ###
source("e2dist.R")
source("SCRi.fn.par1-cheetah_sex.R")
source("scrDataWOeffort.R")

# --- 4. Prepare & Run Model (Wrapped in tryCatch) ---

cat("Running test pipeline...\n")

# Wrap entire execution in tryCatch to catch errors in data prep or model run
tryCatch({

  # 4a. Format Data
  cat("Formatting data...\n")
  scrMaraLionData <- scrData(
    traps = traps,
    captures = captures,
    statespace = statespace,
    Xsex = Xsex,
    Xeff = Xeffort
  )

  # 4b. Run Model
  cat("Running analysis...\n")
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

  cat("Test run completed successfully.\n")

}, error = function(e) {
  cat("Test run FAILED with error:\n")
  print(e)
  # Force non-zero exit code so CI fails
  quit(status = 1)
})
