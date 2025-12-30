# Debug script
cat("Checking dimensions...\n")

# Find a chain directory
chain_dirs <- grep("CH[0-9]+", list.files(), value = TRUE)
if(length(chain_dirs) > 0) {
  ss_file <- list.files(chain_dirs[1], pattern = "SSunscaled", full.names = TRUE)[1]
  ss <- read.csv(ss_file)
  cat("SSunscaled rows (nG):", nrow(ss), "\n")
} else {
  cat("No chain dir found.\n")
}

if(file.exists("Regions_boundaries.csv")) {
  rb <- read.csv("Regions_boundaries.csv")
  cat("Regions_boundaries rows:", nrow(rb), "\n")
  cat("Max field_1 index:", max(rb$field_1, na.rm=TRUE), "\n")
  
  gridVec <- rb$field_1[rb$Value>1]
  cat("Max index in gridVec (>1):", max(gridVec, na.rm=TRUE), "\n")
} else {
  cat("Regions_boundaries.csv missing.\n")
}
