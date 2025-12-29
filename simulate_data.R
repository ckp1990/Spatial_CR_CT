# simulate_data.R
# This script generates dummy data for the MLDmodel11full.R script
# Run this to create the necessary CSV files for testing/verification.

set.seed(123)

# Configuration
n_traps <- 10
n_occasions <- 5
n_animals <- 20
n_captures <- 50

# 1. Create Statespace (new_mask_file.csv)
# Grid of 20x20
x <- seq(0, 2000, length.out = 20)
y <- seq(0, 2000, length.out = 20)
grid <- expand.grid(X = x, Y = y)
grid$Habitat <- 1 # All good habitat
# Add a few bad habitat pixels
grid$Habitat[sample(1:nrow(grid), 10)] <- 0
write.csv(grid, "new_mask_file.csv", row.names = FALSE)
cat("Created new_mask_file.csv\n")

# 2. Create Traps (Traps.csv)
# Traps need to have: TrapID, X, Y, and then one column per occasion indicating if active
# 10 Random traps within the grid
n_traps <- 10
traps <- data.frame(
  TrapID = 1:n_traps,
  X = sample(x, n_traps, replace = TRUE),
  Y = sample(y, n_traps, replace = TRUE)
)

# Add binary columns for active/inactive for each occasion
# MASK<-as.matrix(traps[,4:ncol(traps)]) implies columns 4 onwards are the mask
for (i in 1:n_occasions) {
  traps[[paste0("Active", i)]] <- 1
}

# Add binary columns for active/inactive (assuming 1 session for now, or just 'Active')
# The code `MASK<-as.matrix(traps[,4:ncol(traps)])` implies columns 4 onwards are the mask
traps$Active <- 1
write.csv(traps, "Traps.csv", row.names = FALSE)
cat("Created Traps.csv\n")

# 3. Create Captures (Capture.csv)
captures <- data.frame(
  LOC_ID = sample(1:n_traps, n_captures, replace = TRUE), # Trap ID
  ANIMAL_ID = sample(1:n_animals, n_captures, replace = TRUE),
  SO = sample(1:n_occasions, n_captures, replace = TRUE) # Sampling Occasion
)
write.csv(captures, "Capture.csv", row.names = FALSE)
cat("Created Capture.csv\n")

# 4. Create Sex Data (Sex.csv)
# Sex for each of the 20 animals
sex_data <- data.frame(
  AnimalID = 1:n_animals,
  Sex = sample(c(0, 1), n_animals, replace = TRUE)
)
# The code `Xsex <- sex[,2]` expects the second column to be the sex value
write.csv(sex_data, "Sex.csv", row.names = FALSE)
cat("Created Sex.csv\n")

cat("All dummy data files created successfully.\n")
