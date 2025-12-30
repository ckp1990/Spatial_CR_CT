# few suffixes
Output_number <- "131"
result_id <- "BD2005_gender_sps"
# load up the library
library(coda)
library(mcmcse)
main <- getwd() # this will automatically take the dir where the code stored
# create a dir for storing the result.
result_dir<-function(op_dir){
  # Robust check for "output" in path
  if(!grepl("output", op_dir)){
    if(dir.exists("result")==F){
      dir.create("result")
    }
    return(file.path(op_dir, "result"))
  }
  
  op_split<-strsplit(op_dir,"output")
  setwd(op_split[[1]][1])
  if(dir.exists("result")==F){
    dir.create("result")
  }
  setwd("result/")
  
  # Safety check for accessing second element
  if(length(op_split[[1]]) > 1){
    split_part <- strsplit(op_split[[1]][2],"/")[[1]]
    if(length(split_part) > 1){
       folder_name <- split_part[2]
       if(dir.exists(folder_name)==F){
          dir.create(folder_name)
       }
       file_name<-folder_name
       setwd(file_name)
    }
  }
  return(getwd())
}
result <- result_dir(main)
# back ot op folder
setwd(main)

### Obtain all the MCMC histories. These directory structures obtained dynamically ###

# Robust directory detection looking for Checkpoint folders (usually CH + number)
# Adapting to the patterns seen: (NLD)*(CH1), (MLD)*(CH3), etc.
# We search for directories containing "CH" followed by digits.
chain_dirs <- grep("CH[0-9]+", list.files(), value = TRUE)
chain_dirs <- paste0(chain_dirs, "/")

# Function to load MCMC history from a directory
load_chain <- function(dir_path, main_path) {
  setwd(dir_path)
  # Find the MCMC log file (pattern "MLD")
  file_name <- list.files(pattern = "MLD")
  if (length(file_name) == 0) {
    setwd(main_path)
    warning(paste("No MCMC file found in", dir_path))
    return(NULL)
  }
  hist <- read.csv(file_name[1]) # Take the first match if multiple
  setwd(main_path)
  return(as.mcmc(hist))
}

# Load all chains
raw_chains <- lapply(chain_dirs, load_chain, main_path = main)
# Remove NULLs if any directories failed
raw_chains <- raw_chains[!sapply(raw_chains, is.null)]

if (length(raw_chains) == 0) {
  stop("No valid MCMC chains were loaded!")
}

#### Gelman-Rubin diagnostics ####

## Determine robustness for windowing
# Find the minimum number of iterations across all chains to ensure consistency
min_iter <- min(sapply(raw_chains, nrow))

## Remove beta.behave column (1), X column (maybe 1?), Density (7?)
# Original code removed c(-1, -7). Assuming structure is consistent.
# Set start and end.
start <- 1
end <- min_iter # Use the full available length of the shortest chain
# OR allow user override if they want a specific burn-in, e.g.,
# start <- floor(min_iter * 0.1) # 10% burn-in auto

processed_chains <- lapply(raw_chains, function(ch) {
  # Safety check for column indices
  if (ncol(ch) >= 7) {
    subset_ch <- ch[, c(-1, -7)]
  } else {
    subset_ch <- ch[, -1] # Fallback if fewer columns
  }
  window(subset_ch, start = start, end = end)
})

### Combine chain outputs ###
# Convert list to mcmc.list for coda functions
chainList <- as.mcmc.list(processed_chains)
combinedHist <- do.call(rbind, processed_chains)

## MCMC diagnostics ##
## Multi-chain convergence check using Gelman-Rubin diagnostic
gelmandiag <- gelman.diag(chainList, confidence = FALSE, transform = FALSE, autoburnin = FALSE, multivariate = FALSE)
# capture.output(as.data.frame(gelmandiag[[1]]),file = paste(result,"gelman_30K.txt",sep = "/"))
parameter <- rownames(as.data.frame(gelmandiag[[1]]))
gelman_estimate <- gelmandiag[[1]][, 1]
Gelman_output <- data.frame(Output_number = rep(Output_number, length(parameter)), parameter, gelman_estimate, row.names = NULL)
write.csv(Gelman_output, file = paste(result, "Gelman_output.csv", sep = "/"), row.names = F)

estvssamp(combinedHist)
# capture.output(as.data.frame(mcse.mat(combinedHist)),file="estimate_30k.txt")
## Single chain convergence check using Geweke diagnostic (optional)
# gewekediag <- geweke.diag(histCH1mcmc)

#### Report MCMC diagnostic results.
# For Geweke we want the magnitude (-ve or +ve) for each parameter
# to be less than 1.64.
# For Gelman-Rubin we want Potential Shrink Reduction Factor to b
# less than 1.2 for each parameter.

# gewekediag (optional)

### Summary results. Look for how different median is to the mean. This indicates nature of the posterior distribution. Ideally we would like them to be nearly the same. But OK otherwise too.
mean.model1 <- apply(combinedHist, 2, mean)
sd.model1 <- apply(combinedHist, 2, sd)
mean.model1
sd.model1
## msmc quanlite
qut_2.5 <- mcse.q.mat(combinedHist, c(2.5 / 100))
qut_50 <- mcse.q.mat(combinedHist, c(50 / 100))
qut_97.5 <- mcse.q.mat(combinedHist, c(97.5 / 100))

# 3creating a dataframe later to be saved
result_df <- data.frame(mean.model1, sd.model1)
result_df <- cbind(result_df, se = as.data.frame(mcse.mat(combinedHist))[, 2])
library(dplyr)
result_df <- cbind(result_df, qut_2.5, qut_50, qut_97.5)
names(result_df) <- c(
  "mean", "sd", "mcmcse", "quant_2.5", "quant_2_5_mcmcse", "quant_50", "quant_50_mcmcse",
  "quant_97.5", "quant_97.5_mcmcse"
)
result_df <- cbind(output_number = Output_number, result_id, output_parameter = rownames(result_df), result_df)
write.csv(result_df, paste(result, "Estimate_output.csv", sep = "/"), row.names = F)


source("Parameter_for_traceplot.R")
source("traceplot.R")

trace_plot_parameters <- c("sigma", "sigma2", "lam0", "psi", "psi.sex", "Nsuper", "D.adj") # trace_plot_par_generator(result_df)
for (n in trace_plot_parameters) {
  trace_plot_ckp(chainList, n, result)
}
setwd(main)
gc()
### Applying pixel area correction for LION analysis (this was done to reflect a small discrepancy in GIS pixel area calculations in our specific case) ###

# newDenAdj <- combinedHist[,14]
# mean(newDenAdj)
# sd(newDenAdj)

## Highest posterior density intervals for one of the chains # This piece of code is taken from SPACECAP version 1.1 (Gopalaswamy et al. 2015) ##
## HPDinterval(histCH1mcmc) ##

#### Goodness-of-fit statistics ####

### Obtain all the gof statistics ###
### Obtain all the gof statistics loop ###
gdata_list <- list()
gnew_list <- list()

for (dir_path in chain_dirs) {
  setwd(dir_path)

  # Read gofdata
  gdata_files <- list.files(pattern = "gofdata")
  if (length(gdata_files) > 0) {
    gdata_list[[length(gdata_list) + 1]] <- read.csv(gdata_files[1])
  }

  # Read gofnew
  gnew_files <- list.files(pattern = "gofnew")
  if (length(gnew_files) > 0) {
    gnew_list[[length(gnew_list) + 1]] <- read.csv(gnew_files[1])
  }

  setwd(main)
}

### Combine gdata and gnew ###
if (length(gdata_list) > 0 && length(gnew_list) > 0) {
  gdatacombined <- do.call(rbind, gdata_list)
  gnewcombined <- do.call(rbind, gnew_list)
} else {
  stop("GoF files (gofdata/gofnew) could not be loaded from chain directories.")
}

## Bayesian p-value calculation ##
BayesPval <- mean(gdatacombined[, 2] > gnewcombined[, 2])
# capture.output(BayesPval,file = paste(result,"pval.txt",sep = "/"))
P_value <- data.frame(output_number = Output_number, result_id, p_value = BayesPval)
write.csv(P_value, file = , paste(result, "P_value.csv", sep = "/"), row.names = F)

### Generate pixel-specific density estimates ###

## Obtain all the necessary files pixel-specific density estimates ##
## Obtain all the necessary files pixel-specific density estimates ##
AcCentres_list <- list()
RealIndividuals_list <- list()
SexIndividuals_list <- list()

for (dir_path in chain_dirs) {
  setwd(dir_path)

  # Activity Centres
  ac_files <- list.files(pattern = "AcCentres")
  if (length(ac_files) > 0) AcCentres_list[[length(AcCentres_list) + 1]] <- read.csv(ac_files[1])

  # Real Individuals
  real_files <- list.files(pattern = "RealIndividuals")
  if (length(real_files) > 0) RealIndividuals_list[[length(RealIndividuals_list) + 1]] <- read.csv(real_files[1])

  # Sex Individuals
  sex_files <- list.files(pattern = "SexIndividuals")
  if (length(sex_files) > 0) SexIndividuals_list[[length(SexIndividuals_list) + 1]] <- read.csv(sex_files[1])

  setwd(main)
}

## Combine into combined history ##
AcCentresCombined <- do.call(rbind, AcCentres_list)
RealIndividualsCombined <- do.call(rbind, RealIndividuals_list)
SeXIndividualcombined <- do.call(rbind, SexIndividuals_list)

# male pixel
## since the SexIndividualcombined is 0 and 1 when 0 is female, so multiple with real will give
## of pixel density male tigers.

male_individuals <- RealIndividualsCombined * SeXIndividualcombined

female_individuals <- ((SeXIndividualcombined - 1) * (-1)) * RealIndividualsCombined

## Obtain the un-scaled state space (any one chain is sufficient as it comes from input data) ##
## Obtain the un-scaled state space (using the first available chain) ##
if (length(chain_dirs) > 0) {
  setwd(dir = chain_dirs[1])
  # Using first chain directory dynamically
  ss_files <- list.files(pattern = "SSunscaled")
  SSunscaledCH8 <- read.csv(ss_files[1])
} else {
  stop("No chain directories found to read SSunscaled")
}
nG <- nrow(SSunscaledCH8)

# Set pixel ID of home range centers for phantom animals to zero
indlocs <- AcCentresCombined * male_individuals
indlocnum <- data.matrix(indlocs)
indlocnum <- indlocnum[, -1] ## deleting the row name

# Count the proportion of times each pixel was a home range centre,
#   convert to animals per sq km (here 0.3975 was input data for Lion analysis - so change accordingly)

densVec <- tabulate(indlocnum[, ], nbins = nG) / nrow(indlocs) / 0.336 # make sure first row now used


setwd(main)
maraLionSS <- read.csv("BD_for_Abundance.csv")
pixelDensity <- maraLionSS
pixelDensity$`Pixel Density` <- maraLionSS[, 3]
target_indices <- which(maraLionSS[, 3] > 0)
if(length(target_indices) == length(densVec)){
  pixelDensity$`Pixel Density`[target_indices] <- densVec
} else {
  warning("Dimension mismatch between density vector and map. Skipping pixel density assignment.")
}
names(pixelDensity) <- c("x_coordinate", "y_coordinate", "habitat", "value", "pixel_density")
pixelDensity <- cbind(output_number = Output_number, sex = "1", pixelDensity)
pixelDensity_male <- pixelDensity
# This part is to obtain posterior standard deviations on pixel-specific densities #

# Create an abundance matrix of dimension no. of iterations x total number of grid cells #
abundMatrix <- matrix(data = NA, nrow = nrow(indlocs), ncol = nG)

# Fill up this matrix with abundance counts for each iteration #
for (i in 1:nrow(indlocs)) {
  abundVecTemp <- tabulate(indlocnum[i, ], nbins = nG)
  abundMatrix[i, ] <- abundVecTemp
}

# This part is meant to compute abundances for sub-regions #
# Enter the sequence of grid cell numbers for analysis. This will be a selection of numbers between 1 to nG corresponding to which cells are being analysed. For the entire study area this can simply be 1:nG. In the example below, it indicates that only grid cells 1,3,5,7 are chosen for reporting. This will be according to the sub-region chosen #

boundary_point <- read.csv("Regions_boundaries.csv", header = T)
gridVec <- boundary_point$field_1[boundary_point$Value > 1]

# Obtain total abundance counts for grid cells referenced by gridVec for each iteration #
# Obtain total abundance counts for grid cells referenced by gridVec for each iteration #
if(length(gridVec) > 0 && max(gridVec, na.rm=TRUE) <= ncol(abundMatrix)) {
  abundVecTotal <- rowSums(abundMatrix[,gridVec])
  meanAbund <- mean(abundVecTotal)
  sdAbund <- sd(abundVecTotal)
  Park_inside_abundance_M<-data.frame(output_number=Output_number,result_id,sex="1",mean=meanAbund,sd=sdAbund)
} else {
  warning("Grid indices in Regions_boundaries.csv exceed state space size. Skipping Park Inside Abundance (First Block).")
  Park_inside_abundance_M<-data.frame(output_number=Output_number,result_id,sex="1",mean=NA,sd=NA)
}

## female pixel density
female_individuals <- ((SeXIndividualcombined - 1) * (-1)) * RealIndividualsCombined

## Obtain the unscaled statespace (any one chain is sufficient as it comes from input data) ##
## Obtain the unscaled statespace (using the first available chain) ##
if (length(chain_dirs) > 0) {
  setwd(dir = chain_dirs[1])
  ss_files <- list.files(pattern = "SSunscaled")
  SSunscaledCH8 <- read.csv(ss_files[1])
}
nG <- nrow(SSunscaledCH8)

# Set pixel ID of home range centers for phantom animals to zero
indlocs <- AcCentresCombined * female_individuals
indlocnum <- data.matrix(indlocs)
indlocnum <- indlocnum[, -1] ## deleting the row name

# Count the proportion of times each pixel was a home range centre,
#   convert to animals per sq km (here 0.3975 was input data for Lion analysis - so change accordingly)

densVec <- tabulate(indlocnum[, ], nbins = nG) / nrow(indlocs) / 0.336 # make sure first row now used


setwd(main)
maraLionSS <- read.csv("BD_for_Abundance.csv")
pixelDensity <- maraLionSS
pixelDensity$`Pixel Density` <- maraLionSS[, 3]
pixelDensity$`Pixel Density`[maraLionSS[, 3] > 0] <- densVec
names(pixelDensity) <- c("x_coordinate", "y_coordinate", "habitat", "value", "pixel_density")
pixelDensity_female <- cbind(output_number = Output_number, sex = "0", pixelDensity)
# Generate csv file for pixel densities #

nameoffile3 <- paste(result, "Pixel_density_map.csv", sep = "/")
pixelDensity <- bind_rows(pixelDensity_male, pixelDensity_female)
write.csv(pixelDensity, file = nameoffile3, row.names = F)

# This part is to obtain posterior standard deviations on pixel-specific densities #

# Create an abundance matrix of dimension no. of iterations x total number of grid cells #
abundMatrix <- matrix(data = NA, nrow = nrow(indlocs), ncol = nG)

# Fill up this matrix with abundance counts for each iteration #
for (i in 1:nrow(indlocs)) {
  abundVecTemp <- tabulate(indlocnum[i, ], nbins = nG)
  abundMatrix[i, ] <- abundVecTemp
}

# This part is meant to compute abundances for sub-regions #
# Enter the sequence of grid cell numbers for analysis. This will be a selection of numbers between 1 to nG corresponding to which cells are being analysed. For the entire study area this can simply be 1:nG. In the example below, it indicates that only grid cells 1,3,5,7 are chosen for reporting. This will be according to the sub-region chosen #

boundary_point <- read.csv("Regions_boundaries.csv", header = T)
gridVec <- boundary_point$field_1[boundary_point$Value > 1]

# Obtain total abundance counts for grid cells referenced by gridVec for each iteration #
# Obtain total abundance counts for grid cells referenced by gridVec for each iteration #
if(length(gridVec) > 0 && max(gridVec, na.rm=TRUE) <= ncol(abundMatrix)) {
  abundVecTotal <- rowSums(abundMatrix[,gridVec])
  meanAbund <- mean(abundVecTotal)
  sdAbund <- sd(abundVecTotal)
  Park_inside_abundance_F<-data.frame(output_number=Output_number,result_id,sex="0",mean=meanAbund,sd=sdAbund)
} else {
  warning("Grid indices in Regions_boundaries.csv exceed state space size. Skipping Park Inside Abundance (Female).")
  Park_inside_abundance_F<-data.frame(output_number=Output_number,result_id,sex="0",mean=NA,sd=NA)
}
#capture.output(c(meanAbund,sdAbund),file = paste(result,"Park_inside_abundance_female.txt",sep = "/"))
Park_inside_abundance<-bind_rows(Park_inside_abundance_M,Park_inside_abundance_F)
write.csv(Park_inside_abundance,file = paste(result,"Park_inside_abundance.csv",sep = "/"),row.names = F)


## Generate pair-wise plots. This will be useful for assessing estimation covariances and parameter redundancies (if any) owing to poor sample sizes ##
jpeg(file = paste(result, "Pair_plot.jpeg", sep = "/"), width = 800, height = 800, units = "px")
library(PerformanceAnalytics)
chart.Correlation(combinedHist[, c("sigma", "sigma2", "lam0", "psi", "psi.sex", "D.adj")], histogram = TRUE, pch = 19)
dev.off()
library(correlation)
dat <- as.data.frame(combinedHist[, c("sigma", "sigma2", "lam0", "psi", "psi.sex", "D.adj")])

co_plot_table <- correlation(dat,
  include_factors = TRUE, method = "pearson"
)
co_plot_table <- co_plot_table[, c(1, 2, 3)]
co_plot_table <- cbind(output_number = Output_number, co_plot_table)
write.csv(co_plot_table, file = paste(result, "co_plot_table.csv", sep = "/"), row.names = F)
## outside abundance#########################################################
###############################################################################
# male pixel
## since the SexIndividualcombined is 0 and 1 when 0 is female, so multiple with real will give
## of pixel density male tigers.

male_individuals <- RealIndividualsCombined * SeXIndividualcombined

female_individuals <- ((SeXIndividualcombined - 1) * (-1)) * RealIndividualsCombined

## Obtain the un-scaled state space (any one chain is sufficient as it comes from input data) ##
## Obtain the un-scaled state space ##
if (length(chain_dirs) > 0) {
  setwd(dir = chain_dirs[1])
  ss_files <- list.files(pattern = "SSunscaled")
  SSunscaledCH8 <- read.csv(ss_files[1])
}
nG <- nrow(SSunscaledCH8)

# Set pixel ID of home range centers for phantom animals to zero
indlocs <- AcCentresCombined * male_individuals
indlocnum <- data.matrix(indlocs)
indlocnum <- indlocnum[, -1] ## deleting the row name

# Count the proportion of times each pixel was a home range centre,
#   convert to animals per sq km (here 0.3975 was input data for Lion analysis - so change accordingly)

densVec <- tabulate(indlocnum[, ], nbins = nG) / nrow(indlocs) / 0.336 # make sure first row now used


setwd(main)
maraLionSS <- read.csv("BD_for_Abundance.csv")
pixelDensity <- maraLionSS
pixelDensity$`Pixel Density` <- maraLionSS[, 3]
target_indices <- which(maraLionSS[, 3] > 0)
if(length(target_indices) == length(densVec)){
  pixelDensity$`Pixel Density`[target_indices] <- densVec
} else {
  warning("Dimension mismatch between density vector and map. Skipping pixel density assignment (Block 2).")
}
names(pixelDensity) <- c("x_coordinate", "y_coordinate", "habitat", "value", "pixel_density")
pixelDensity <- cbind(output_number = Output_number, sex = "1", pixelDensity)
pixelDensity_male <- pixelDensity
# This part is to obtain posterior standard deviations on pixel-specific densities #

# Create an abundance matrix of dimension no. of iterations x total number of grid cells #
abundMatrix <- matrix(data = NA, nrow = nrow(indlocs), ncol = nG)

# Fill up this matrix with abundance counts for each iteration #
for (i in 1:nrow(indlocs)) {
  abundVecTemp <- tabulate(indlocnum[i, ], nbins = nG)
  abundMatrix[i, ] <- abundVecTemp
}

# This part is meant to compute abundances for sub-regions #
# Enter the sequence of grid cell numbers for analysis. This will be a selection of numbers between 1 to nG corresponding to which cells are being analysed. For the entire study area this can simply be 1:nG. In the example below, it indicates that only grid cells 1,3,5,7 are chosen for reporting. This will be according to the sub-region chosen #

boundary_point <- read.csv("Regions_boundaries.csv", header = T)
gridVec <- boundary_point$field_1[boundary_point$Value == 1]

# Obtain total abundance counts for grid cells referenced by gridVec for each iteration #
# Obtain total abundance counts for grid cells referenced by gridVec for each iteration #
if(length(gridVec) > 0 && max(gridVec, na.rm=TRUE) <= ncol(abundMatrix)) {
  abundVecTotal <- rowSums(abundMatrix[,gridVec])
  meanAbund <- mean(abundVecTotal)
  sdAbund <- sd(abundVecTotal)
  Park_inside_abundance_M<-data.frame(output_number=Output_number,result_id,sex="1",mean=meanAbund,sd=sdAbund)
} else {
  warning("Grid indices in Regions_boundaries.csv exceed state space size. Skipping Park Inside Abundance (Male - Outside).")
  Park_inside_abundance_M<-data.frame(output_number=Output_number,result_id,sex="1",mean=NA,sd=NA)
}

## female pixel density
female_individuals <- ((SeXIndividualcombined - 1) * (-1)) * RealIndividualsCombined

## Obtain the unscaled statespace (any one chain is sufficient as it comes from input data) ##
## Obtain the unscaled statespace ##
if (length(chain_dirs) > 0) {
  setwd(dir = chain_dirs[1])
  ss_files <- list.files(pattern = "SSunscaled")
  SSunscaledCH8 <- read.csv(ss_files[1])
}
nG <- nrow(SSunscaledCH8)

# Set pixel ID of home range centers for phantom animals to zero
indlocs <- AcCentresCombined * female_individuals
indlocnum <- data.matrix(indlocs)
indlocnum <- indlocnum[, -1] ## deleting the row name

# Count the proportion of times each pixel was a home range centre,
#   convert to animals per sq km (here 0.3975 was input data for Lion analysis - so change accordingly)

densVec <- tabulate(indlocnum[, ], nbins = nG) / nrow(indlocs) / 0.336 # make sure first row now used


setwd(main)
maraLionSS <- read.csv("BD_for_Abundance.csv")
pixelDensity <- maraLionSS
pixelDensity$`Pixel Density` <- maraLionSS[, 3]
pixelDensity$`Pixel Density`[maraLionSS[, 3] > 0] <- densVec
names(pixelDensity) <- c("x_coordinate", "y_coordinate", "habitat", "value", "pixel_density")
pixelDensity_female <- cbind(output_number = Output_number, sex = "0", pixelDensity)
# Generate csv file for pixel densities #

# This part is to obtain posterior standard deviations on pixel-specific densities #
# Create an abundance matrix of dimension no. of iterations x total number of grid cells #
abundMatrix <- matrix(data = NA, nrow = nrow(indlocs), ncol = nG)

# Fill up this matrix with abundance counts for each iteration #
for (i in 1:nrow(indlocs)) {
  abundVecTemp <- tabulate(indlocnum[i, ], nbins = nG)
  abundMatrix[i, ] <- abundVecTemp
}

# This part is meant to compute abundances for sub-regions #
# Enter the sequence of grid cell numbers for analysis. This will be a selection of numbers between 1 to nG corresponding to which cells are being analysed. For the entire study area this can simply be 1:nG. In the example below, it indicates that only grid cells 1,3,5,7 are chosen for reporting. This will be according to the sub-region chosen #

boundary_point <- read.csv("Regions_boundaries.csv", header = T)
gridVec <- boundary_point$field_1[boundary_point$Value == 1]

# Obtain total abundance counts for grid cells referenced by gridVec for each iteration #
abundVecTotal <- rowSums(abundMatrix[, gridVec])

# Obtain posterior mean and standard deviations of the sub-region
meanAbund <- mean(abundVecTotal)
sdAbund <- sd(abundVecTotal)
Park_inside_abundance_F <- data.frame(output_number = Output_number, result_id, sex = "0", mean = meanAbund, sd = sdAbund)
# capture.output(c(meanAbund,sdAbund),file = paste(result,"Park_inside_abundance_female.txt",sep = "/"))
Park_inside_abundance <- bind_rows(Park_inside_abundance_M, Park_inside_abundance_F)
write.csv(Park_inside_abundance, file = paste(result, "Park_outside_abundance.csv", sep = "/"), row.names = F)
########################
rm(list = ls())
gc()
