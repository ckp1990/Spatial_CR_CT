# SCR Analysis Workflow Guide

This document provides a step-by-step guide for ecologists to run the Spatial Capture-Recapture (SCR) analysis using the provided R scripts.

## 1. Prerequisites

Before running the analysis, ensure you have R installed. The following R packages are required:

- `coda`
- `mcmcse`
- `doParallel`
- `foreach`
- `PerformanceAnalytics`
- `correlation`
- `dplyr`

You can install them using:
```R
install.packages(c("coda", "mcmcse", "doParallel", "foreach", "PerformanceAnalytics", "correlation", "dplyr"))
```

## 2. Data Preparation

Ensure your input data files are formatted correctly and placed in the project root directory.

### Required Files

1.  **Traps File** (`Traps.csv` by default)
    *   **Structure**: The first column is the Trap ID. Columns 2 and 3 are X and Y coordinates. Subsequent columns represent binary activity status (1 = active, 0 = inactive) for each sampling occasion.
    *   **Example**:
        ```csv
        TrapID,X,Y,Occasion1,Occasion2,Occasion3
        1,500000,200000,1,1,1
        2,500050,200050,1,0,1
        ```

2.  **Captures File** (`Capture.csv` by default)
    *   **Structure**: Must contain at least three columns representing the Individual ID, Sampling Occasion, and Trap Location ID.
    *   **Expected Headers** (based on script logic): `ANIMAL_ID`, `SO` (Sample Occasion), `LOC_ID`.
    *   **Example**:
        ```csv
        ANIMAL_ID,SO,LOC_ID
        Leo1,1,10
        Leo1,2,12
        Leo2,1,5
        ```

3.  **Sex File** (`Sex.csv` by default)
    *   **Structure**: Typically two columns. The script expects the sex information (0 or 1) to be in the second column.
    *   **Example**:
        ```csv
        ID,Sex
        Leo1,1
        Leo2,0
        ```
    *   *Note*: Ensure `1` and `0` typically map to Male/Female consistent with your analysis interpretation (e.g., in this code, often 1=Male, 0=Female).

4.  **Statespace/Mask File** (`new_mask_file.csv` by default)
    *   **Structure**: Defines the study area grid. Columns usually represent X coordinate, Y coordinate, and Habitat suitability (1 = suitable, 0 = unsuitable).
    *   **Example**:
        ```csv
        X,Y,Habitat
        500000,200000,1
        500500,200000,1
        501000,200000,0
        ```

## 3. Running the Analysis

The main control script is `MLDmodel11full.R`.

### Step 3.1: Configuration

Open `MLDmodel11full.R` and modify the configuration section at the top to match your specific analysis needs:

```R
config <- list(
  model_number = "MyAnalysisRun",  # Unique identifier for this run
  n_chains = 3,                    # Number of parallel chains
  n_iter = 120000,                 # Total iterations
  burn_in = 20000,                 # Burn-in period
  thining_rate = 1,                # Thinning rate
  nz = 100,                        # Number of augmented zeros

  # Input filenames
  file_statespace = "new_mask_file.csv",
  file_traps = "Traps.csv",
  file_captures = "Capture.csv",
  file_sex = "Sex.csv",

  # ... other model parameters
)
```

### Step 3.2: Execution

Run the entire script `MLDmodel11full.R` in R or RStudio.
*   The script will verify if all required files exist.
*   It will launch parallel chains using `doParallel`.
*   Results will be saved in a new directory named like `NLDresultsModel_[ModelNumber]_[Timestamp]CH[ChainID]`.

## 4. Post-Processing & Results

After the MCMC chains finish, use `Analysis_code_gender_sps.R` to process the results.

### Step 4.1: Setup

Open `Analysis_code_gender_sps.R`.
*   Ensure the `output_number` or directory logic points to your results. The script dynamically looks for folders with "CH" in their name (e.g., `CH1`, `CH2`).
*   Ensure you are running this from the parent directory containing the `NLDresults...` folders.

### Step 4.2: Execution

Run `Analysis_code_gender_sps.R`. This script will:
1.  **Load Chains**: Automatically detect and load MCMC logs from the result directories.
2.  **Diagnostics**: Calculate Gelman-Rubin statistics to check for convergence.
    *   Output: `Gelman_output.csv` in the `result/` folder.
3.  **Estimates**: Calculate mean, SD, and quantiles for model parameters.
    *   Output: `Estimate_output.csv`.
4.  **Density Mapping**: Generate pixel-specific density maps for males and females.
    *   Output: `Pixel_density_map.csv`.
5.  **Abundance**: Calculate abundance inside specified regions (defined in `Regions_boundaries.csv`).
    *   Output: `Park_inside_abundance.csv` and `Park_outside_abundance.csv`.

## 5. Output Interpretation

Check the `result/` directory (created automatically) for:
*   **Convergence**: Look at `Gelman_output.csv`. Values close to 1 (e.g., < 1.1) indicate good convergence.
*   **Estimates**: `Estimate_output.csv` provides the population parameters (sigma, lambda, density, etc.).
*   **Maps**: Use `Pixel_density_map.csv` in GIS software to visualize spatial density.

---
**Note**: Always ensure your working directory is set to the project root before starting.
