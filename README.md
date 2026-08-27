# Atmospheric River Compound Flood Analysis – U.S. West Coast

## Overview

This repository contains the workflow used to identify, stratify, and model compound rainfall–water-level extreme events along the U.S. West Coast. The methodology combines Peak-Over-Threshold (POT) event extraction, Atmospheric River (AR) stratification, Generalized Pareto Distribution (GPD) modeling, copula-based dependence modeling, and multivariate design event generation.

The workflow separates compound events into AR and non-AR populations and develops population-specific statistical models to estimate design events across multiple annual exceedance probabilities (AEPs).

---

## Workflow

### 1. POT Event Extraction

The first step identifies rainfall POT, nTWL POT, and compound rainfall–nTWL POT extreme events from hourly time series. For the main analysis, compound rainfall–nTWL POT extreme events are used

#### Main Script

- `Main.m`

#### Functions

- `find_POT_compound_RF_WLNC.m`
- `find_POT_oneway_RF.m`
- `find_POT_oneway_NTTWL.m`

#### Procedure

1. Calculate daily rainfall totals.
2. Estimate rainfall and water-level thresholds using percentile-based methods.
3. Estimate Stable Extremal Dependence (SED) parameters.
4. Identify:
   - Compound rainfall–water-level events
   - Rainfall-only extreme events
   - Water-level-only extreme events
5. Decluster events using a moving temporal window.
6. Save extracted POT events for all sites.

#### Outputs

- `POT_both_extreme.mat`
- `POT_RF_only_extreme.mat`
- `POT_NTTWL_only_extreme.mat`

---

### 2. Atmospheric River Stratification

Compound events are classified as either Atmospheric River (AR) or non-AR events using the ERA5 Global Atmospheric River Catalog.

#### Main Script

- `Main.m`

#### Procedure

1. Load compound POT events.
2. Load ERA5 Atmospheric River catalog.
3. Match rainfall peak times with AR occurrences.
4. Apply temporal and spatial search criteria around each watershed centroid.
5. Classify each event as:
   - AR event
   - Non-AR event

#### Outputs

- `AR_data.mat`
- `non_AR_data.mat`

---

### 3. Copula Selection

Dependence structures between rainfall and water level are modeled separately for AR and non-AR populations.

#### Function

- `Best_Copula.R`

#### Procedure

1. Convert observations to pseudo-observations.
2. Compute Kendall's Tau.
3. Fit candidate copula families using the VineCopula package.
4. Select the optimal copula using AIC.
5. Return the best copula family for each population.


### 4. Marginal Distribution Modeling

Marginal distributions are modeled independently for rainfall and water level using Generalized Pareto Distributions (GPD).

#### Functions

- `GPD_Threshold_Tails.R`
- `fit_gpd_cdf.R`
- `inv_gpd_cdf.R`

#### Procedure

1. Estimate optimal thresholds.
2. Fit GPD tail models.
3. Construct cumulative distribution functions.
4. Generate inverse CDFs for simulation.

---

### 5. Design Event Generation

Design events are generated separately for AR and non-AR populations and combined through a mixed-population framework.(Maduwantha et al., 2024)

#### Main Script

- `Multisite_Copula_fitting_and_sim_Splitted_SS_and_R2.R`

#### Supporting Function

- `Design_Event_2D_Multi_Pop_splitt.R`

#### Procedure

For each watershed:

1. Load AR and non-AR event populations.
2. Fit marginal distributions.
3. Select optimal copula families.
4. Generate AEP surfaces.
5. Estimate most-likely design events.
6. Generate Monte Carlo ensembles of design events.



## Data Requirements

- Rainfall
- Storm surge
- 2% Runup (Significant wave height, wave period, beach slope)
- Atmospheric River Catalog (Spatial footprints)



## Outputs

### Event Databases

- `POT_both_extreme.mat`
- `POT_RF_only_extreme.mat`
- `POT_NTTWL_only_extreme.mat`
- `AR_data.mat`
- `non_AR_data.mat`
- Selected copula families
- GPD thresholds
- GPD parameters
- AEP surfaces
- Most likely design events
- Ensemble simulations

Saved in:

- `AEP_Grid_all.mat`
- `All_simulated_data.rds`


## Software Requirements

### MATLAB

Required toolboxes:

- Statistics and Machine Learning Toolbox
- Mapping Toolbox
- NetCDF support

### R

Required packages:

```r
VineCopula
VGAM
tweedie
statmod
truncnorm
MASS
RColorBrewer
ks
R.matlab
here
ismev
dplyr
texmex
MultiHazard
parallel
doParallel
foreach
```

---

## Author

Pravin Maduwantha
