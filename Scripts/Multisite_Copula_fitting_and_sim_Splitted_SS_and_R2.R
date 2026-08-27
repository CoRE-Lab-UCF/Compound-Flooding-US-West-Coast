setwd("")
library(VineCopula)
library(VGAM)
library(tweedie)
library(statmod)
library(truncnorm)
library(MASS)
library(VGAM)
library(tweedie)
library(statmod)
library(truncnorm)
library(RColorBrewer)
library(MASS)
library(ks)
library(R.matlab)
library(here)
library(ismev)
library(dplyr)
library(texmex)
library(MultiHazard)
library(parallel)
library(doParallel)
library(foreach)
source("Best_Copula.R")


nSites<-77


# loading the ful timeseries data
Data <- readMat("../Creating_WLONC/Time_series_data.mat")

AR_data <- readMat("../Stratification/AR_data.mat")
AR_events <- lapply(AR_data$AR.data, `[[`, 1)
AR_events <- lapply(AR_events, function(x){colnames(x) <- c("Time_RF", "RF","Time_Hmo", "Hmo")
  x
})


non_AR_data <- readMat("../Stratification/non_AR_data.mat")
non_AR_events <- lapply(non_AR_data$non.AR.data, `[[`, 1)
non_AR_events <- lapply(non_AR_events, function(x){colnames(x) <- c("Time_RF", "RF","Time_Hmo", "Hmo")
  x
})


#source("Migpd_Fit.R")
source("Design_Event_2D_Multi_Pop_splitt.R")
source("fit_gpd_cdf.R")
source("inv_gpd_cdf.R")
source("GPD_Threshold_Tails.R")


n_ensembles <- 1000


# pre-allocate lists outside your loop
AEP_P1_list       <- vector("list", length = nSites)  # each: [nGrid x 2]
AEP_P2_list       <- vector("list", length = nSites)  # each: [nGrid x 2]
Grid_list         <- vector("list", length = nSites)  # each: [nGrid x 2]
MostLikely_list   <- vector("list", length = nSites)  # each: list of 1x2 data.frames per RP
Ensemble_list     <- vector("list", length = nSites)  # each: list of 1x2 data.frames per RP
All_data_list     <- vector("list", length = nSites)  # each: list of 1x2 data.frames per RP

thres_AR_Hmo_all <- vector("list", length = nSites)  # each: list of 1x2 data.frames per RP
thres_AR_RF_all <- vector("list", length = nSites)  # each: list of 1x2 data.frames per RP
thres_non_AR_Hmo_all <- vector("list", length = nSites)  # each: list of 1x2 data.frames per RP
thres_non_AR_RF_all <- vector("list", length = nSites)  # each: list of 1x2 data.frames per RP



for (i in 1:77) { 

WL_above_normal <- Data$Time.series.data[[i]][[1]][, 5] + Data$Time.series.data[[i]][[1]][, 6]
  
full_timeseries <- data.frame(
  Hmo  = WL_above_normal,
  RF = Data$Time.series.data[[i]][[1]][, 2]
)
  
  
# 24-hr rolling accumulation (assumes hourly data, no missing hours)
k <- 24
rf <- full_timeseries$RF
rf24 <- rep(NA_real_, length(rf))
rf24[k:length(rf)] <- vapply(k:length(rf), function(j) sum(rf[(j-k+1):j], na.rm = TRUE), 0.0)

full_timeseries$RF <- rf24
full_timeseries$RF[is.na(full_timeseries$RF)] <- 0



temp_AR_events_data     <- AR_events[[i]]
temp_non_AR_events_data  <- non_AR_events[[i]]


# adding smaller values for zeros
eps <- 1e-8

temp_AR_events_data[temp_AR_events_data == 0]           <- eps
temp_non_AR_events_data[temp_non_AR_events_data == 0]   <- eps


#


# finding number of years
# MATLAB datenum -> Date
toDate <- function(x) as.Date(x, origin = "0000-01-01")
yrs_AR    <- diff(range(toDate(temp_AR_events_data[,  "Time_Hmo"]), na.rm = TRUE))
yrs_nonAR <- diff(range(toDate(temp_non_AR_events_data[, "Time_Hmo"]), na.rm = TRUE))

years_available <- ceiling(as.numeric(max(c(yrs_AR, yrs_nonAR))) / 365.25)


#thres_AR_Hmo <- GPD_Threshold_Solari_mod(Event=temp_AR_events_data,Data=full_timeseries,Min_Quantile = 0.9,var_select="Hmo")
#thres_AR_RF <- GPD_Threshold_Solari_mod(Event=temp_AR_events_data,Data=full_timeseries,Min_Quantile = 0.9,var_select="RF")
#thres_non_AR_Hmo <- GPD_Threshold_Solari_mod(Event=temp_non_AR_events_data,Data=full_timeseries,Min_Quantile = 0.8,var_select="Hmo")
#thres_non_AR_RF <- GPD_Threshold_Solari_mod(Event=temp_non_AR_events_data,Data=full_timeseries,Min_Quantile = 0.8,var_select="RF")

thres_AR_Hmo <- GPD_Threshold_Tails(data_raw=temp_AR_events_data,record_length=years_available,baseline_event_RP= 0.25,threshold_probabilities = seq(0.5, 0.99, by = 0.005),parallel= TRUE,cores= 4,var_select="Hmo",plot= FALSE)
thres_AR_RF <- GPD_Threshold_Tails(data_raw=temp_AR_events_data,record_length=years_available,baseline_event_RP= 0.25,threshold_probabilities = seq(0.5, 0.99, by = 0.005),parallel= TRUE,cores= 4,var_select="RF",plot= FALSE)
thres_non_AR_Hmo <- GPD_Threshold_Tails(data_raw=temp_non_AR_events_data,record_length=years_available,baseline_event_RP= 0.25,threshold_probabilities = seq(0.5, 0.99, by = 0.005),parallel= TRUE,cores= 4,var_select="Hmo",plot= FALSE)
thres_non_AR_RF <- GPD_Threshold_Tails(data_raw=temp_non_AR_events_data,record_length=years_available,baseline_event_RP= 0.25,threshold_probabilities = seq(0.5, 0.99, by = 0.005),parallel= TRUE,cores= 4,var_select="RF",plot= FALSE)

if (length(thres_AR_Hmo$Candidate_Thres$thresh) == 0) {
  thres_AR_Hmo$Candidate_Thres$thresh <- sort(temp_AR_events_data[,4])[2]
}
if (length(thres_AR_RF$Candidate_Thres$thresh) == 0) {
  thres_AR_RF$Candidate_Thres$thresh <- sort(temp_AR_events_data[,2])[2]
}
if (length(thres_non_AR_Hmo$Candidate_Thres$thresh) == 0) {
  thres_non_AR_Hmo$Candidate_Thres$thresh <- sort(temp_non_AR_events_data[,4])[2]
}
if (length(thres_non_AR_RF$Candidate_Thres$thresh) == 0) {
  thres_non_AR_RF$Candidate_Thres$thresh <- sort(temp_non_AR_events_data[,2])[2]
}


AR_GPD_cdf_Hmo<-fit_gpd_cdf(data=temp_AR_events_data[ ,3:4],threshold=thres_AR_Hmo$Candidate_Thres$thresh,n_years=years_available,plot=FALSE)
AR_GPD_cdf_RF<-fit_gpd_cdf(data=temp_AR_events_data[ ,1:2],threshold=thres_AR_RF$Candidate_Thres$thresh,n_years=years_available,plot=FALSE)
non_AR_GPD_cdf_Hmo<-fit_gpd_cdf(data=temp_non_AR_events_data[ ,3:4],threshold=thres_non_AR_Hmo$Candidate_Thres$thresh,n_years=years_available,plot=FALSE)
non_AR_GPD_cdf_RF<-fit_gpd_cdf(data=temp_non_AR_events_data[ ,1:2],threshold=thres_non_AR_RF$Candidate_Thres$thresh,n_years=years_available,plot=FALSE)






#Creating arrays for copula fitting (only the value)
AR_data_for_cop <- data.frame(temp_AR_events_data[,4], temp_AR_events_data[,2])
colnames(AR_data_for_cop)<-c("Hmo","RF")


#Creating arrays for copula fitting (only the value)
non_AR_data_for_cop <- data.frame(temp_non_AR_events_data[,4], temp_non_AR_events_data[,2])
colnames(non_AR_data_for_cop)<-c("Hmo","RF")



# finding the best copula
Copulas <- Best_Copula(Con_var_1= AR_data_for_cop, Con_var_2=non_AR_data_for_cop)
AR_copula_Family <- Copulas$copula_Var1_Family
non_AR_copula_Family <- Copulas$copula_Var2_Family


# --- 2. Run the function ---
out <- Design_Event_2D_Multi_Pop_splitt(
  Data               = full_timeseries,          # full hourly time series, cols: Rainfall, OsWL
  Data_Con1          = AR_data_for_cop,        # conditioned on Rainfall, population 1
  Data_Con4          = non_AR_data_for_cop,        # conditioned on OsWL, population 2
  Copula_Family_pop1 = AR_copula_Family,                # Clayton copula for pop1; NA = auto-select
  Copula_Family_pop2 = non_AR_copula_Family,                # Gumbel copula for pop2
  Con1               = "Rainfall",
  Con2               = "Hmo",
  GPD_P1_var1        = AR_GPD_cdf_Hmo,
  GPD_P1_var2        = AR_GPD_cdf_RF,
  GPD_P2_var1        = non_AR_GPD_cdf_Hmo,
  GPD_P2_var2        = non_AR_GPD_cdf_RF,
  Rate_Con1          = NA,               # auto-computed from data length and mu
  Rate_Con4          = NA,
  # AEP(s)
  RP = c(0.5, 0.2, 0.1, 0.02, 0.01),
  
  # (optional) sampling + plotting controls
  mu = 365.25*24,
  Isoline_Probs = "Sample",
  N = 10^4,
  Sim_Max = 10,
  y_lim_max=200,
  x_lim_max=20,
  N_Ensemble = n_ensembles,
  
  # defining the grid (two times the max values will be used)
  Grid_x_min = 0,
  Grid_x_max = 20,
  Grid_y_min = 0,
  Grid_y_max = 300,
  Grid_x_interval = 0.05,
  Grid_y_interval = 0.5
)



AEP_P1_list[[i]]     <- as.matrix(out$AEP_P1)
AEP_P2_list[[i]]     <- as.matrix(out$AEP_P2)
MostLikely_list[[i]] <- out$MostLikelyEvent
Ensemble_list[[i]] <- out$Ensemble
All_data_list[[i]] <-out

# Saving the thresholds
thres_AR_Hmo_all[i] <- thres_AR_Hmo$Candidate_Thres
thres_AR_RF_all[i] <- thres_AR_RF$Candidate_Thres
thres_non_AR_Hmo_all[i] <- thres_non_AR_Hmo$Candidate_Thres
thres_non_AR_RF_all[i] <- thres_non_AR_RF$Candidate_Thres

Grid_list[[i]]<- as.matrix(out$Grid)
cat("finished for site", i)
}



saveRDS(All_data_list, file = "All_simulated_data.rds")

# dimensions

nGrid  <- nrow(Grid_list[[1]])

# (A) Grid as a numeric array [nGrid x 2 x nSites]
Grid_arr <- array(NA_real_, dim = c(nGrid, 2, nSites))
for (ii in seq_len(nSites)) {
  Grid_arr[ , , ii] <- Grid_list[[ii]]
}

# (B) AEP as numeric arrays [nGrid x 2 x nSites]
AEP_P1_arr <- array(NA_real_, dim = c(nGrid, 2, nSites))
AEP_P2_arr <- array(NA_real_, dim = c(nGrid, 2, nSites))
for (ii in seq_len(nSites)) {
  AEP_P1_arr[ , , ii] <- AEP_P1_list[[ii]]
  AEP_P2_arr[ , , ii] <- AEP_P2_list[[ii]]
}

# (C) MostLikelyEvent as numeric array [nRP x 2 x nSites] + RP vector
# Prefer: you know RP used in the run
RP_vec <- c(2,5, 10, 50, 100)   # <-- change if you used different RPs
nRP <- length(RP_vec)

MLE_arr <- array(NA_real_, dim = c(nRP, 2, nSites)) # [RP x (Hmo,RF) x site]

for (ii in seq_len(nSites)) {
  # if your MostLikelyEvent list is named "10","50","100" use by name;
  # if not, fall back to position.
  ml <- MostLikely_list[[ii]]
  
  for (kk in seq_len(nRP)) {
    key <- as.character(RP_vec[kk])
    
    df <- NULL
    if (!is.null(names(ml)) && key %in% names(ml)) {
      df <- ml[[key]]        # by name
    } else if (length(ml) >= kk) {
      df <- ml[[kk]]         # by position
    }
    
    if (!is.null(df) && is.data.frame(df) && nrow(df) >= 1) {
      # expects columns Hmo and RF
      MLE_arr[kk, , ii] <- as.numeric(df[1, c("Hmo","RF")])
    }
  }
}



# For ensemble of data
nRP <- length(RP_vec)

nPts <- n_ensembles  # number of ensembles

# [RP x point x (Hmo,RF) x site]
ENS_arr <- array(NA_real_, dim = c(nRP, nPts, 2, nSites),dimnames = list(RP   = as.character(RP_vec),pt   = seq_len(nPts),var  = c("Hmo","RF"),site = NULL))




for (ii in seq_len(nSites)) {
  
  ens <- Ensemble_list[[ii]]   # this is out$Ensemble for site ii
  
  for (kk in seq_len(nRP)) {
    key <- as.character(RP_vec[kk])
    
    df <- NULL
    if (!is.null(names(ens)) && key %in% names(ens)) {
      df <- ens[[key]]         # by name (your case)
    } else if (length(ens) >= kk) {
      df <- ens[[kk]]          # fallback by position
    }
    
    if (!is.null(df) && is.data.frame(df)) {
      
      # enforce expected columns + numeric
      tmp <- df[, c("Hmo","RF"), drop = FALSE]
      tmp$Hmo <- as.numeric(tmp$Hmo)
      tmp$RF  <- as.numeric(tmp$RF)
      
      if (nrow(tmp) != nPts) {
        warning(sprintf("Site %d RP %s has %d points (expected %d). Saving min().",
                        ii, key, nrow(tmp), nPts))
      }
      
      m <- min(nPts, nrow(tmp))
      ENS_arr[kk, seq_len(m), , ii] <- as.matrix(tmp[seq_len(m), c("Hmo","RF")])
    }
  }
}







# ---- write .mat ----
# Use v7.3 if objects are large (it usually helps for big arrays)
writeMat("AEP_Grid_all.mat",
  AEP_P1 = AEP_P1_arr,
  AEP_P2 = AEP_P2_arr,
  Grid   = Grid_arr,
  MLE    = MLE_arr,
  RP     = RP_vec,
  Ensemble = ENS_arr,
  thres_AR_Hmo_all=thres_AR_Hmo_all,
  thres_AR_RF_all=thres_AR_RF_all,
  thres_non_AR_Hmo_all=thres_non_AR_Hmo_all,
  thres_non_AR_RF_all=thres_non_AR_RF_all
)
