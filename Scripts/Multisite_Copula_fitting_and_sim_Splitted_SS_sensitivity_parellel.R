start_time <- Sys.time()

# Packages
library(VineCopula)
library(VGAM)
library(tweedie)
library(statmod)
library(truncnorm)
library(MASS)
library(ks)
library(R.matlab)
library(ismev)
library(dplyr)
library(texmex)
library(MultiHazard)
library(parallel)
library(doParallel)
library(foreach)

# Analysis functions
source("Best_Copula.R")
source("Design_Event_2D_Multi_Pop_splitt.R")
source("fit_gpd_cdf.R")
source("inv_gpd_cdf.R")
source("GPD_Threshold_Tails.R")


# -------------------------------------------------------------------------
# Analysis settings
# -------------------------------------------------------------------------

n_boot <- 500
n_ensembles <- 1000

# Representative sites used for the bootstrap sensitivity analysis
sites <- c(16, 39, 76)

# AEP levels
aep_levels <- c(0.5, 0.2, 0.1, 0.02, 0.01)

# Number of attempts allowed when resampling produces incompatible
# GPD shape combinations
max_attempts <- 15

# Leave one processor available for other system processes
n_cores <- max(1, parallel::detectCores() - 1)


# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------

matlab_datenum_to_date <- function(x) {
  as.Date(x - 719529, origin = "1970-01-01")
}


# -------------------------------------------------------------------------
# Load data
# -------------------------------------------------------------------------

Data <- readMat("Time_series_data.mat")

AR_data <- readMat("AR_data.mat")
AR_events <- lapply(AR_data$AR.data, `[[`, 1)

AR_events <- lapply(AR_events, function(x) {
  colnames(x) <- c("Time_RF", "RF", "Time_Hmo", "Hmo")
  x
})

non_AR_data <- readMat("non_AR_data.mat")
non_AR_events <- lapply(non_AR_data$non.AR.data, `[[`, 1)

non_AR_events <- lapply(non_AR_events, function(x) {
  colnames(x) <- c("Time_RF", "RF", "Time_Hmo", "Hmo")
  x
})


# -------------------------------------------------------------------------
# Parallel cluster
# -------------------------------------------------------------------------

cl <- makeCluster(n_cores)
registerDoParallel(cl)

cat("Parallel cluster started with", n_cores, "cores.\n")


# -------------------------------------------------------------------------
# Bootstrap sensitivity analysis
# -------------------------------------------------------------------------

for (i in sites) {

  set.seed(1234)

  original_AR_events <- AR_events[[i]]
  original_non_AR_events <- non_AR_events[[i]]

  # Full nTWL and rainfall time series
  site_data <- Data$Time.series.data[[i]][[1]]

  WL_above_normal <- site_data[, 5] + site_data[, 6]

  full_timeseries <- data.frame(
    Hmo = WL_above_normal,
    RF  = site_data[, 2]
  )

  # 24-hour accumulated rainfall
  k <- 24
  rf <- full_timeseries$RF

  rf24 <- rep(NA_real_, length(rf))

  rf24[k:length(rf)] <- vapply(
    k:length(rf),
    function(j) sum(rf[(j - k + 1):j], na.rm = TRUE),
    numeric(1)
  )

  full_timeseries$RF <- rf24
  full_timeseries$RF[is.na(full_timeseries$RF)] <- 0

  # Length of available event record
  yrs_AR <- diff(
    range(
      matlab_datenum_to_date(original_AR_events[, "Time_Hmo"]),
      na.rm = TRUE
    )
  )

  yrs_nonAR <- diff(
    range(
      matlab_datenum_to_date(original_non_AR_events[, "Time_Hmo"]),
      na.rm = TRUE
    )
  )

  years_available <- ceiling(
    as.numeric(max(c(yrs_AR, yrs_nonAR))) / 365.25
  )

  # Site-specific output folder
  output_dir <- paste0("Bootstrap_site_", i)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  cat(
    "\n----------------------------------------\n",
    "Site:", i, "\n",
    "Bootstrap samples:", n_boot, "\n",
    "Parallel cores:", n_cores, "\n",
    "----------------------------------------\n"
  )


  # -----------------------------------------------------------------------
  # Bootstrap loop
  # -----------------------------------------------------------------------

  Bootstrap_results <- foreach(
    b = seq_len(n_boot),

    .packages = c(
      "VineCopula",
      "VGAM",
      "tweedie",
      "statmod",
      "truncnorm",
      "MASS",
      "ks",
      "ismev",
      "dplyr",
      "texmex",
      "MultiHazard"
    ),

    .export = c(
      "original_AR_events",
      "original_non_AR_events",
      "full_timeseries",
      "years_available",
      "i",
      "n_ensembles",
      "aep_levels",
      "output_dir",
      "max_attempts",
      "Best_Copula",
      "Design_Event_2D_Multi_Pop_splitt",
      "fit_gpd_cdf",
      "inv_gpd_cdf",
      "GPD_Threshold_Tails"
    )

  ) %dopar% {

    attempt <- 0

    repeat {

      attempt <- attempt + 1

      # Reproducible but independent seed for each bootstrap and retry
      set.seed(1234 + b * 10000L + attempt)

      # Resample AR and non-AR events independently
      temp_AR_events_data <- original_AR_events[
        sample(
          seq_len(nrow(original_AR_events)),
          size = nrow(original_AR_events),
          replace = TRUE
        ),
        ,
        drop = FALSE
      ]

      temp_non_AR_events_data <- original_non_AR_events[
        sample(
          seq_len(nrow(original_non_AR_events)),
          size = nrow(original_non_AR_events),
          replace = TRUE
        ),
        ,
        drop = FALSE
      ]


      # -------------------------------------------------------------------
      # GPD threshold selection
      # -------------------------------------------------------------------

      thres_AR_Hmo <- GPD_Threshold_Tails(
        data_raw = temp_AR_events_data,
        record_length = years_available,
        baseline_event_RP = 0.25,
        threshold_probabilities = seq(0.5, 0.99, by = 0.005),
        parallel = FALSE,
        cores = 1,
        var_select = "Hmo",
        plot = FALSE
      )

      thres_AR_RF <- GPD_Threshold_Tails(
        data_raw = temp_AR_events_data,
        record_length = years_available,
        baseline_event_RP = 0.25,
        threshold_probabilities = seq(0.5, 0.99, by = 0.005),
        parallel = FALSE,
        cores = 1,
        var_select = "RF",
        plot = FALSE
      )

      thres_non_AR_Hmo <- GPD_Threshold_Tails(
        data_raw = temp_non_AR_events_data,
        record_length = years_available,
        baseline_event_RP = 0.25,
        threshold_probabilities = seq(0.5, 0.99, by = 0.005),
        parallel = FALSE,
        cores = 1,
        var_select = "Hmo",
        plot = FALSE
      )

      thres_non_AR_RF <- GPD_Threshold_Tails(
        data_raw = temp_non_AR_events_data,
        record_length = years_available,
        baseline_event_RP = 0.25,
        threshold_probabilities = seq(0.5, 0.99, by = 0.005),
        parallel = FALSE,
        cores = 1,
        var_select = "RF",
        plot = FALSE
      )


      # Use the second-lowest observation if no candidate threshold is found
      if (length(thres_AR_Hmo$Candidate_Thres$thresh) == 0) {
        thres_AR_Hmo$Candidate_Thres$thresh <-
          sort(temp_AR_events_data[, 4])[2]
      }

      if (length(thres_AR_RF$Candidate_Thres$thresh) == 0) {
        thres_AR_RF$Candidate_Thres$thresh <-
          sort(temp_AR_events_data[, 2])[2]
      }

      if (length(thres_non_AR_Hmo$Candidate_Thres$thresh) == 0) {
        thres_non_AR_Hmo$Candidate_Thres$thresh <-
          sort(temp_non_AR_events_data[, 4])[2]
      }

      if (length(thres_non_AR_RF$Candidate_Thres$thresh) == 0) {
        thres_non_AR_RF$Candidate_Thres$thresh <-
          sort(temp_non_AR_events_data[, 2])[2]
      }


      # -------------------------------------------------------------------
      # GPD fitting
      # -------------------------------------------------------------------

      AR_GPD_cdf_Hmo <- fit_gpd_cdf(
        data = temp_AR_events_data[, 3:4],
        threshold = thres_AR_Hmo$Candidate_Thres$thresh,
        n_years = years_available,
        plot = FALSE
      )

      AR_GPD_cdf_RF <- fit_gpd_cdf(
        data = temp_AR_events_data[, 1:2],
        threshold = thres_AR_RF$Candidate_Thres$thresh,
        n_years = years_available,
        plot = FALSE
      )

      non_AR_GPD_cdf_Hmo <- fit_gpd_cdf(
        data = temp_non_AR_events_data[, 3:4],
        threshold = thres_non_AR_Hmo$Candidate_Thres$thresh,
        n_years = years_available,
        plot = FALSE
      )

      non_AR_GPD_cdf_RF <- fit_gpd_cdf(
        data = temp_non_AR_events_data[, 1:2],
        threshold = thres_non_AR_RF$Candidate_Thres$thresh,
        n_years = years_available,
        plot = FALSE
      )


      # Repeat the resampling when AR and non-AR tails have the specified
      # incompatible shape combination
      bad_Hmo <- (
        non_AR_GPD_cdf_Hmo$xi > 0 &&
        AR_GPD_cdf_Hmo$xi < 0
      )

      bad_RF <- (
        non_AR_GPD_cdf_RF$xi > 0 &&
        AR_GPD_cdf_RF$xi < 0
      )

      condition_resolved <- !(bad_Hmo || bad_RF)

      if (condition_resolved || attempt >= max_attempts) {
        break
      }
    }


    # ---------------------------------------------------------------------
    # Copula fitting
    # ---------------------------------------------------------------------

    AR_data_for_cop <- data.frame(
      Hmo = temp_AR_events_data[, 4],
      RF  = temp_AR_events_data[, 2]
    )

    non_AR_data_for_cop <- data.frame(
      Hmo = temp_non_AR_events_data[, 4],
      RF  = temp_non_AR_events_data[, 2]
    )

    Copulas <- Best_Copula(
      Con_var_1 = AR_data_for_cop,
      Con_var_2 = non_AR_data_for_cop
    )

    AR_copula_Family <- Copulas$copula_Var1_Family
    non_AR_copula_Family <- Copulas$copula_Var2_Family


    # ---------------------------------------------------------------------
    # Joint AEP analysis
    # ---------------------------------------------------------------------

    # Redirect plotting produced by the design-event function
    tmp_pdf <- tempfile(fileext = ".pdf")
    pdf(tmp_pdf)

    out <- Design_Event_2D_Multi_Pop_splitt(
      Data               = full_timeseries,
      Data_Con1          = AR_data_for_cop,
      Data_Con4          = non_AR_data_for_cop,
      Copula_Family_pop1 = AR_copula_Family,
      Copula_Family_pop2 = non_AR_copula_Family,
      Con1               = "Rainfall",
      Con2               = "Hmo",
      GPD_P1_var1        = AR_GPD_cdf_Hmo,
      GPD_P1_var2        = AR_GPD_cdf_RF,
      GPD_P2_var1        = non_AR_GPD_cdf_Hmo,
      GPD_P2_var2        = non_AR_GPD_cdf_RF,
      Rate_Con1          = NA,
      Rate_Con4          = NA,
      RP                 = aep_levels,
      mu                 = 365.25 * 24,
      Isoline_Probs      = "Sample",
      N                  = 10^4,
      Sim_Max            = 10,
      y_lim_max          = 300,
      x_lim_max          = 20,
      N_Ensemble         = n_ensembles,
      Grid_x_min         = 0,
      Grid_x_max         = 20,
      Grid_y_min         = 0,
      Grid_y_max         = 300,
      Grid_x_interval    = 0.05,
      Grid_y_interval    = 0.5
    )

    dev.off()
    unlink(tmp_pdf)


    # ---------------------------------------------------------------------
    # Store bootstrap results
    # ---------------------------------------------------------------------

    result_b <- list(
      bootstrap = b,
      site = i,

      Grid = as.matrix(out$Grid),
      AEP_P1 = as.matrix(out$AEP_P1),
      AEP_P2 = as.matrix(out$AEP_P2),

      MostLikelyEvent = out$MostLikelyEvent,
      Ensemble = out$Ensemble,

      AR_Copula_Code = AR_copula_Family,
      Non_AR_Copula_Code = non_AR_copula_Family,

      AR_Copula_Name = BiCopName(
        AR_copula_Family,
        short = FALSE
      ),

      Non_AR_Copula_Name = BiCopName(
        non_AR_copula_Family,
        short = FALSE
      ),

      AR_Hmo_threshold =
        thres_AR_Hmo$Candidate_Thres$thresh,

      AR_RF_threshold =
        thres_AR_RF$Candidate_Thres$thresh,

      non_AR_Hmo_threshold =
        thres_non_AR_Hmo$Candidate_Thres$thresh,

      non_AR_RF_threshold =
        thres_non_AR_RF$Candidate_Thres$thresh,

      AR_Hmo_scale = AR_GPD_cdf_Hmo$sigma,
      AR_Hmo_shape = AR_GPD_cdf_Hmo$xi,

      AR_RF_scale = AR_GPD_cdf_RF$sigma,
      AR_RF_shape = AR_GPD_cdf_RF$xi,

      non_AR_Hmo_scale = non_AR_GPD_cdf_Hmo$sigma,
      non_AR_Hmo_shape = non_AR_GPD_cdf_Hmo$xi,

      non_AR_RF_scale = non_AR_GPD_cdf_RF$sigma,
      non_AR_RF_shape = non_AR_GPD_cdf_RF$xi,

      attempts = attempt,
      condition_resolved = condition_resolved
    )

    saveRDS(
      result_b,
      file = file.path(
        output_dir,
        sprintf("site_%d_bootstrap_%03d.rds", i, b)
      )
    )

    result_b
  }


  # -----------------------------------------------------------------------
  # Save all bootstrap results for the site
  # -----------------------------------------------------------------------

  Grid_ref <- Bootstrap_results[[1]]$Grid

  for (b in seq_along(Bootstrap_results)) {

    if (!isTRUE(all.equal(
      Grid_ref,
      Bootstrap_results[[b]]$Grid
    ))) {
      warning(
        paste(
          "Grid changed in bootstrap",
          b,
          "for site",
          i
        )
      )
    }
  }

  saveRDS(
    list(
      Site_ID = i,
      Grid = Grid_ref,
      Bootstrap_results = Bootstrap_results
    ),
    file = paste0(
      "Bootstrap_All_Results_site_",
      i,
      ".rds"
    )
  )

  cat(
    "Finished",
    n_boot,
    "bootstrap samples for site",
    i,
    "\n"
  )
}


# -------------------------------------------------------------------------
# Finish
# -------------------------------------------------------------------------

stopCluster(cl)

end_time <- Sys.time()
runtime <- end_time - start_time

cat("\nTotal runtime:", as.character(runtime), "\n")