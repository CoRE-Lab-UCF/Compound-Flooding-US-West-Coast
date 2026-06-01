GPD_Threshold_Tails <- function(
    data_raw,
    record_length,
    baseline_event_RP = 0.25,
    threshold_probabilities = seq(0.5, 0.99, by = 0.01),
    parallel = TRUE,
    cores = 4,
    plot = TRUE,
    var_select=NULL
) {
  
  #' Threshold selection method for univariate extremes
  #'
  #' 'thresh_select_function' selects a constant threshold which best captures the tail of the data via a Generalised Pareto distribution.
  #'
  #' @author Tom Collings et al. 
  #'
  #' @param data A numeric vector. The data for which to select a threshold.
  #' @param thresh A numeric vector of proposed thresholds to test.
  #' @param threshold_probabilities A numeric vector of the corresponding threshold probabilities. Can be null
  #' @param baseline_probability A probability value above which we define the tail of the data.
  #' @param ppy A positive value denoting the number of observations per year
  #' @param k  A positive integer denoting the number of bootstraps.
  #' @param m A positive integer denoting the number of equally-spaced probabilities at which to evaluate quantiles.
  #' @param parallel A TRUE/FALSE value specifying whether or not to run the code in parallel. 
  #' @param cores A positive integer denoting the number of cores for running in parallel. 
  #'
  #' @returns A list containing the chosen threshold, the parameters of the fitted GPD, the number of observations above the chosen thresholds and the metric values 'd' corresponding to each proposed threshold.
  #'
  #' @examples
  #' set.seed(12345)
  #' data_test1 <- rgpd(1000, shape = 0.1, scale=0.5, mu=1)
  #' thresholds1 <- quantile(data_test1,seq(0,0.95,by=0.05))
  #' baseline_probability <- 0.5
  #' ppy <- 365
  #' (example <- thresh_select_function(data = data_test, thresh = thresholds,baseline_probability = baseline_probability,ppy = ppy))
  
  thresh_select_function <- function(data, thresh, threshold_probabilities=NULL, baseline_probability,ppy,k = 100, m = 500,parallel = FALSE,cores=NULL){
    
    # Check inputs are valid
    if (!is.numeric(data)) stop("data must be a vector")
    if (!is.numeric(thresh)) stop("thresh to be tested needs to be a vector")
    if (!is.numeric(threshold_probabilities) & !is.null(threshold_probabilities)) stop("threshold_probabilities need to be a vector or NULL")
    if (!is.numeric(baseline_probability) | length(baseline_probability)!=1 | baseline_probability < 0 | baseline_probability > 1) stop("baseline_probability must be a single value between 0 and 1")
    if (k <= 0 | k %% 1 != 0) stop("Number of bootstrapped samples (k) must be a positive integer")
    if (m <= 0 | m %% 1 != 0) stop("Number of equally spaced probabilities (m) must be a positive integer")
    if (!is.numeric(ppy) | length(ppy) != 1) stop("ppy must be a single numeric value")
    if (!is.logical(parallel)) stop("parallel can only take the values 'TRUE' or 'FALSE'")
    
    #length of data 
    n = length(data)
    
    if(is.null(threshold_probabilities)){
      ##### convert thresholds to empirical exceedance probabilities ######
      # Create the empirical CDF from the data
      ecdf_func <- ecdf(data)
      
      # Evaluate the empirical CDF at the selected quantiles
      threshold_probabilities <- ecdf_func(thresh)
    }
    
    #select probabilities of interest
    if (baseline_probability>1 - (10/n)) stop("Baseline event probability is too large - please decrease.")
    quantile_probs <- seq(baseline_probability,1 - (10/n),length=m) 
    
    if(parallel == FALSE){ #if parallel is set to false, we run in a for loop, which takes longer 
      
      meandistances <- xis <- sigmas <- num_excess <- numeric(length(thresh))
      
      for (i in 1:length(thresh)) { #for each threshold
        u <- thresh[i] #select threshold
        excess <- data[data > u] - u #exceedences over threshold
        num_excess[i] <- length(excess) #number of exceedences over threshold
        threshold_probability <- threshold_probabilities[i] #probability of the threshold compared to the empirical distribution
        
        if (threshold_probability <= 1 - (1/ppy) ) { #We only fit the model if the threshold is less than the 1 in 1 year empirical return level
          mle0 <- mean(excess)
          init.fit <- optim(GPD_LL, z = excess, par = c(mle0,0.1), control = list(fnscale = -1))
          xis[i] <- init.fit$par[[2]]
          sigmas[i] <- init.fit$par[[1]]
          distances <- numeric(k)
          
          #computing theoretical quantiles. This remains constant over bootstrapping 
          exceedance_probs <- 1 - (1-quantile_probs)/(1-threshold_probability)
          
          # mask to remove any exceedance probabilities that are less than 0 or equal to 1 (i.e., not valid)
          # note that this accounts for the cases when the baseline event is lower than the threshold - we just estimate quantiles above the threshold in t
          mask <- exceedance_probs > 0 & exceedance_probs < 1
          
          for (j in 1:k) { #for each bootstrap
            X <- sample(excess, num_excess[i], replace = TRUE) #sample with replacement from the exceedences
            mle <- mean(X)
            ifelse(xis[i] < 0, pars_init <-  c(mle, 0.1) ,pars_init <- c(sigmas[i], xis[i]) )
            gpd.fit <- optim(GPD_LL, z = X, par = pars_init, control = list(fnscale = -1)) #fit the gpd to the sample
            
            ############ using matched quantiles (from callum's code)
            
            theoretical_quantiles <- qgpd(exceedance_probs[mask], scale = gpd.fit$par[[1]], shape = gpd.fit$par[[2]])
            empirical_quants <- quantile(X,probs=exceedance_probs[mask]) #to get the empirical quantiles for this iteration
            
            #find the errors between empirical and theoretical quantiles
            errors <- abs(empirical_quants - theoretical_quantiles) 
            
            #calculate mean errors - EQD
            distances[j] <- (1/length(errors)) * sum(errors) # errors over the threshold is inheriant in the quantile selection
          }
          
          meandistances[i] <- mean(distances) #find mean across all bootstraps
          
        } else{
          meandistances[i] <- NA
        }
      }
      
      chosen_index <- which.min(meandistances)
      chosen_threshold <- thresh[chosen_index]
      chosen_threshold_prob <- threshold_probabilities[chosen_index]
      xi <- xis[chosen_index]
      sigma <- sigmas[chosen_index]
      len <- num_excess[chosen_index]
      
    } else {
      
      if(is.null(cores)){
        cores = floor(0.5*detectCores())
      } else {
        if (!is.numeric(cores) | length(cores) != 1 | cores %% 1 != 0) stop("cores must be a single integer value, or NULL")
        if (cores > detectCores()-1) stop("cores is too high for your CPU setup. Please reduce")
      }
      
      cl <- makeCluster(cores)
      
      clusterExport(cl, list("data","GPD_LL","quantile_probs","k","m","qgpd","ppy"), envir = environment())
      
      TS_info <- parApply(cl, rbind(thresh,threshold_probabilities),2, function(xcol) {
        
        u <- xcol[1]
        excess <- data[data > u] - u #exceedences over threshold
        num_excess <- length(excess) #number of exceedences over threshold
        threshold_probability <- xcol[2] #probability of the threshold compared to the empirical distribution
        
        if (threshold_probability <= 1 - (1/ppy)) { #We only fit the model if the threshold is less than the 1 in 1 year empirical return level
          mle0 <- mean(excess)
          init.fit <- optim(GPD_LL, z = excess, par = c(mle0,0.1), control = list(fnscale = -1))
          xi <- init.fit$par[[2]]
          sigma <- init.fit$par[[1]]
          distances <- numeric(k)
          
          #computing theoretical quantiles. This remains constant over bootstrapping 
          exceedance_probs <- 1 - (1-quantile_probs)/(1-threshold_probability)
          
          # mask to remove any exceedance probabilities that are less than 0 or equal to 1 (i.e., not valid)
          # note that this accounts for the cases when the baseline event is lower than the threshold - we just estimate quantiles above the threshold in t
          mask <- exceedance_probs > 0 & exceedance_probs < 1
          
          for (j in 1:k) { #for each bootstrap
            X <- sample(excess, num_excess, replace = TRUE) #sample with replacement from the exceedences
            mle <- mean(X)
            ifelse(xi < 0, pars_init <-  c(mle, 0.1) ,pars_init <- c(sigma, xi) )
            gpd.fit <- optim(GPD_LL, z = X, par = pars_init, control = list(fnscale = -1)) #fit the gpd to the sample
            
            ############ using matched quantiles (from callum's code)
            
            theoretical_quantiles <- qgpd(exceedance_probs[mask], scale = gpd.fit$par[[1]], shape = gpd.fit$par[[2]])
            empirical_quants <- quantile(X,probs=exceedance_probs[mask]) #to get the empirical quantiles for this iteration
            
            #find the errors between empirical and theoretical quantiles
            errors <- abs(empirical_quants - theoretical_quantiles) 
            
            #calculate mean errors - EQD
            distances[j] <- (1/length(errors)) * sum(errors) # errors over the threshold is inheriant in the quantile selection
            
          }
          
          meandistance <- mean(distances) #find mean across all bootstraps
        } else{
          meandistance <- NA
          xi <- NA
          sigma <- NA
          num_excess <- NA
        }
        
        return(c(num_excess,sigma,xi,meandistance))
      })
      
      # Stop parallel processing
      stopCluster(cl)
      
      TS_info = t(TS_info)
      
      colnames(TS_info) = c("n_exc","sigmas","xis","dists")
      
      TS_info = as.data.frame(TS_info)
      
      meandistances <- TS_info$dists
      chosen_index <- which.min(meandistances)
      chosen_threshold <- thresh[chosen_index]
      chosen_threshold_prob <- threshold_probabilities[chosen_index]
      xi <- TS_info$xis[chosen_index]
      sigma <- TS_info$sigmas[chosen_index]
      len <- TS_info$n_exc[chosen_index]
      
    }
    
    result <- list(
      thresh = chosen_threshold, 
      threshold_prob = chosen_threshold_prob,
      par = c(sigma,xi), 
      num_excess = len, 
      dists = meandistances
    )
    
    return(result)
    
  }
  
  
  #=====================================================================
  # Functions for Generalised Pareto Distribution.
  # Added option to use nu parameterisation
  # checks that param values are valid
  #=====================================================================
  # pgpd
  # qgpd
  # dgpd
  # rgpd
  # GPD_LL
  #=====================================================================
  #' Generalised Pareto Distribution
  #'
  #' Cumulative density function of the GPD specified in terms of (sigma,xi) or (nu,xi).
  #' Improvement over evir function as it returns an error if the (implied) shape
  #' parameter is non-positive. Also properly handles cases where and xi=0 or p<mu.
  #'
  #' @author Zak Varty
  #'
  #' @param q vector of quantiles.
  #' @param shape shape parameter (xi)
  #' @param scale scale parameter (sigma)
  #' @param nu  alternative scale parameter: nu = sigma/(1+xi)
  #' @param mu  location parameter
  #' @param skip_checks logical. Speed up evaluation by skipping checks on inputs? (Beware!)
  #' @return Probability of the GPD X<=q
  #' @importFrom stats pexp
  #' @examples
  #' pgpd(q = c(-1,1.5,3), shape = 1, scale = 1)
  #' pgpd(q = 1.5, shape = c(0,-1), scale = c(0.1,1))
  #' @export
  pgpd <- function(q, shape, scale = NULL, nu = NULL, mu = 0, skip_checks = FALSE){
    
    if (!skip_checks) {
      # one and only one of {nu, scale} may be specified
      if (is.null(scale) & is.null(nu)) {
        stop('Define one of the parameters nu or scale.')
      }
      if (!is.null(scale) & !is.null(nu)) {
        stop('Define only one of the parameters nu and scale.')
      }
      # Calculate scale from nu if required
      if (!is.null(nu) & is.null(scale)) {
        scale <- nu / (1 + shape)
        if (any(scale <= 0)) {
          stop('Implied scale parameter(s) must be positive.')
        }
        
      }
      # Check that scale value(s) are positive
      if (any(scale <= 0)) {
        stop('Scale parameter(s) must be positive.')
      }
      
      # Ensure q, scale, shape and mu are of same length.
      if (length(scale) == 1 & length(q) > 1) {
        scale <- rep(scale, length(q))
      }
      if (length(shape) == 1 & length(q) > 1) {
        shape <- rep(shape, length(q))
      }
      if (length(mu) == 1 & length(q) > 1) {
        mu <- rep(mu, length(q))
      }
    } else {
      if (!is.null(nu) & is.null(scale)) {
        scale <- nu / (1 + shape)
      }
    }
    #calculate probabilities
    p <- (1 - (1 + (shape * (q - mu))/scale)^(-1/shape))
    #correct probabilities below mu or above upper end point
    p[q < mu] <- 0
    p[(shape < 0) & (q >= (mu - scale/shape))] <- 1
    
    #correct probabilities where xi = 0
    if (any(abs(shape) < 1e-10)) {
      #ex <- which(shape ==0)
      ex <- which(abs(shape) < 1e-10)
      p[ex] <- pexp(q = q[ex] - mu[ex], rate = 1 / scale[ex])
    }
    
    return(p)
  }
  
  #' Generalised Pareto Distribution
  #'
  #' Cumulative density function of the GPD specified in terms of (sigma,xi) or (nu,xi).
  #' Improvement over evir function as it returns an error if the (implied) shape
  #' parameter is non-positive. Also properly handles cases where and xi=0 or p is not a valid
  #' probability.
  #'
  #' @author Zak Varty
  #'
  #' @param p vector of quantiles.
  #' @param shape shape parameter (xi)
  #' @param scale scale parameter (sigma)
  #' @param nu  alternative scale parameter: nu = sigma/(1+xi)
  #' @param mu  location parameter
  #' @return Probability of the GPD X<=x
  #' @examples
  #' qgpd(p = 0.5, shape = 0.5, scale = 0.5)
  #' \dontrun{ qgpd(p = -0.1, shape = 0, scale = 1, mu = 0.1) }
  #' @export
  qgpd <- function(p, shape, scale = NULL, nu = NULL, mu = 0){
    # one and only one of {nu, scale} may be specified
    if (is.null(scale) & is.null(nu)) {
      stop('Define one of the parameters nu or scale.')
    }
    if (!is.null(scale) & !is.null(nu)) {
      stop('Define only one of the parameters nu and scale.')
    }
    
    # Probabilities must all be positive
    if (!all((p >= 0) & (p <= 1))) {
      stop('Probabilities p must be in the range [0,1].')
    }
    
    # Calculate scale from nu if required
    if (!is.null(nu) & is.null(scale)) {
      scale <- nu / (1 + shape)
      if (any(scale <= 0)) {
        stop('Implied scale parameter(s) must be positive.')
      }
      
    }
    
    # Check that scale value(s) are positive
    if (any(scale <= 0)) {
      stop('Scale parameter(s) must be positive.')
    }
    # Ensure p, scale, shape and mu are of same length.
    if (length(scale) == 1 & length(p) > 1) {
      scale <- rep(scale, length(p))
    }
    if (length(shape) == 1 & length(p) > 1) {
      shape <- rep(shape, length(p))
    }
    if (length(mu) == 1 & length(p) > 1) {
      mu <- rep(mu, length(p))
    }
    
    #calculate quantiles
    q <- mu + (scale/shape) * ((1 - p)^(-shape) - 1)
    
    #correct quantiles where xi = 0
    #ex <- which(shape ==0)
    if (any(abs(shape) < 1e-10)) {
      ex <- which(abs(shape) < 1e-10)
      q[ex] <- mu[ex] + stats::qexp(p = p[ex],rate = 1/scale[ex])
    }
    return(q)
  }
  
  #' Generalised Pareto Distribution
  #'
  #' Density function of the GPD specified in terms of (sigma,xi) or (nu,xi).
  #' Improvement over evir function as it returns an error if the (implied) shape
  #' parameter is non-positive. Also properly handles cases where and xi=0 or x is
  #' outside of the domain of the given distribution.
  #'
  #' @author Zak Varty
  #'
  #' @param x vector of values as which to evaluate density.
  #' @param shape shape parameter (xi)
  #' @param scale scale parameter (sigma)
  #' @param nu  alternative scale parameter
  #' @param mu  location parameter
  #' @param log  locical. Return log
  #' @return density of the GPD at x
  #' @examples
  #' dgpd(x = c(-1,0.5,1,1.9,5),shape = -0.5, scale = 1)
  #' @export
  #'
  dgpd <- function(x, shape, scale = NULL, nu = NULL, mu = 0, log = FALSE){
    # one and only one of {nu, scale} may be specified
    if (is.null(scale) & is.null(nu)) {
      stop('Define one of the parameters nu or scale.')
    }
    if (!is.null(scale) & !is.null(nu)) {
      stop('Define only one of the parameters nu and scale.')
    }
    
    # Calculate scale from nu if required
    if (!is.null(nu) & is.null(scale)) {
      scale <- nu / (1 + shape)
      if (any(scale <= 0)) {
        stop('Implied scale parameter(s) must be positive.')
      }
    }
    
    # Check that scale value(s) are positive
    if (any(scale <= 0)) {
      stop('Scale parameter(s) must be positive.')
    }
    # Ensure x, scale, shape and mu are of same length.
    if (length(scale) == 1 & length(x) > 1) {
      scale <- rep(scale, length(x))
    }
    if (length(shape) == 1 & length(x) > 1) {
      shape <- rep(shape, length(x))
    }
    if (length(mu) == 1 & length(x) > 1) {
      mu <- rep(mu, length(x))
    }
    
    if (log == FALSE) {
      out <- (scale^(-1)) * pmax((1 + shape * (x - mu)/scale),0)^((-1/shape) - 1)
      # amend values below threshold
      out[which(x < mu)] <- 0
      # amend values above upper endpoint (if it exists)
      out[which((shape < 0) & (x >= (mu - scale/shape)))] <- 0
      # amend values where xi = 0 (if they exist)
      if (any(abs(shape < 1e-10))) {
        ex <- which(abs(shape) < 1e-10)
        out[ex] <- stats::dexp(x = x[ex] - mu[ex], rate = 1/scale[ex])
      }
    } else {
      out <-  -log(scale) + ((-1/shape) - 1)*log(pmax((1 + shape * (x - mu)/scale),0))
      # amend values below threshold
      out[which(x < mu)] <- -Inf
      # amend values above upper endpoint (if it exists)
      out[which((shape < 0) & (x >= (mu - scale/shape)))] <- -Inf
      # amend values where xi = 0 (if they exist)
      if (any(abs(shape) < 1e-10)) {
        ex <- which(abs(shape) < 1e-10)
        out[ex] <- stats::dexp(x = x[ex] - mu[ex], rate = 1 / scale[ex],log = TRUE)
      }
    }
    return(out)
  }
  
  #' Generalised Pareto Distribution
  #'
  #' Sample the GPD specified in terms of (sigma,xi) or (nu,xi).
  #' Improvement over evir function as it returns an error if the (implied) shape
  #' parameter is non-positive. Also properly handles cases where and xi=0.
  #'
  #' @author Zak Varty
  #'
  #' @param n sample size.
  #' @param shape shape parameter (xi).
  #' @param scale scale parameter (sigma).
  #' @param nu  alternative scale parameter.
  #' @param mu  location parameter.
  #' @return Random sample from generalised pareto distirbution.
  #'
  #' @examples
  #' rgpd(n = 100, shape = 0, scale = 1:100)
  #' @export
  rgpd <- function(n, shape, scale = NULL, nu = NULL, mu = 0){
    ## Input checks
    # one and only one of {nu, scale} may be specified
    if (is.null(scale) & is.null(nu)) {
      stop('Define one of the parameters nu or scale.')
    }
    if (!is.null(scale) & !is.null(nu)) {
      stop('Define only one of the parameters nu and scale.')
    }
    # Calculate scale from nu if required
    if (!is.null(nu) & is.null(scale)) {
      scale <- nu / (1 + shape)
      if (any(scale <= 0)) {
        stop('Implied scale parameter(s) must be positive.')
      }
    }
    # Check that scale value(s) are positive
    if (any(scale <= 0)) {
      stop('Scale parameter(s) must be positive.')
    }
    # Ensure q, scale, shape and mu are of same length.
    if ((length(scale) == 1) & (n > 1)) {
      scale <- rep(scale, n)
    }
    if ((length(shape) == 1) & (n > 1)) {
      shape <- rep(shape, n)
    }
    if ((length(mu) == 1) & (n > 1)) {
      mu <- rep(mu, n)
    }
    
    #simulate sample
    sample <- mu + (scale/shape) * ((1 - stats::runif(n))^(-shape) - 1)
    #correct sample values where xi = 0
    #ex <- which(shape ==0)
    if (any(abs(shape) < 1e-10)) {
      ex <- which(abs(shape) < 1e-10)
      sample[ex] <- mu[ex] +
        stats::rexp(n = length(ex),rate = 1/scale[ex])
    }
    return(sample)
  }
  
  #' Generalised Pareto log-likelihood
  #'
  #' @author Conor Murphy
  #'
  #' @param par A numeric vector of parameter values of length 2.
  #' @param z A numeric vector of excesses of some threshold.
  #'
  #' @returns A numeric value of the log-likeihood.
  #'
  #' @examples
  #' test1 <- rgpd(1000, shape = 0.1, scale=0.5, mu=1)
  #' excess <- test1[test1>1.5] - 1.5
  #' GPD_LL(par=c(1,0.4), z=excess)
  #' @export
  GPD_LL <- function(par, z){
    sigma <- par[1]
    xi <- par[2]
    if (sigma > 0) {
      if (abs(xi) < 1e-10) {
        return(-length(z) * log(sigma) - ((1 / sigma) * sum(z)))
      }
      else {
        if (all(1 + (xi * z) / sigma > 0)) {
          return(-(length(z) * log(sigma)) - ((1 + 1 / xi)*(sum(log(1 + (xi * z) / sigma)))))
        }
        else{
          return(-1e6)
        }
      }
    }
    else{
      return(-1e7)
    }
  }
  
  
  
  
  
  
  
  
  
  
  #################### Above are helpers #################
  
  # If data is a data.frame, extract sea_level column
  event_col      <- switch(var_select, "Hmo" = 4, "RF" = 2)
  data <- data_raw[, event_col]
  
  
  
  # Compute derived parameters
  points_per_year      <- length(data) / record_length
  while ((1 - 1 / (baseline_event_RP * points_per_year)) < 0) {
    baseline_event_RP <- baseline_event_RP + 0.1
  }
  baseline_probability <- 1 - 1 / (baseline_event_RP * points_per_year)
  thresholds           <- unname(quantile(data, probs = threshold_probabilities))
  
  
  # Run threshold selection algorithm
  thresh_select <- thresh_select_function(
    data                    = data,
    thresh                  = thresholds,
    threshold_probabilities = threshold_probabilities,
    baseline_probability    = baseline_probability,
    parallel                = parallel,
    cores                   = cores,
    ppy                     = points_per_year
  )
  

  # Optional plotting
  if (plot) {
    par(mfrow = c(1, 2), mgp = c(2.5, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
    
    valid_idx <- !is.na(thresh_select$dists)
    
    plot(thresholds, thresh_select$dists,
         xlim = c(min(thresholds), max(thresholds[valid_idx])),
         type = "l", lwd = 3, col = "grey",
         xlab = expression(u), ylab = expression(tilde(d)(u)),
         main = "Distance metric vs threshold",
         cex.lab = 1.3, cex.axis = 1.5, cex.main = 1.5)
    abline(v = thresh_select$thresh, lwd = 4, col = "grey1")
    
    plot(threshold_probabilities, thresh_select$dists,
         xlim = c(min(threshold_probabilities), max(threshold_probabilities[valid_idx])),
         type = "l", lwd = 3, col = "grey",
         xlab = expression(p), ylab = expression(tilde(d)(u)),
         main = "Distance metric vs non-exceedance probability",
         cex.lab = 1.3, cex.axis = 1.5, cex.main = 1.5)
    abline(v = thresh_select$threshold_prob, lwd = 4, col = "grey1")
  }
  
  # Return results
  invisible(list(
    Candidate_Thres           = thresh_select,
    thresholds              = thresholds,
    threshold_probabilities = threshold_probabilities,
    points_per_year         = points_per_year,
    baseline_probability    = baseline_probability
  ))
}