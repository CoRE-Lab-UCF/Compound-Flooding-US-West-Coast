# =============================================================================
# Inverse CDF: map probability -> real-scale value
# =============================================================================

# -----------------------------------------------------------------------------
#' inv_gpd_cdf
#'
#' Given a fitted GPD model (from fit_gpd_cdf) and one or more CDF
#' probabilities, returns the corresponding real-scale values and plots
#' original values vs back-transformed values for visual validation.
#'
#'   Below threshold  ->  empirical quantile (linear interpolation via approx)
#'   Above threshold  ->  analytical GPD inverse:
#'       x = u + (sigma/xi) * [ ((1-p)/(1-F_u))^(-xi) - 1 ]
#'   Special case xi ~ 0:
#'       x = u - sigma * log( (1-p)/(1-F_u) )
#'
#' @param p           Numeric vector of probabilities in (0, 1).
#' @param fit         List returned by fit_gpd_cdf().  Must contain:
#'                      $sigma, $xi, $threshold, $F_threshold, $values
#' @param col_orig    Colour for original values in the comparison plot
#'                    (default "steelblue").
#' @param col_back    Colour for back-transformed values in the comparison plot
#'                    (default "firebrick").
#' @param main        Plot title (default auto-generated).
#'
#' @return Numeric vector of real-scale values, same length as p.
# -----------------------------------------------------------------------------
inv_gpd_cdf <- function(p,
                        fit,
                        col_orig = "steelblue",
                        col_back = "firebrick",
                        main     = NULL,
                        plot     = TRUE) {
  
  # -- Input checks -----------------------------------------------------------
  if (any(p <= 0 | p >= 1))
    stop("All probabilities in `p` must be strictly between 0 and 1.")
  if (is.null(fit$values))
    stop("`fit$values` is missing. Re-run fit_gpd_cdf() to get an updated fit object.")
  
  sigma  <- fit$sigma
  xi     <- fit$xi
  u      <- fit$threshold
  F_u    <- fit$F_threshold
  values <- fit$values          # pulled directly from the fit object
  
  # -- Split p into below / above threshold -----------------------------------
  below  <- p <= F_u
  above  <- p >  F_u
  result <- numeric(length(p))
  
  # -- Below threshold: linear interpolation of empirical quantile ------------
  if (any(below)) {
    emp_sorted <- sort(values[values <= u])
    n_total    <- length(values)
    emp_probs  <- seq_along(emp_sorted) / n_total
    
    result[below] <- approx(x    = emp_probs,
                            y    = emp_sorted,
                            xout = p[below],
                            rule = 2)$y
  }
  
  # -- Above threshold: analytical GPD inverse --------------------------------
  if (any(above)) {
    p_above <- p[above]
    
    if (abs(xi) < 1e-10) {
      result[above] <- u - sigma * log((1 - p_above) / (1 - F_u))
    } else {
      result[above] <- u + (sigma / xi) *
        (((1 - p_above) / (1 - F_u))^(-xi) - 1)
    }
  }
  
  # -- Plot: original values vs back-transformed values (optional) -----------
  if (plot) {
    if (is.null(main)) main <- "Original vs Back-Transformed Values"
    
    n_pts   <- length(p)
    x_index <- seq_len(n_pts)
    
    # Sort both series by the back-transformed value for a clean comparison
    ord          <- order(result)
    orig_sorted  <- sort(values)[seq_len(n_pts)]   # same-length sorted original
    back_sorted  <- result[ord]
    
    y_range <- range(c(orig_sorted, back_sorted), na.rm = TRUE)
    
    plot(x_index, orig_sorted,
         type = "l",
         col  = col_orig,
         lwd  = 2,
         xlab = "Rank",
         ylab = "Value",
         main = main,
         ylim = y_range,
         las  = 1)
    
    abline(h = pretty(y_range), col = "grey92", lty = 1)
    abline(v = pretty(x_index), col = "grey92", lty = 1)
    box()
    
    lines(x_index, orig_sorted, col = col_orig, lwd = 2)   # redraw over grid
    lines(x_index, back_sorted, col = col_back, lwd = 2, lty = 2)
    
    abline(h   = u,
           col = "grey40", lty = 3, lwd = 1.2)
    text(x      = 1,
         y      = u,
         labels = sprintf("u = %.4g", u),
         pos    = 4, col = "grey35", cex = 0.85)
    
    legend("topleft",
           legend = c("Original values (sorted)",
                      "Back-transformed from CDF"),
           col    = c(col_orig, col_back),
           lty    = c(1, 2),
           lwd    = c(2, 2),
           bty    = "n", cex = 0.88)
    
  } # end if (plot)
  
  invisible(result)
}