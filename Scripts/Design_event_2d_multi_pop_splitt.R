Design_Event_2D_Multi_Pop_splitt <- function(Data, Data_Con1, Data_Con4, Copula_Family_pop1, Copula_Family_pop2, Con1="Rainfall", Con2="OsWL", GPD_P1_var1, GPD_P1_var2, GPD_P2_var1, GPD_P2_var2, Rate_Con1=NA, Rate_Con4=NA, mu=365.25*24, Grid_x_min=NA, Grid_x_max=NA, Grid_y_min=NA, Grid_y_max=NA, Grid_x_interval=NA, Grid_y_interval=NA, RP, x_lab="Rainfall (mm)", y_lab="O-sWL (mNGVD 29)", x_lim_min=0, x_lim_max=10, y_lim_min=0, y_lim_max=500, Isoline_Probs="Sample", N=10^6, N_Ensemble=0, Sim_Max=10){
  
  # --- Preliminaries ---
  Isoline         <- vector(mode="list", length=length(RP)); names(Isoline)         <- RP
  Contour         <- vector(mode="list", length=length(RP)); names(Contour)         <- RP
  Ensemble        <- vector(mode="list", length=length(RP)); names(Ensemble)        <- RP
  MostLikelyEvent <- vector(mode="list", length=length(RP)); names(MostLikelyEvent) <- RP
  FullDependence  <- vector(mode="list", length=length(RP)); names(FullDependence)  <- RP
  
  if(inherits(Data[,1], c("Date","factor","POSIXct"))){ Data <- Data[,-1] }
  
  # --- Local GPD CDF transform ---
  transform_to_gpd_cdf <- function(new_values, fit_result) {
    threshold <- fit_result$threshold
    sigma     <- fit_result$sigma
    xi        <- fit_result$xi
    F_u       <- fit_result$F_threshold
    cdf_out   <- numeric(length(new_values))
    for(i in seq_along(new_values)){
      x <- new_values[i]
      if(x <= threshold){
        cdf_out[i] <- approx(x=fit_result$cdf_x, y=fit_result$cdf_values, xout=x, rule=2)$y
      } else {
        G <- pgpd(x - threshold, sigma=sigma, xi=xi)
        cdf_out[i] <- F_u + (1 - F_u) * G
      }
    }
    return(cdf_out)
  }
  
  # --- Grid ---
  Grid_x_min      <- ifelse(is.na(Grid_x_min),      min(Data[,1], na.rm=TRUE), Grid_x_min)
  Grid_x_max      <- 2 * ifelse(is.na(Grid_x_max),  max(Data[,1], na.rm=TRUE), Grid_x_max)
  Grid_y_min      <- ifelse(is.na(Grid_y_min),      min(Data[,2], na.rm=TRUE), Grid_y_min)
  Grid_y_max      <- 2 * ifelse(is.na(Grid_y_max),  max(Data[,2], na.rm=TRUE), Grid_y_max)
  Grid_x_interval <- ifelse(is.na(Grid_x_interval), 2,   Grid_x_interval)
  Grid_y_interval <- ifelse(is.na(Grid_y_interval), 0.1, Grid_y_interval)
  
  var1  <- seq(Grid_x_min, Grid_x_max, Grid_x_interval)
  var2  <- seq(Grid_y_min, Grid_y_max, Grid_y_interval)
  Pgrid <- expand.grid(var1, var2)
  
  # --- Column indices ---
  con1 <- which(names(Data) == Con1)
  con2 <- which(names(Data) == Con2)
  
  # --- Occurrence rates ---
  time.period <- nrow(Data[which(!is.na(Data[,1]) & !is.na(Data[,2])),]) / mu
  if(is.na(Rate_Con1)) Rate_Con1 <- nrow(Data_Con1) / time.period
  if(is.na(Rate_Con4)) Rate_Con4 <- nrow(Data_Con4) / time.period
  
  # --- Copula fit & simulation: Population 1 ---
  obj1      <- BiCopSelect(pobs(Data_Con1[,1]), pobs(Data_Con1[,2]), familyset=Copula_Family_pop1, selectioncrit="AIC", indeptest=FALSE, level=0.05, weights=NA, rotations=TRUE, se=FALSE, presel=TRUE, method="mle")
  cop.sim1  <- BiCopSim(round(N * Rate_Con1 / (Rate_Con1 + Rate_Con4), 0), obj1)  # <<< renamed from 'sample'
  cop.sample.pop.1 <- data.frame(Var1=inv_gpd_cdf(p=cop.sim1[,1], fit=GPD_P1_var1, plot=FALSE),
                                 Var2=inv_gpd_cdf(p=cop.sim1[,2], fit=GPD_P1_var2, plot=FALSE))
  
  # --- Copula fit & simulation: Population 2 ---
  obj4      <- BiCopSelect(pobs(Data_Con4[,1]), pobs(Data_Con4[,2]), familyset=Copula_Family_pop2, selectioncrit="AIC", indeptest=FALSE, level=0.05, weights=NA, rotations=TRUE, se=FALSE, presel=TRUE, method="mle")
  cop.sim4  <- BiCopSim(round(N * Rate_Con4 / (Rate_Con1 + Rate_Con4), 0), obj4)  # <<< renamed from 'sample'
  cop.sample.pop.2 <- data.frame(Var1=inv_gpd_cdf(p=cop.sim4[,1], fit=GPD_P2_var1, plot=FALSE),
                                 Var2=inv_gpd_cdf(p=cop.sim4[,2], fit=GPD_P2_var2, plot=FALSE))
  
  cop.sample <- rbind(cop.sample.pop.1, cop.sample.pop.2)
  
  # --- Transform grid to (0,1) and evaluate copula CDFs ---
  pop1.x.u  <- transform_to_gpd_cdf(Pgrid[,1], GPD_P1_var1)
  pop1.y.u  <- transform_to_gpd_cdf(Pgrid[,2], GPD_P1_var2)
  UU1       <- BiCopCDF(pop1.x.u, pop1.y.u, obj1)
  P1 <- 1 - pop1.x.u - pop1.y.u + UU1
  # old version
  # AEP.pop.1 <- matrix(P1 / (1/Rate_Con1), nrow=length(var1))
  
  # Poison assumption
  AEP.pop.1 <-matrix(1 - exp(-Rate_Con1 * P1),nrow=length(var1))
  
  
  pop2.x.u  <- transform_to_gpd_cdf(Pgrid[,1], GPD_P2_var1)
  pop2.y.u  <- transform_to_gpd_cdf(Pgrid[,2], GPD_P2_var2)
  UU4       <- BiCopCDF(pop2.x.u, pop2.y.u, obj4)
  P2 <- 1 - pop2.x.u - pop2.y.u + UU4
  # old version
  #AEP.pop.2 <- matrix(P2 / (1/Rate_Con4), nrow=length(var1))
  
  # Poison assumption
  AEP.pop.2 <-matrix(1 - exp(-Rate_Con4*P2),nrow=length(var1))
  
  
  
  #AEP.pop.1[AEP.pop.1 > 1] <- 1
  #AEP.pop.2[AEP.pop.2 > 1] <- 1
  RP_Grid <- 1 - (1 - AEP.pop.1) * (1 - AEP.pop.2)
  

  # --- Isolines & design events ---
  for(k in 1:length(RP)){
    iso          <- contourLines(var1, var2, RP_Grid, levels=RP[k])
    Isoline[[k]] <- data.frame(as.numeric(unlist(iso[[1]][2])), as.numeric(unlist(iso[[1]][3])))
    colnames(Isoline[[k]]) <- c(names(Data)[1], names(Data)[2])
    Iso <- Isoline[[k]]
    
    remove <- which(cop.sample[,1] > Sim_Max*max(Data[,1], na.rm=TRUE) | cop.sample[,2] > Sim_Max*max(Data[,2], na.rm=TRUE))
    if(length(remove) > 1) cop.sample <- cop.sample[-remove,]
    prediction <- kde(x=cop.sample, eval.points=Iso)$estimate
    Contour[[k]] <- (prediction - min(prediction)) / (max(prediction) - min(prediction))
    
    MostLikelyEvent.AND <- data.frame(as.numeric(Iso[which(prediction==max(prediction, na.rm=TRUE)), 1]),
                                      as.numeric(Iso[which(prediction==max(prediction, na.rm=TRUE)), 2]))
    colnames(MostLikelyEvent.AND) <- c(names(Data)[1], names(Data)[2])
    MostLikelyEvent[[k]] <- MostLikelyEvent.AND
    
    FullDependence.AND <- data.frame(max(Iso[,1]), max(Iso[,2]))
    colnames(FullDependence.AND) <- c(names(Data)[1], names(Data)[2])
    FullDependence[[k]] <- FullDependence.AND
    
    idx <- which(prediction > 0)
    #ensemble.AND <- Iso[sample(idx, size=N_Ensemble, replace=TRUE, prob=prediction[idx]),]  # <<< 'sample' now unambiguous
    
    ######################## added
    
    # Arc-length parameterization
    dx <- diff(Iso[,1]); dy <- diff(Iso[,2])
    seg_len <- c(0, cumsum(sqrt(dx^2 + dy^2)))
    
    # Interpolate prediction onto seg_len points (all isoline points, not just idx)
    pred_weights <- pmax(prediction, 0)  # ensure non-negative
    pred_weights <- pred_weights / sum(pred_weights)  # normalize
    
    # Sample continuous positions weighted by prediction
    t_samples <- sample(seg_len[-length(seg_len)], size=N_Ensemble, replace=TRUE, 
                        prob=pred_weights[-length(pred_weights)])  # match lengths
    
    t_samples <- t_samples + runif(N_Ensemble, 0, min(diff(seg_len)))
    
    # Interpolate coordinates
    ensemble.AND <- data.frame(
      Var1 = approx(seg_len, Iso[,1], xout=t_samples, rule=2)$y,
      Var2 = approx(seg_len, Iso[,2], xout=t_samples, rule=2)$y
    )
    colnames(ensemble.AND) <- c(names(Data)[1], names(Data)[2])
    
    colnames(ensemble.AND) <- c(names(Data)[1], names(Data)[2])
    Ensemble[[k]] <- data.frame(ensemble.AND)
  }
  
  # --- Plot ---
  x_min <- ifelse(is.na(x_lim_min), min(na.omit(Data[,con1])), x_lim_min)
  x_max <- ifelse(is.na(x_lim_max), max(na.omit(Data[,con1])), x_lim_max)
  y_min <- ifelse(is.na(y_lim_min), min(na.omit(Data[,con2])), y_lim_min)
  y_max <- ifelse(is.na(y_lim_max), max(na.omit(Data[,con2])), y_lim_max)
  
  par(mfrow=c(1,1))
  plot(Data[,1], Data[,2], xlim=c(x_min,x_max), ylim=c(y_min,y_max), col="Light Grey", xlab=x_lab, ylab=y_lab, cex.lab=1.5, cex.axis=1.5)
  points(Data_Con1[,1], Data_Con1[,2], col="Red",  pch=4, cex=1.5)
  points(Data_Con4[,1], Data_Con4[,2], col="Blue", pch=1, cex=1.5)
  for(k in 1:length(RP)){
    points(Isoline[[k]][,1], Isoline[[k]][,2], col=rev(heat.colors(150))[20:120][1+100*Contour[[k]]], lwd=3, pch=16, cex=1.75)
    if(N_Ensemble > 0) points(Ensemble[[k]][,1], Ensemble[[k]][,2], col=1, lwd=2, pch=16, cex=1)
    points(MostLikelyEvent[[k]][,1], MostLikelyEvent[[k]][,2], pch=18, cex=1.75)
    text(MostLikelyEvent[[k]][,1], MostLikelyEvent[[k]][,2], paste(RP[k]), col="White", cex=0.5)
    points(FullDependence[[k]][,1], FullDependence[[k]][,2], pch=17, cex=1.75)
    text(FullDependence[[k]][,1], FullDependence[[k]][,2], paste(RP[k]), col="White", cex=0.5)
  }
  
  # --- Output ---
  res <- list("GPD_P1_var1"=GPD_P1_var1,"GPD_P1_var2"=GPD_P1_var2,"GPD_P2_var1"=GPD_P2_var1,"GPD_P2_var2"=GPD_P2_var2,"AEP_P1"=AEP.pop.1, "AEP_P2"=AEP.pop.2, "Grid"=Pgrid, "FullDependence"=FullDependence, "MostLikelyEvent"=MostLikelyEvent, "Ensemble"=Ensemble, "Isoline"=Isoline, "Contour"=Contour)
  return(res)
}