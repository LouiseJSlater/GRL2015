

#  Slater, L.J., Singer, M.B. and Kirchner, J.W., 2015. 
#  Hydrologic versus geomorphic drivers of trends in flood hazard. 
#  Geophysical Research Letters, 42(2), pp.370-376.

library(minpack.lm)

###################################
#                                 #
#        IRLS COEFFS              #
#                                 #
###################################

IRLS.exp.unbiasedmean <- function(Y, X, type="Cauchy") { 
  Cauchy <- function(x,MAR) { 
    if (MAR==0) w <- ifelse(x==0, 1, 0)
    else w <- 1/(1+(x/(3.536*MAR))^2) 
    return(w)
  }
  if (type!="Cauchy") stop("IRLS stopped: no valid weight type specified.  Only 'Cauchy' is implemented here!")
  wt <- rep(1,length(Y))
  delta_x <- X - mean(X)
  mean_y <- weighted.mean(Y, wt)
  fit <- nlsLM(Y/mean_y ~ exp(a*delta_x)/weighted.mean(exp(a*delta_x),wt), 
               start=list(a=0),  
               #start=list(a =0.0),  
               weights=wt, na.action="na.omit", control = nls.lm.control (maxiter=1000))
  wt_chg <- 999.0
  iter <- 0
  while ( (max(wt_chg,na.rm=TRUE) > 0.01) & !all(iter>10, summary(fit)$r.squared>0.999) 
  ){
    if (iter>1000) stop("IRLS stopped: more than 1000 interations, sorry!")   
    iter <- iter+1
    old_wt <- wt
    mean_y <- weighted.mean(Y,wt)   
    slope <- summary(fit)$coefficients[1,1]   
    resid <- Y - mean_y*(exp(slope*delta_x)/mean(exp(slope*delta_x)))    
    abs_resid <- abs(resid)
    abs_resid_nonzero <- ifelse(Y==0, NA, abs_resid)
    MAR <- median(abs_resid_nonzero, na.rm=TRUE)  
    if (MAR==0.0) stop("IRLS stopped. Solution has collapsed: median absolute residual is zero!")
    wt <- Cauchy(abs_resid,MAR)   
    fit <- nlsLM(Y/mean_y ~ exp(a*delta_x)/weighted.mean(exp(a*delta_x), wt), 
                 start=list(a=0),  
                 #start=list(a=0.0),
                 weights=wt, na.action="na.omit", control = nls.lm.control (maxiter=1000))      
    wt_chg <- abs(wt-old_wt)   
  } 
  return(fit)
}


# trends <- tapply(1:length(df$site), df$site, function (idx) { 
#      fit <- summary(IRLS.exp.unbiasedmean(
#      df$Qhat_extracted_RIs[idx], as.numeric(df$Date)[idx],type="Cauchy"))})


###################################
#                                 #
#        IRLS FITS                #
#                                 #
###################################

IRLS.fits <- function(Y, X, type="Cauchy") { 
  Cauchy <- function(x,MAR) { 
    if (MAR==0) w <- ifelse(x==0, 1, 0)
    else w <- 1/(1+(x/(3.536*MAR))^2) 
    return(w)
  }
  if (type!="Cauchy") stop("IRLS stopped: no valid weight type specified.  Only 'Cauchy' is implemented here!")
  wt <- rep(1,length(Y))
  delta_x <- X - mean(X)
  mean_y <- weighted.mean(Y, wt)
  fit <- nlsLM(Y/mean_y ~ exp(a*delta_x)/weighted.mean(exp(a*delta_x),wt), 
               start=list(a=0),  
               #start=list(a=0.0), 
               weights=wt, na.action="na.omit", control = nls.lm.control (maxiter=1000))
  wt_chg <- 999.0
  iter <- 0
  while ( (max(wt_chg,na.rm=TRUE) > 0.01) & !all(iter>10, summary(fit)$r.squared>0.999) 
  ){
    if (iter>1000) stop("IRLS stopped: more than 1000 interations, sorry!")   
    iter <- iter+1
    old_wt <- wt
    mean_y <- weighted.mean(Y,wt)   
    slope <- summary(fit)$coefficients[1,1]   
    resid <- Y - mean_y*(exp(slope*delta_x)/mean(exp(slope*delta_x)))    
    abs_resid <- abs(resid)
    abs_resid_nonzero <- ifelse(Y==0, NA, abs_resid)
    MAR <- median(abs_resid_nonzero, na.rm=TRUE)   
    if (MAR==0.0) stop("IRLS stopped. Solution has collapsed: median absolute residual is zero!")
    wt <- Cauchy(abs_resid,MAR)   
    fit <- nlsLM(Y/mean_y ~ exp(a*delta_x)/weighted.mean(exp(a*delta_x), wt), 
                 start=list(a=0),  
                 #start = list(a = 0.0), 
                 weights=wt, na.action="na.omit", control = nls.lm.control (maxiter=1000))      
    wt_chg <- abs(wt-old_wt)   
  } 
  qq <- list("fit"=fit, "wt"=wt)
  return(qq)
}


# fits <- tapply(1:length(df$site), df$site, function (idx) { 
#          fit <- fitted((IRLS.fits(df$Qhat_extracted_RIs[idx],as.numeric(df$Date)[idx],type="Cauchy"))$fit)*
#                 weighted.mean(df$Qhat_extracted_RIs[idx],IRLS.fits(df$Qhat_extracted_RIs[idx], 
#                                         as.numeric(df$Date)[idx],type="Cauchy")$wt)  })
# 
# fits <- tapply(1:length(site), site, function (idx) { 
#   fit <- fitted((IRLS.fits(y[idx], x[idx],type="Cauchy"))$fit)*
#   weighted.mean(y[idx],IRLS.fits(y[idx], x[idx],type="Cauchy")$wt)  })


###################################
#                                 #
#        FF TREND                 #
#                                 #
###################################

# coefficients 
FF_ExpoUnb <- function(Y,X) { 
  delta_x <- X-mean(X)
  mean_y <- mean(Y)
  fit <- nlsLM(Y/mean_y~exp(a*delta_x)/mean(exp(a*delta_x)), 
               start=list(a=0), 
               #start=list(a=0.0),
               control=nls.lm.control(maxiter=1000))
  return(fit)  
}


# Ftrend <- tapply(1:length(mergeFF$site), mergeFF$site, 
#                function(idx) { fit <- summary(FF_ExpoUnb(mergeFF$FF_POT[idx], mergeFF$Year[idx]))})
# FFeffect <-ldply(lapply(Ftrend, function(x)coef(x)), rbind) 


# list of fits for printing
FF_fits <- function(Y,X) { 
  delta_x <- X-mean(X)
  mean_y <- mean(Y)
  fit <- nlsLM(Y/mean_y~exp(a*delta_x)/mean(exp(a*delta_x)), 
               start=list(a=0), 
               #start=list(a=0.0), 
               control=nls.lm.control(maxiter=1000))
  qq <- list("fit"=fit, "wt"=mean_y)
  return(qq)
}

# fits <- tapply(1:length(mergeFF$site), mergeFF$site, 
#              function (idx) { fit <- fitted((FF_fits(mergeFF$Q_POT[idx], mergeFF$Year[idx]))$fit)*
#                FF_fits(mergeFF$Q_POT[idx], mergeFF$Year[idx])$wt  })



#Monte Carlo P values
MonteP <- function(Y,X) { 
  delta_x <- X-mean(X)
  mean_y <- mean(Y)
  fit <- nlsLM(Y/mean_y~exp(a*delta_x)/mean(exp(a*delta_x)), 
               start=list(a=0), 
               #start=list(a=0.0), 
               control=nls.lm.control(maxiter=1000))
  slope.obs <- summary(fit)$coefficients[1]
  #now the shuffled bit
  teststat <- rep(NA, 1000)
  slope <- rep(NA, 1000)
  for(i in 1:1000) {
    ySHUFFLE <- sample(Y)
    SHUFFLE.nls <- nlsLM(ySHUFFLE/mean_y~exp(a*delta_x)/mean(exp(a*delta_x)), 
                         start=list(a=0), 
                         #start=list(a=0.0), 
                         control=nls.lm.control(maxiter=1000))
    slope[i] <- summary(SHUFFLE.nls)$coefficients[1] 
  }   
  result <- sum(abs(slope)>=abs(slope.obs))/1000 # count the number of abs slopes that are greater than the abs measured one
  return(result)
}

# MonteP <- tapply(1:length(mergeFF$site), mergeFF$site, 
#               function(idx) { result <- MonteP(mergeFF$Q_POT[idx], mergeFF$Year[idx])})
# as.data.frame(MonteP)


