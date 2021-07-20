######################################
####################################
# Define convenience functions
# to calculate correlations of ratios
######################################
######################################
library(MASS) # for mvrnorm

#~~~~~~~~~~~~~~~~~~~~~~~~
# Pearson's approximation 
#~~~~~~~~~~~~~~~~~~~~~~~~
# for expected value of cor(xn/xd, yn/yd)

PearsonApprox <- function(xn, xd, 
                          yn, yd){
  x <- xn # numerator of first ratio
  w <- xd # denominator of first ratio
  y <- yn # numerator of second ratio
  z <- yd # denominator of second ratio
  
  r.yx <- cor(y,x)
  r.yw <- cor(y,w)
  r.xz <- cor(x,z)
  r.zw <- cor(z,w)
  r.yz <- cor(y,z)
  r.xw <- cor(x,w)
  v.x <- sd(x)/mean(x)
  v.y <- sd(y)/mean(y)
  v.w <- sd(w)/mean(w)
  v.z <- sd(z)/mean(z)
  
  r.yz.xw <- (r.yx*v.y*v.x 
              - r.yw*v.y*v.w 
              - r.xz*v.x*v.z 
              + r.zw*v.z*v.w) /
    ( sqrt(v.y^2 + v.z^2 - 2*r.yz*v.y*v.z) * 
        sqrt(v.x^2 + v.w^2 - 2*r.xw*v.x*v.w) )
  return(r.yz.xw)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Permutation-based correlation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (rather than using asymptotic estimator)
# Adapted from perm.cor.test of jmuOutlier package
# https://rdrr.io/cran/jmuOutlier/src/R/perm.cor.test.R

cor.test.ratios <- function(xn, xd,
                            yn, yd,
                            alternative = "two.sided", 
                            method = "pearson", 
                            ci.probs = c(0.025, 0.975),
                            num.sim = 59999 ) {
  if (!(alternative %in% c("two.sided", "less", "greater")))
    stop("'alternative' must be 'two.sided', 'less', or 'greater'.")
  if (any( 
    c(length(xn), length(xd), length(yn), length(yd)) != length(xn)))
    stop("All vectors must have the same length.")
  if (method %in% c("pearson", "spearman")) {
    test.stat.obs <- cor(xn/xd, yn/yd, method = method)
    test.stat.PearsonApprox <- NULL
    if(method %in% 'pearson'){
      test.stat.PearsonApprox <- PearsonApprox(xn, xd, yn, yd)
    }
    # Note that only the numerators are shuffled.
    # If the denominators are also randomized, 
    # then the expected correlation is zero.
    test.stat <- replicate(num.sim, 
                           cor(sample(xn) / xd, sample(yn) / yd,
                               method = method))
    test.stat.mean <- mean(test.stat)
    conf.int <- quantile(test.stat, probs = ci.probs, na.rm = TRUE)
    if (alternative == "two.sided")
      p.value <- mean(abs(test.stat) >= abs(test.stat.obs))
    if (alternative == "less")
      p.value <- mean(test.stat <= test.stat.obs)
    if (alternative == "greater")
      p.value <- mean(test.stat >= test.stat.obs)
  }
  out <- list(
    summary=list(observed.statistic = test.stat.obs,
                 conf.int = conf.int,
                 statistic.mean = test.stat.mean, 
                 statistic.PearsonApprox = test.stat.PearsonApprox,
                 alternative = alternative, 
                 p.value = p.value,
                 method = method,
                 num.sim = num.sim),
    stat.vals = list(stat.vals = test.stat))
  print(out$summary)
  return(out)
}

