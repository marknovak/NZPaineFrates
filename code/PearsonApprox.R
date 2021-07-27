######################################
####################################
# Define convenience functions
# to calculate correlations of ratios
######################################
######################################
require(MASS) # for mvrnorm

#~~~~~~~~~~~~~~~~~~~~~~~~
# Pearson's approximation 
#~~~~~~~~~~~~~~~~~~~~~~~~
# for expected value of cor(xn/xd, yn/yd)

PearsonApprox <- function(xn, xd, 
                          yn, yd,
                          log10.transf = FALSE){
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
  
  out <- list(estimate = r.yz.xw,
              r = c(r.xnxd = r.xw,
                    r.ynyd = r.yz,
                    r.xnyd = r.xz,
                    r.ynxd = r.yw,
                    r.xnyn = r.yx,
                    r.xdyd = r.zw),
              cv = c(v.xn = v.x, 
                     v.xd = v.w,
                     v.yn = v.y,
                     v.yd = v.z))
  return(out)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Permutation-based correlation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Adapted from perm.cor.test of jmuOutlier package
# https://rdrr.io/cran/jmuOutlier/src/R/perm.cor.test.R

cor.test.ratios <- function(xn, xd,
                            yn, yd,
                            log10.transf = FALSE,
                            alternative = "two.sided", 
                            method = "pearson", 
                            ci.probs = c(0.025, 0.975),
                            num.sim = 99999,
                            keep.sim.vals = FALSE) {
  
  if (!(alternative %in% c("two.sided", "less", "greater")))
    stop("'alternative' must be 'two.sided', 'less', or 'greater'.")
  
  if (any( 
    c(length(xn), length(xd), length(yn), length(yd)) != length(xn)))
    stop("All vectors must have the same length.")
  
  cc <- complete.cases(xn, xd, yn, yd)
  xn <- xn[cc]
  xd <- xd[cc]
  yn <- yn[cc]
  yd <- yd[cc]
  n <- sum(cc)
  
  if (method %in% c("pearson", "spearman")) {
    if(log10.transf){
      test.stat.obs <- cor.test(log10(xn/xd), log10(yn/yd), method = method)}
    else{
      test.stat.obs <- cor.test(xn/xd, yn/yd, method = method)}
    
    pearson.approx <- NULL
    if(method %in% 'pearson' & !log10.transf){
      pearson.approx <- PearsonApprox(xn, xd, yn, yd)
    }
    # Note that only the numerators are shuffled.
    # If the denominators are also randomized, 
    # then the expected correlation is zero
    if(log10.transf){
      test.stat <- replicate(num.sim, 
                             cor(log10(sample(xn) / xd), log10(sample(yn) / yd),
                                 method = method))}
    else{
      test.stat <- replicate(num.sim, 
                             cor(sample(xn) / xd, sample(yn) / yd,
                                 method = method))}
    
    test.stat.mean <- mean(test.stat)
    conf.int <- quantile(test.stat, probs = ci.probs, na.rm = TRUE)
    
    if (alternative == "two.sided")
      p.value <- mean(abs(test.stat) >= abs(test.stat.obs$estimate))
    if (alternative == "less")
      p.value <- mean(test.stat <= test.stat.obs$estimate)
    if (alternative == "greater")
      p.value <- mean(test.stat >= test.stat.obs$estimate)
  }
  if(!keep.sim.vals){test.stat <- NA}
  
  out <- list(
    analytic = test.stat.obs,
    pearson.approx = pearson.approx,
    perm = list(estimate = test.stat.mean, 
                conf.int = conf.int,
                alternative = alternative, 
                p.value = p.value,
                permutations = num.sim,
                stat.vals = test.stat),
    method = list(sample.size = n,
                  method = method,
                  log10.transf = log10.transf))

  print(out)
  return(out)
}

