library(MASS) # for mvrnorm
library(boot) # for boot(strap)

set.seed(2)

S <- 10
cov <- 0.8 # Covariance between vectors of detection times



# # Generate random proportions
# p1 <- runif(S)
#   p1 <- p1/sum(p1)
# p2 <- runif(S)
#   p2 <- p2/sum(p2)
# 
# # Generate log-normal distributed correlated data
# mu <- c(d1 = 1, d2 = 1)
# sigma<-rbind(c(1, cov), c(cov, 1))
# d <- exp(mvrnorm(n = S, mu = mu, Sigma = sigma))


p1 <- rnorm(S, 10, 0.5)
p2 <- rnorm(S, 10, 0.5)
d <- rnorm(S,10,0.5)
d <- cbind(d1=d, d2=d)

d1 <- dat$d1
d2 <- dat$d2



dat <- data.frame(p1, p2, d)
pairs(dat)
cor(dat)

f <- data.frame(f1 = dat$p1 / dat$d1,
                f2 = dat$p2 / dat$d2)

plot(f)
plot(f, log = 'xy')
cor.test(f$f1, f$f2)

###############
# Pearson's approximation for expected value of cor(xn/xd, yn/yd)
PearApprox <- function(xn, xd, 
                       yn, yd){
  x <- xn
  w <- xd
  y <- yn
  z <- yd

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

PearApprox(xn = dat$p1,
           xd = dat$d1, 
           yn = dat$p2, 
           yd = dat$d2)


###############
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
    test.stat.exp <- NULL
    if(method %in% 'pearson'){
      test.stat.exp <- PearApprox(xn, xd, yn, yd)
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
                 statistic.exp = test.stat.exp,
                 alternative = alternative, 
                 p.value = p.value,
                 method = method,
                 num.sim = num.sim),
            stat.vals = list(stat.vals = test.stat))
  print(out$summary)
  return(out)
}

####################################

cor.test(f$f1, f$f2)
cor.test(dat$p1/dat$d1, dat$p2/dat$d2)
out <- cor.test.ratios(x = dat$p1, 
                       xd = dat$d1, 
                       yn = dat$p2, 
                       yd = dat$d2,
                       method = 'pearson')



hist(out$stat.vals$stat.vals, breaks=200)
abline(v=out$summary$statistic.mean, col='green')
abline(v=out$summary$statistic.exp, col='red')



cor.obs <- cor(p1 / d1, p2 / d2)
cor.exp <- PearApprox(p1, d1, p2, d2)

num.sim <- 9999
cor.spurr <- replicate(num.sim, cor(sample(p1) / d1, 
                                    sample(p2) / d2))
cor.nonspurr <- replicate(num.sim, cor(sample(p1) / sample(d1), 
                                       sample(p2) / sample(d2)))


cor.obs
mean(cor.spurr)
mean(cor.nonspurr)

hist(cor.nonspurr, breaks = 200, xlim = c(-1,1), col = 'grey')
hist(cor.spurr, breaks = 200, xlim = c(-1,1))
  abline(v = mean(cor.spurr), col = 'red')
  abline(v = cor.obs, col = 'blue')
  abline(v = cor.exp, col = 'green')
