####################################
####################################
# Compare Pearson's approximation
# to permutation-based estimates
####################################
####################################
library(MASS) # for mvrnorm
source('PearsonApprox.R') # load functions

set.seed(2)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Multivariate-normal random variables
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# independent numerators 
# but correlated denominators
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

S <- 100 # Sample size

xn <- rnorm(S, 10, 0.5)
yn <- rnorm(S, 10, 0.5)
cov <- 0.9 # Covariance between denominator variables
d <- data.frame(mvrnorm(n = S, 
                        mu = c(xd = 10, yd = 10), 
                        Sigma = rbind(c(1, cov), c(cov, 1))))
xd <- d$xd
yd <- d$yd

pairs(cbind(xn, yn, xd, yd))
cor(cbind(xn, yn, xd, yd))

# Ratios are correlated (despite non-correlated numerators)
plot(xn/xd, yn/yd)
cor(xn/xd, yn/yd)

# Pearson's expected correlation
PearsonApprox(xn, xd, yn, yd)

# Permutation based (exact)
out <- cor.test.ratios(xn, xd, yn, yd)

# Plot distribution of correlations
hist(out$stat.vals$stat.vals, breaks = 200)
abline(v = out$summary$statistic.mean, col='green', lwd = 8)
abline(v = out$summary$statistic.PearsonApprox, col='red', lwd = 5)
legend('topright', 
       legend = c('Mean', "Pearson's approximation"),
       col = c('green','red'),
       lty = 1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Apply to feeding-rate-like variables
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# independent numerators (proportions observed feeding)
# but correlated denominators (detection times)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set.seed(2)

S <- 100 # Sample size

# Random proportions
p1 <- runif(S)
xn <- p1/sum(p1)
p2 <- runif(S)
yn <- p2/sum(p2)

# Log-normal correlated detection times
cov <- 0.9 # Covariance
d <- exp(
  data.frame(mvrnorm(n = S, 
                     mu = c(xd = 1, yd = 1), 
                     Sigma = rbind(c(1, cov), c(cov, 1)))))
xd <- d$xd
yd <- d$yd


pairs(cbind(xn, yn, xd, yd))
cor(cbind(xn, yn, xd, yd))

# Ratios are correlated
plot(xn/xd, yn/yd)
cor(xn/xd, yn/yd)

# Pearson's expected correlation
PearsonApprox(xn, xd, yn, yd)

# Permutation based (exact)
out <- cor.test.ratios(xn, xd, yn, yd)

# Plot distribution of correlations
hist(out$stat.vals$stat.vals, breaks = 200)
abline(v = out$summary$statistic.mean, col='green', lwd = 8)
abline(v = out$summary$statistic.PearsonApprox, col='red', lwd = 5)
legend('topleft', 
       legend = c('Mean', "Pearson's approximation"),
       col = c('green','red'),
       lty = 1)

####################################
