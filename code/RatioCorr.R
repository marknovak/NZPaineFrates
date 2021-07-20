#################################################
#################################################
# Spurious versus non-spurious ratio correlations
#################################################
#################################################
library(MASS) # for mvrnorm
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Multivariate-normal random variables
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

cor.obs <- cor(xn/xd, yn/yd)

num.sim <- 9999

# It's a spurious correlation when drawing inference about the numerators
# (randomizing both xn and yn is technically unnecessary; one would suffice)
cor.spur <- replicate(num.sim, cor(sample(xn) / xd, 
                                   sample(yn) / yd))
# It's *not* a spurious correlation when drawing inference about the ratios
cor.nonspur <- replicate(num.sim, cor(sample(xn) / sample(xd), 
                                      sample(yn) / sample(yd)))

cor.obs
mean(cor.spur)
mean(cor.nonspur)

hist(cor.nonspur, breaks = 200, xlim = c(-1,1))
hist(cor.spur, breaks = 200, xlim = c(-1,1), add = TRUE)
abline(v = cor.obs, col = 'blue')
abline(v = mean(cor.spurr), col = 'red')
legend('topleft',
       legend = c('Observed','Expected'),
       lty = 1,
       col = c('blue','red'))

# Conclusion: In the context of comparing feeding rates, 
#  we are drawing inference about the correlation of the ratios
#  (not the numerator diet proportions). 
# A correlation of zero, i.e. mean(cor.nonspur), 
#  is thus the appropriate null hypothesis.
# If we were drawing inference on the numerator diet proportions, 
#  then the non-zero correlation, i.e. mean(cor.spur),
#  is the appropriate null hypothesis.  


