#################################################
#################################################
# Spurious versus non-spurious ratio correlations
#################################################
#################################################
require(MASS) # for mvrnorm
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Multivariate-normal random variables
# independent numerators 
# but correlated denominators
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
S <- 100 # Sample size
xn <- rnorm(S, 10, 0.5) # first numerator
yn <- rnorm(S, 10, 0.5) # second numerator
cov <- 0.9 # Covariance between denominator variables
d <- data.frame(mvrnorm(n = S, 
                        mu = c(xd = 10, yd = 10), 
                        Sigma = rbind(c(1, cov), c(cov, 1))))
xd <- d$xd # first denominator
yd <- d$yd # second denominator

# Only the denominators are correlated
pairs(cbind(xn, yn, xd, yd))
cor(cbind(xn, yn, xd, yd))

# The ratios are correlated
cor.obs <- cor(xn/xd, yn/yd)

# Contrast spurious versus non-spurious inferences be permuting either 
# just the numerators or both the numerators and denominators

num.sim <- 9999 # number of permutations

# It's a spurious correlation when drawing inference about the numerators
cor.spur <- replicate(num.sim, cor(sample(xn) / xd, 
                                         (yn) / yd,
                                   method = 'pearson'))
# But *not* a spurious correlation when drawing inference about the ratios
cor.nonspur <- replicate(num.sim, cor(sample(xn/xd), 
                                             yn/yd,
                                      method = 'pearson'))

cor.obs
mean(cor.spur)
mean(cor.nonspur)

# Inspect distributions of observed correlations
pdf('../figs/CorrelationOfRatios.pdf', height = 3, width = 5)
par(
  cex = 0.8,
  cex.axis = 0.9,
  cex.lab = 1,
  tcl = -0.2,
  mar = c(3, 4, 1, 1),
  mgp = c(2, 0.3, 0),
  las = 1,
  yaxs = 'i'
)
h1 <- hist(cor.nonspur, breaks = 200, xlim = c(-1,1), 
           main = '', xlab = 'Correlation')
  abline(v = cor.obs, col = 'blue', lwd = 2)
  abline(v = mean(cor.spur), col = 'red', lwd = 2)
h2 <- hist(cor.spur, breaks = 200, xlim = c(-1,1), 
           add = TRUE)

  legend('topleft',
         legend = c('Observed',"Mean (Expected)"),
         lty = 1,
         col = c('blue','red'),
         bty = 'n')
dev.off()
# Conclusion: In the context of comparing feeding rates, 
#  we are drawing inference about the correlation of the ratios
#  (not the numerator diet proportions). 
# A correlation of zero, mean(cor.nonspur), 
#  is thus the appropriate null hypothesis.
# If we were drawing inference on the numerator diet proportions, 
#  then the non-zero correlation, mean(cor.spur),
#  is the appropriate null hypothesis.  


