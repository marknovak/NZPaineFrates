#################################################
#################################################
# Spurious versus non-spurious ratio correlations
#################################################
#################################################
require(MASS) # for mvrnorm
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simple contrived example(s)
xn <- c(0.3, 0.2, 0.4) # first numerator (diet proportions)
yn <- c(0.2, 0.3, 0.1) # second numerator (diet proportions)
xd <- c(2.5, 1, 5) # first denominator (detection times)
yd <- c(1, 0.5, 3.5) # second denominator (detection times)

xr <- xn/xd
yr <- yn/yd

corrs = round(c(n = cor(xn, yn),
                d = cor(xd, yd),
                r = cor(xr,yr)),
             3)

cv <- function(x){sd(x)/mean(x)} # coefficient of variation
cvs <- round(c(xn = cv(xn),
               yn = cv(yn),
               xd = cv(xd),
               yd = cv(yd)),
             3)

Spp <- paste('Prey', c(1,2,3))

pdf('../figs/Example.pdf', width = 7.5, height = 2.5)
par(mfrow = c(1,3), pty = 's', las = 1, mar = c(3.25, 3.25, 1.5, 0.1),
    cex.lab = 1.25, tcl = -0.2, mgp = c(2, 0.3, 0))

plot(xn, yn, 
     pch = 21, bg = 'grey', 
     xlim = c(0, 0.5), ylim = c(0, 0.5),
     xlab = 'Survey 1',
     ylab = 'Survey 2')
  title(main = 'Diet proportions')
  box(lwd = 1)
  text(xn, yn, Spp,
       adj = c(0.9,-1))
  legend('topright',
         bty = 'n',
         cex = 1.1,
         legend = c(as.expression(bquote(italic(r) == .(corrs['n']))),
                    as.expression(bquote(italic(cv[1]) == .(cvs['xn']))),
                    as.expression(bquote(italic(cv[2]) == .(cvs['yn'])))))
  mtext('A', 3, adj = 0)

plot(xd, yd, 
     pch = 21, bg = 'grey', 
     xlim = c(0, 6), ylim = c(0, 5),
     xlab = 'Survey 1',
     ylab = 'Survey 2')
  title(main = 'Detection times')
  text(xd, yd, Spp,
       adj = c(0.8,-1))
  legend('topleft',
         bty = 'n',
         cex = 1.1,
         legend = c(as.expression(bquote(italic(r) == .(corrs['d']))),
                    as.expression(bquote(italic(cv[1]) == .(cvs['xd']))),
                    as.expression(bquote(italic(cv[2]) == .(cvs['yd'])))))
  mtext('B', 3, adj = 0)

plot(xr, yr, 
     pch = 21, bg = 'grey', 
     xlim = c(0.05, 0.25), ylim = c(0, 0.7),
     xlab = 'Survey 1',
     ylab = 'Survey 2')
  title(main = 'Feeding rates')
  box(lwd = 1)
  text(xr, yr, Spp,
       adj = c(0,-1))
  legend('topleft',
         inset = 0,
         bty = 'n',
         cex = 1.1,
         legend = c(as.expression(bquote(italic(r) == .(corrs['r'])))))
  mtext('C', 3, adj = 0)
  
dev.off()

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


