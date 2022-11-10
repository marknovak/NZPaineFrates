#########################################################################
#########################################################################
#########################################################################
# Calculate and plot correlations
# and
# plot feeding rates versus prey abundances
#########################################################################
#########################################################################
rm(list = ls())
options(stringsAsFactors = F)

library(MASS)
source('PearsonApprox.R') # load PearsonApprox() & cor.test.ratio()
library(sfsmisc) # for eaxis
library(dplyr)
library(tidyr)
library(Hmisc) # for LaTeX table export
options(xdvicmd='open')
library(xtable) # for LaTeX summary table export

############################
data <- read.csv('../data/derived/NZ-1969_2004-tab_Summarized.csv')

############################
# Feeding rates
# f_i = n_i/n * 1/h_i
data$fi <- data$pi /data$h.mean

# Attack rates
# (only able to calculate for species that also showed up in abundance surveys)
# a_i = n_i/n_0 * 1/(h_i*N_i)
data$ai <- data$pi0 / (data$h.mean * data$N.mean)

sum(!is.na(data$fi))
sum(!is.na(data$ai))

############################
# Time-periods side-by-side and drop rows without feeding rate estimates
# (including Not Feeding)
dat <- data %>%
  pivot_wider(id_cols = c('Site','Prey'),
              names_from = 'Year',
              values_from = c('n.obs','h.mean',
                              'pi','pi0',
                              'N.mean',
                              'fi','ai')) %>%
  drop_na('fi_1969', 'fi_2004')

nrow(dat)
unique(dat$Prey)
############################
# Define function to place asterisk(s) for "significance"
signif <- function(x){
  p <- x$p.value
  out <-ifelse(p < 0.001, '***',
        ifelse(p < 0.01, '**',
        ifelse(p < 0.05, '*',
        ifelse(p < 0.1, '',
                 '   '))))
  return(out)
}

# Define function to calculate Mean Deviation and Mean Absolute Deviation
# as well as overall variable mean
# (on both natural and log10 scale)
deviation <- function(x){
  y <- x[,2]
  x <- x[,1]
  out <- c(mean = mean(c(x, y)),
              MD = mean(x - y),
              MAD = mean(abs(x - y)),
              mean.log10 = mean(log10(c(x,y))),
              MLD = mean(log10(x/y)),
              MALD = mean(abs(log10(x/y))),
              tenMLD = 10^mean(log10(x/y)),
              tenMALD = 10^mean(abs(log10(x/y))))
  # print(out)
  return(out)
}

boot.deviation <- function(x, n = 10000){
  boots <- replicate(n, deviation(x[sample(1:nrow(x), replace = TRUE),]))
  out <- cbind(
    mean = apply(boots, 1, mean),
    t(apply(boots, 1, quantile, probs = c(0.025, 0.975)))
  )
  print(out)
  return(out)
}

############################
pdf('../figs/Paine-Comparisons.pdf',
    width = 17/2.54,
    height = 8.5/2.54)
par(
  mfrow = c(2, 3),
  pty = 's',
  cex = 0.7,
  cex.axis = 0.9,
  cex.lab = 1,
  tcl = -0.2,
  mar = c(3, 4, 1, 1),
  mgp = c(2, 0.3, 0),
  las = 1
)
#~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~
# All feeding survey sites
#~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~
xymag <- c(0.5, 1.5) # reduce/increase xylims
#~~~~~~~~~~~~~~
# Feeding rates
#~~~~~~~~~~~~~~
f.r <- cor.test.ratios(dat$pi_1969, dat$h.mean_1969,
                       dat$pi_2004, dat$h.mean_2004)
flog.r <- cor.test.ratios(dat$pi_1969, dat$h.mean_1969,
                          dat$pi_2004, dat$h.mean_2004,
                          log10.transf = TRUE)
f.rs <- cor.test.ratios(dat$pi_1969, dat$h.mean_1969,
                        dat$pi_2004, dat$h.mean_2004,
                        method = 'spearman')

f.dev <- boot.deviation(cbind(dat$fi_1969, dat$fi_2004))

f.range <- range(c(dat$fi_1969, dat$fi_2004), na.rm = TRUE)
xylim <-  f.range * xymag

matplot(
  spread(dplyr::select(dat, Site, Prey, fi_1969), Site, fi_1969)[,-1],
  spread(dplyr::select(dat, Site, Prey, fi_2004), Site, fi_2004)[,-1],
  pch = 25:21,
  col = 'black',
  bg = 'grey',
  log = 'xy',
  xlab = '',
  ylab = '2004 feeding rate',
  axes = FALSE,
  xlim = xylim,
  ylim = xylim
)
title(xlab = '1968-9 feeding rate', line = 1.4)
abline(0, 1, lty = 2, col = 'grey70')
ats <- 10 ^ seq(-6, -1)
eaxis(1, at = ats)
eaxis(2, at = ats)
box(lwd = 1)
legend(
  'topleft',
  legend = c('Leigh (ER)', 'Leigh (TR)'),
  pch = 25:24,
  col = 'black',
  pt.bg = 'grey',
  inset = 0.0,
  cex = 0.7,
  bty = 'n'
)
legend(
  'bottomright',
  legend = c(
    as.expression(bquote(italic(r) == 
                           paste(.(round(f.r$analytic$estimate, 2)), 
                                 .(signif(f.r$analytic))))), 
    as.expression(bquote(italic(r)[10] == 
                           paste(.(round(flog.r$analytic$estimate, 2)), 
                                 .(signif(flog.r$analytic))))),
    as.expression(bquote(italic(r)[s] == 
                           paste(.(round(f.rs$analytic$estimate, 2)), 
                                 .(signif(f.rs$analytic)))))
    ),
  bty = 'n',
  inset = 0.0,
  y.intersp = 1,
  cex = 0.7
)
mtext('(a)', 3, adj = 0, line = 0.1, cex = 0.7)

#~~~~~~~~~~~~~~~~~
# Diet proportions
#~~~~~~~~~~~~~~~~~
fp.r <- cor.test(dat$pi_1969, dat$pi_2004, 
                     use = 'complete.obs')
fplog.r <- cor.test(log10(dat$pi_1969), log10(dat$pi_2004),
                        use = 'complete.obs')
fp.rs <- cor.test(dat$pi_1969, dat$pi_2004, 
                       use = 'complete.obs', method = 'spearman')

fp.dev <- boot.deviation(cbind(dat$pi_1969, dat$pi_2004))

fp.range <- range(c(dat$pi_1969, dat$pi_2004), na.rm = TRUE)

fp.mean <- apply(cbind(dat$pi_1969, dat$pi_2004), 2, 
                function(x){mean(x)})
fp.sd <- apply(cbind(dat$pi_1969, dat$pi_2004), 2, 
              function(x){sd(x)})
fp.cv <- apply(cbind(dat$pi_1969, dat$pi_2004), 2, 
              function(x){sd(x)/mean(x)})

xylim <-  fp.range * xymag

matplot(
  spread(dplyr::select(dat, Site, Prey, pi_1969), Site, pi_1969)[,-1],
  spread(dplyr::select(dat, Site, Prey, pi_2004), Site, pi_2004)[,-1],
  pch = 25:21,
  col = 'black',
  bg = 'grey',
  log = 'xy',
  xlab = '',
  ylab = '2004 diet proportion',
  axes = FALSE,
  xlim = xylim,
  ylim = xylim
)
title(xlab = '1968-9 diet proportion', line = 1.4)
abline(0, 1, lty = 2, col = 'grey70')
ats <- 10 ^ seq(-6, -1)
eaxis(1, at = ats)
eaxis(2, at = ats)
box(lwd = 1)
legend(
  'topleft',
  legend = rev(c('Leigh (WR)', 'Rangitoto', 'Whangaparaoa')),
  pch = rev(23:21),
  col = 'black',
  pt.bg = 'grey',
  inset = 0.0,
  cex = 0.7,
  bty = 'n'
)
legend(
  'bottomright',
  legend = c(
    as.expression(bquote(italic(r) == 
                           paste(.(round(fp.r$estimate, 2)), 
                                 .(signif(fp.r))))), 
    as.expression(bquote(italic(r)[10] == 
                           paste(.(round(fplog.r$estimate, 2)), 
                                 .(signif(fplog.r))))),
    as.expression(bquote(italic(r)[s] == 
                           paste(.(round(fp.rs$estimate, 2)), 
                                 .(signif(fp.rs)))))
  ),
  bty = 'n',
  inset = 0.0,
  y.intersp = 1,
  cex = 0.7
)
mtext('(b)', 3, adj = 0, line = 0.1, cex = 0.7)


#~~~~~~~~~~~~~~~~~~~~~~~~~
# Detection/handling times
#~~~~~~~~~~~~~~~~~~~~~~~~~
h.r <- cor.test(dat$h.mean_1969, dat$h.mean_2004, 
           use = 'complete.obs')
hlog.r <- cor.test(log10(dat$h.mean_1969), log10(dat$h.mean_2004), 
                       use = 'complete.obs')
h.rs <- cor.test(dat$h.mean_1969, dat$h.mean_2004, 
                      use = 'complete.obs', method = 'spearman')

h.dev <- boot.deviation(cbind(dat$h.mean_1969, dat$h.mean_2004))

h.range <- range(c(dat$h.mean_1969, dat$h.mean_2004), na.rm = TRUE)

h.range * 24 # Range in hours
h.mean <- apply(cbind(dat$h.mean_1969, dat$h.mean_2004), 2, 
              function(x){mean(x)})
h.sd <- apply(cbind(dat$h.mean_1969, dat$h.mean_2004), 2, 
              function(x){sd(x)})
h.cv <- apply(cbind(dat$h.mean_1969, dat$h.mean_2004), 2, 
              function(x){sd(x)/mean(x)})

xylim <- h.range * xymag

matplot(
  spread(dplyr::select(dat, Site, Prey, h.mean_1969), Site, h.mean_1969)[,-1],
  spread(dplyr::select(dat, Site, Prey, h.mean_2004), Site, h.mean_2004)[,-1],
  pch = 25:21,
  col = 'black',
  bg = 'grey',
  log = 'xy',
  xlab = '',
  ylab = '2004 detection time',
  axes = FALSE,
  xlim = xylim,
  ylim = xylim
)
title(xlab = '1968-9 detection time', line = 1.4)
abline(0, 1, lty = 2, col = 'grey70')
ats <- 10 ^ seq(-1, 4)
eaxis(1, at = ats)
eaxis(2, at = ats)
box(lwd = 1)
legend(
  'bottomright',
  legend = c(
      as.expression(bquote(italic(r) == 
                             paste(.(round(h.r$estimate, 2)), 
                                   .(signif(h.r))))), 
      as.expression(bquote(italic(r)[10] == 
                             paste(.(round(hlog.r$estimate, 2)), 
                                   .(signif(hlog.r))))),
      as.expression(bquote(italic(r)[s] == 
                             paste(.(round(h.rs$estimate, 2)), 
                                   .(signif(h.rs)))))
    ),
  bty = 'n',
  inset = 0.0,
  y.intersp = 1,
  cex = 0.7
)
mtext('(c)', 3, adj = 0, line = 0.1, cex = 0.7)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# All feeding survey AND abundance survey sites
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
datsub <- dat %>%  drop_na('N.mean_1969', 'N.mean_2004')

#~~~~~~~~~~~~~~
# Feeding rates
#~~~~~~~~~~~~~~
f.r.sub <- cor.test.ratios(datsub$pi_1969, datsub$h.mean_1969,
                           datsub$pi_2004, datsub$h.mean_2004)
flog.r.sub <- cor.test.ratios(datsub$pi_1969, datsub$h.mean_1969,
                              datsub$pi_2004, datsub$h.mean_2004,
                              log10.transf = TRUE)
f.rs.sub <- cor.test.ratios(datsub$pi_1969, datsub$h.mean_1969,
                            datsub$pi_2004, datsub$h.mean_2004,
                            method = 'spearman')

f.dev.sub <- boot.deviation(cbind(datsub$fi_1969, datsub$fi_2004))

f.range.sub <- range(c(datsub$fi_1969, datsub$fi_2004), na.rm = TRUE)
xylim <- f.range.sub * xymag

matplot(
  spread(dplyr::select(datsub, Site, Prey, fi_1969), Site, fi_1969)[,-1],
  spread(dplyr::select(datsub, Site, Prey, fi_2004), Site, fi_2004)[,-1],
  pch = 23:21,
  col = 'black',
  bg = 'grey',
  log = 'xy',
  xlab = '',
  ylab = '2004 feeding rate',
  axes = FALSE,
  xlim = xylim,
  ylim = xylim
)
title(xlab = '1968-9 feeding rate', line = 1.4)
abline(0, 1, lty = 2, col = 'grey70')
ats <- 10 ^ seq(-6, -1)
eaxis(1, at = ats)
eaxis(2, at = ats)
box(lwd = 1)
legend(
  'bottomright',
  legend = c(
    as.expression(bquote(italic(r) == 
                           paste(.(round(f.r.sub$analytic$estimate, 2)), 
                                 .(signif(f.r.sub$analytic))))), 
    as.expression(bquote(italic(r)[10] == 
                           paste(.(round(flog.r.sub$analytic$estimate, 2)), 
                                 .(signif(flog.r.sub$analytic))))),
    as.expression(bquote(italic(r)[s] == 
                           paste(.(round(f.rs.sub$analytic$estimate, 2)), 
                                 .(signif(f.rs.sub$analytic)))))
  ),
  bty = 'n',
  inset = 0.0,
  y.intersp = 1,
  cex = 0.7
)
mtext('(d)', 3, adj = 0, line = 0.1, cex = 0.7)


#~~~~~~~~~~~
# Proportions
#~~~~~~~~~~~
fp.r.sub <- cor.test(datsub$pi_1969, datsub$pi_2004, 
                 use = 'complete.obs')
fplog.r.sub <- cor.test(log10(datsub$pi_1969), log10(datsub$pi_2004),
                    use = 'complete.obs')
fp.rs.sub <- cor.test(datsub$pi_1969, datsub$pi_2004, 
                  use = 'complete.obs', method = 'spearman')

fp.dev.sub <- boot.deviation(cbind(datsub$pi_1969, datsub$pi_2004))

fp.range.sub <- range(c(datsub$pi_1969, datsub$pi_2004), na.rm = TRUE)

#~~~~~~~~~~~
# Detection times
#~~~~~~~~~~~
h.r.sub <- cor.test(datsub$h.mean_1969, datsub$h.mean_2004, 
                use = 'complete.obs')
hlog.r.sub <- cor.test(log10(datsub$h.mean_1969), log10(datsub$h.mean_2004), 
                   use = 'complete.obs')
h.rs.sub <- cor.test(datsub$h.mean_1969, datsub$h.mean_2004, 
                 use = 'complete.obs', method = 'spearman')

h.dev.sub <- boot.deviation(cbind(datsub$h.mean_1969, datsub$h.mean_2004))

h.range.sub <- range(c(datsub$h.mean_1969, datsub$h.mean_2004), na.rm = TRUE)


#~~~~~~~~~~~
# Abundances
#~~~~~~~~~~~
N.r.sub <-  cor.test(datsub$N.mean_1969, datsub$N.mean_2004,
                 use = 'complete.obs')
Nlog.r.sub <- cor.test(log10(datsub$N.mean_1969), log10(datsub$N.mean_2004),
                   use = 'complete.obs')
N.rs.sub <- cor.test(datsub$N.mean_1969, datsub$N.mean_2004,
                 use = 'complete.obs', method = 'spearman')

N.dev.sub <- boot.deviation(cbind(datsub$N.mean_1969, datsub$N.mean_2004))

N.range.sub <- range(c(datsub$N.mean_1969, datsub$N.mean_2004), na.rm = TRUE)

xylim <- N.range.sub * xymag

matplot(
  spread(dplyr::select(datsub, Site, Prey, N.mean_1969), Site, N.mean_1969)[,-1],
  spread(dplyr::select(datsub, Site, Prey, N.mean_2004), Site, N.mean_2004)[,-1],
  pch = 23:21,
  col = 'black',
  bg = 'grey',
  log = 'xy',
  xlab = '',
  ylab = '2004 abundance',
  axes = FALSE,
  xlim = xylim,
  ylim = xylim
)
title(xlab = '1968-9 abundance', line = 1.4)
abline(0, 1, lty = 2, col = 'grey70')
ats <- 10 ^ seq(-6, 4)
eaxis(1, at = ats)
eaxis(2, at = ats)
box(lwd = 1)
legend(
  'bottomright',
  legend = c(
    as.expression(bquote(italic(r) == 
                           paste(.(round(N.r.sub$estimate, 2)), 
                                 .(signif(N.r.sub))))), 
    as.expression(bquote(italic(r)[10] == 
                           paste(.(round(Nlog.r.sub$estimate, 2)), 
                                 .(signif(Nlog.r.sub))))),
    as.expression(bquote(italic(r)[s] == 
                           paste(.(round(N.rs.sub$estimate, 2)), 
                                 .(signif(N.rs.sub)))))
  ),
  bty = 'n',
  inset = 0.0,
  y.intersp = 1,
  cex = 0.7
)
mtext('(e)', 3, adj = 0, line = 0.1, cex = 0.7)


#~~~~~~~~~~~~~
# Attack rates
#~~~~~~~~~~~~~
a.r.sub <- cor.test(datsub$ai_1969, datsub$ai_2004,
                use = 'complete.obs')
alog.r.sub <- cor.test(log10(datsub$ai_1969), log10(datsub$ai_2004), 
           use = 'complete.obs')
a.rs.sub <- cor.test(datsub$ai_1969, datsub$ai_2004, 
                 use = 'complete.obs', method = 'spearman')

a.dev.sub <- boot.deviation(cbind(datsub$ai_1969, datsub$ai_2004))

a.range.sub <- range(c(datsub$ai_1969, datsub$ai_2004), na.rm = TRUE)
xylim <- a.range.sub * xymag

matplot(
  spread(dplyr::select(datsub, Site, Prey, ai_1969), Site, ai_1969)[,-1],
  spread(dplyr::select(datsub, Site, Prey, ai_2004), Site, ai_2004)[,-1],
  pch = 23:21,
  col = 'black',
  bg = 'grey',
  log = 'xy',
  xlab = '',
  ylab = '2004 attack rate',
  axes = FALSE,
  xlim = xylim,
  ylim = xylim
)
title(xlab = '1968-9 attack rate', line = 1.4)
abline(0, 1, lty = 2, col = 'grey70')
ats <- 10 ^ seq(-6, -1)
eaxis(1, at = ats)
eaxis(2, at = ats)
box(lwd = 1)
legend(
  'bottomright',
  legend = c(
    as.expression(bquote(italic(r) == 
                           paste(.(round(a.r.sub$estimate, 2)), 
                                 .(signif(a.r.sub))))), 
    as.expression(bquote(italic(r)[10] == 
                           paste(.(round(alog.r.sub$estimate, 2)), 
                                 .(signif(alog.r.sub))))),
    as.expression(bquote(italic(r)[s] == 
                           paste(.(round(a.rs.sub$estimate, 2)), 
                                 .(signif(a.rs.sub)))))
  ),
  bty = 'n',
  inset = 0.0,
  y.intersp = 1,
  cex = 0.7
)
mtext('(f)', 3, adj = 0, line = 0.1, cex = 0.7)
dev.off()



##########################################################
# Compare relationships between feeding rate and abundance
##########################################################
datFN <- data %>%
  drop_na(N.mean, fi)
# These data include 5 unpaired observations
# (i.e. where only one of us observed a given prey species)

regFNYearInt <- lm(log10(fi) ~ log10(N.mean) * Year, data = datFN)
  summary(regFNYearInt) # Year & intxn not signif
regFNYear <- lm(log10(fi) ~ log10(N.mean) + Year, data = datFN)
  summary(regFNYear) # Year not signif
regFN <- lm(log10(fi) ~ log10(N.mean), data = datFN)
  summary(regFN) # Intercept and slope signif


latex(
  xtable(summary(regFNYearInt)),
  file='../tables/Paine-FNYearInt.tex',
  # rowname = NULL,
  rowlabel = '',
  label = 'tab:FNYint',
  # center = 'centering',
  first.hline.double = FALSE,
  caption="Summary table for the regression of 
  prey-specific feeding rate 
  on prey-specific abundance,
  time period (\\emph{Year}), 
  and their interaction.",
  digits=3,
  where="!htbp"
)
latex(
  xtable(summary(regFNYear)),
  file='../tables/Paine-FNYear.tex',
  # rowname = NULL,
  rowlabel = '',
  label = 'tab:FNY',
  # center = 'centering',
  first.hline.double = FALSE,
  caption="Summary table for the regression of 
  prey-specific feeding rate 
  on prey-specific abundance and time period (\\emph{Year}).",
  digits=3,
  where="!htbp"
)
latex(
  xtable(summary(regFN)),
  file='../tables/Paine-FN.tex',
  # rowname = NULL,
  rowlabel = '',
  label = 'tab:FN',
  # center = 'centering',
  first.hline.double = FALSE,
  caption="Summary table for the regression of 
  prey-specific feeding rate 
  on prey-specific abundance.",
  digits=3,
  where="!htbp"
)

xlim <- range(datFN$N[!is.na(datFN$fi)], na.rm = TRUE)
ylim <- range(datFN$fi, na.rm = TRUE)

pdf('../figs/Paine-FuncResp.pdf',
    width = 8.5/2.54,
    height = 8.5*0.75/2.54)
par(
  cex = 0.8,
  cex.axis = 0.9,
  cex.lab = 1,
  tcl = -0.2,
  mar = c(3, 4, 1, 1),
  mgp = c(2, 0.3, 0),
  las = 1
)

  plot(datFN$N.mean, datFN$fi,
       log = 'xy',
       type = 'n',
       axes = FALSE,
       xlab = expression(paste('Prey abundance, ', italic(N[i]))),
       ylab = expression(paste('Feeding rate, ', italic(f[i]))))
  points(datFN$N[datFN$Year==1969], datFN$f[datFN$Year==1969],
         pch = 21,
         bg = 'grey30')
  points(datFN$N[datFN$Year==2004], datFN$f[datFN$Year==2004],
         pch = 23,
         bg = 'grey70')
  ats <- 10 ^ seq(-6, 5)
  eaxis(1, at = ats)
  eaxis(2, at = ats)
  box(lwd = 1)
  legend('topleft',
         legend = c('1968-9','2004'),
         pch = c(21, 23),
         pt.bg = c('grey30', 'grey70'),
         inset = 0,
         bty = 'n'
         )
  clip(xlim[1], xlim[2], ylim[1], ylim[2])
  abline(regFN, untf = F)
  
dev.off()



#########################################################################
#########################################################################
#########################################################################
