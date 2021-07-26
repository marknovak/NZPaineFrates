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
library(dplyr)
library(sfsmisc) # for eaxis
############################

dat <- read.csv('../data/derived/NZ-1969_2004-tab_Summarized.csv')

############################
# Feeding rates
# f_i=n_i/n * 1/h_i
dat$fi <- dat$pi /dat$h.mean

# Attack rates
# (only able to calculate for species that also showed up in abundance surveys)
# a_i= n_i/n_0 * 1/(h_i*N_i)
dat$ai <- dat$pi0 / (dat$h.mean * dat$N.mean)

sum(!is.na(dat$fi))
sum(!is.na(dat$ai))

############################
# Time-periods side-by-side and drop rows without feeding rate estimates
# (including Not Feeding)
dat <- pivot_wider(data = dat, 
            id_cols = c('Site','Prey'),
            names_from = 'Year',
            values_from = c('n.obs','h.mean','pi','pi0','N.mean','fi','ai')) %>%
       drop_na('fi_1969', 'fi_2004')

############################
# Define convenience function to place asterisk(s) for "significance"
signif <- function(x){
  p <- x$p.value
  out <-ifelse(p < 0.001, '***',
        ifelse(p < 0.01, '**',
        ifelse(p < 0.05, '*',
        ifelse(p < 0.1, '',
                 ' ns'))))
  return(out)
}

############################
pdf('../figs/Paine-Comparisons.pdf',
    width = 8,
    height = 4)
par(
  mfrow = c(2, 3),
  pty = 's',
  cex = 0.8,
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
xymag <- c(0.5, 1.5)
#~~~~~~~~~~~~~~
# Feeding rates
#~~~~~~~~~~~~~~
xylim <- range(c(dat$fi_1969, dat$fi_2004), na.rm = TRUE) * xymag
f.rho.all <-
  cor.test(dat$fi_1969, dat$fi_2004,
      use = 'complete.obs')
flog.rho.all <-
  cor.test(log10(dat$fi_1969), log10(dat$fi_2004), 
      use = 'complete.obs')
f.rhos.all <-
  cor.test(dat$fi_1969, dat$fi_2004,
      use = 'complete.obs', method = 'spearman')
matplot(
  spread(select(dat, Site, Prey, fi_1969), Site, fi_1969)[,-1],
  spread(select(dat, Site, Prey, fi_2004), Site, fi_2004)[,-1],
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
                           paste(.(round(f.rho.all$estimate, 2)), 
                                 .(signif(f.rho.all))))), 
    as.expression(bquote(italic(r)[10] == 
                           paste(.(round(flog.rho.all$estimate, 2)), 
                                 .(signif(flog.rho.all))))),
    as.expression(bquote(italic(r)[s] == 
                           paste(.(round(f.rhos.all$estimate, 2)), 
                                 .(signif(f.rhos.all)))))
    ),
  bty = 'n',
  inset = 0.0,
  y.intersp = 1,
  cex = 0.7
)

#~~~~~~~~~~~~~~~~~
# Diet proportions
#~~~~~~~~~~~~~~~~~
xylim <- range(c(dat$pi_1969, dat$pi_2004), na.rm = TRUE) * xymag
fp.rho.all <-
  cor.test(dat$pi_1969, dat$pi_2004, 
      use = 'complete.obs')
fplog.rho.all <-
  cor.test(log10(dat$pi_1969), log10(dat$pi_2004),
      use = 'complete.obs')
fp.rhos.all <-
  cor.test(dat$pi_1969, dat$pi_2004, 
      use = 'complete.obs', method = 'spearman')
matplot(
  spread(select(dat, Site, Prey, pi_1969), Site, pi_1969)[,-1],
  spread(select(dat, Site, Prey, pi_2004), Site, pi_2004)[,-1],
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
                           paste(.(round(fp.rho.all$estimate, 2)), 
                                 .(signif(fp.rho.all))))), 
    as.expression(bquote(italic(r)[10] == 
                           paste(.(round(fplog.rho.all$estimate, 2)), 
                                 .(signif(fplog.rho.all))))),
    as.expression(bquote(italic(r)[s] == 
                           paste(.(round(fp.rhos.all$estimate, 2)), 
                                 .(signif(fp.rhos.all)))))
  ),
  bty = 'n',
  inset = 0.0,
  y.intersp = 1,
  cex = 0.7
)

#~~~~~~~~~~~~~~~~~~~~~~~~~
# Detection/handling times
#~~~~~~~~~~~~~~~~~~~~~~~~~
xylim <- range(c(dat$h.mean_1969, dat$h.mean_2004), na.rm = TRUE) * xymag
h.rho.all <-
  cor.test(dat$h.mean_1969, dat$h.mean_2004, 
      use = 'complete.obs')
hlog.rho.all <-
  cor.test(log10(dat$h.mean_1969), log10(dat$h.mean_2004), 
      use = 'complete.obs')
h.rhos.all <-
  cor.test(dat$h.mean_1969, dat$h.mean_2004, 
      use = 'complete.obs', method = 'spearman')
matplot(
  spread(select(dat, Site, Prey, h.mean_1969), Site, h.mean_1969)[,-1],
  spread(select(dat, Site, Prey, h.mean_2004), Site, h.mean_2004)[,-1],
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
                             paste(.(round(h.rho.all$estimate, 2)), 
                                   .(signif(h.rho.all))))), 
      as.expression(bquote(italic(r)[10] == 
                             paste(.(round(hlog.rho.all$estimate, 2)), 
                                   .(signif(hlog.rho.all))))),
      as.expression(bquote(italic(r)[s] == 
                             paste(.(round(h.rhos.all$estimate, 2)), 
                                   .(signif(h.rhos.all)))))
    ),
  bty = 'n',
  inset = 0.0,
  y.intersp = 1,
  cex = 0.7
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# All feeding survey AND abundance survey sites
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
datsub <- dat %>%  drop_na('N.mean_1969', 'N.mean_2004')

#~~~~~~~~~~~~~~
# Feeding rates
#~~~~~~~~~~~~~~
xylim <- range(c(datsub$fi_1969, datsub$fi_2004), na.rm = TRUE) * xymag
f.rho <-
  cor.test(datsub$fi_1969, datsub$fi_2004, 
      use = 'complete.obs')
flog.rho <-
  cor.test(log10(datsub$fi_1969), log10(datsub$fi_2004),
      use = 'complete.obs')
f.rhos <-
  cor.test(datsub$fi_1969, datsub$fi_2004, 
      use = 'complete.obs', method = 'spearman')
matplot(
  spread(select(datsub, Site, Prey, fi_1969), Site, fi_1969)[,-1],
  spread(select(datsub, Site, Prey, fi_2004), Site, fi_2004)[,-1],
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
                           paste(.(round(f.rho$estimate, 2)), 
                                 .(signif(f.rho))))), 
    as.expression(bquote(italic(r)[10] == 
                           paste(.(round(flog.rho$estimate, 2)), 
                                 .(signif(flog.rho))))),
    as.expression(bquote(italic(r)[s] == 
                           paste(.(round(f.rhos$estimate, 2)), 
                                 .(signif(f.rhos)))))
  ),
  bty = 'n',
  inset = 0.0,
  y.intersp = 1,
  cex = 0.7
)

#~~~~~~~~~~~
# Abundances
#~~~~~~~~~~~
xylim <- range(c(datsub$N.mean_1969, datsub$N.mean_2004), 
               na.rm = TRUE)  * xymag
N.rho <-
  cor.test(datsub$N.mean_1969, datsub$N.mean_2004, 
      use = 'complete.obs')
Nlog.rho <-
  cor.test(log10(datsub$N.mean_1969), log10(datsub$N.mean_2004),
      use = 'complete.obs')
N.rhos <-
  cor.test(datsub$N.mean_1969, datsub$N.mean_2004,
      use = 'complete.obs', method = 'spearman')
matplot(
  spread(select(datsub, Site, Prey, N.mean_1969), Site, N.mean_1969)[,-1],
  spread(select(datsub, Site, Prey, N.mean_2004), Site, N.mean_2004)[,-1],
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
                           paste(.(round(N.rho$estimate, 2)), 
                                 .(signif(N.rho))))), 
    as.expression(bquote(italic(r)[10] == 
                           paste(.(round(Nlog.rho$estimate, 2)), 
                                 .(signif(Nlog.rho))))),
    as.expression(bquote(italic(r)[s] == 
                           paste(.(round(N.rhos$estimate, 2)), 
                                 .(signif(N.rhos)))))
  ),
  bty = 'n',
  inset = 0.0,
  y.intersp = 1,
  cex = 0.7
)

#~~~~~~~~~~~~~
# Attack rates
#~~~~~~~~~~~~~
xylim <- range(c(datsub$ai_1969, datsub$ai_2004), na.rm = TRUE) * xymag
a.rho <-
  cor.test(datsub$ai_1969, datsub$ai_2004,
      use = 'complete.obs')
alog.rho <-
  cor.test(log10(datsub$ai_1969), log10(datsub$ai_2004), 
      use = 'complete.obs')
a.rhos <-
  cor.test(datsub$ai_1969, datsub$ai_2004, 
      use = 'complete.obs', method = 'spearman')
matplot(
  spread(select(datsub, Site, Prey, ai_1969), Site, ai_1969)[,-1],
  spread(select(datsub, Site, Prey, ai_2004), Site, ai_2004)[,-1],
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
                           paste(.(round(a.rho$estimate, 2)), 
                                 .(signif(a.rho))))), 
    as.expression(bquote(italic(r)[10] == 
                           paste(.(round(alog.rho$estimate, 2)), 
                                 .(signif(alog.rho))))),
    as.expression(bquote(italic(r)[s] == 
                           paste(.(round(a.rhos$estimate, 2)), 
                                 .(signif(a.rhos)))))
  ),
  bty = 'n',
  inset = 0.0,
  y.intersp = 1,
  cex = 0.7
)

dev.off()



##########################################################
# Compare relationships between feeding rate and abundance
##########################################################

# Separate models
reg69 <- lm(log10(datsub$fi_1969) ~ log10(datsub$N.mean_1969))
reg04 <- lm(log10(datsub$fi_2004) ~ log10(datsub$N.mean_2004))

# Single model
datr <- data.frame(
  f = c(datsub$fi_1969, datsub$fi_2004),
  N = c(datsub$N.mean_1969, datsub$N.mean_2004),
  Year = rep.int(c('1969', '2004'), 
                 c(nrow(datsub), nrow(datsub))),
  Prey = datsub$Prey
)


regFNYearInt <- lm(log10(f) ~ log10(N) * Year, data = datr)
summary(regFNYearInt) # Year & intxn not signif
regFNYear <- lm(log10(f) ~ log10(N) + Year, data = datr)
summary(regFNYear) # Year not signif
regFN <- lm(log10(f) ~ log10(N), data = datr)
summary(regFN) # Intercept and slope signif

latex(
  xtable(summary(regFNYearInt)),
  file='../tables/Paine-FNYearInt.tex',
  # rowname = NULL,
  rowlabel = '',
  label = 'tab:FNYint',
  # center = 'centering',
  first.hline.double = FALSE,
  caption="Summary table for the regression of prey-specific feeding rate 
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
  caption="Summary table for the regression of prey-specific feeding rate 
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
  caption="Summary table for the regression of prey-specific feeding rate 
  on prey-specific abundance.",
  digits=3,
  where="!htbp"
)

xlim <- range(datr$N[!is.na(datr$f)], na.rm = TRUE)
ylim <- range(datr$f, na.rm = TRUE)

pdf('../figs/Paine-FuncResp.pdf',
    width = 4,
    height = 3)
par(
  cex = 0.8,
  cex.axis = 0.9,
  cex.lab = 1,
  tcl = -0.2,
  mar = c(3, 4, 1, 1),
  mgp = c(2, 0.3, 0),
  las = 1
)

  plot(datr$N, datr$f,
       log = 'xy',
       type = 'n',
       axes = FALSE,
       xlab = expression(paste('Abundance, ', italic(N[i]))),
       ylab = expression(paste('Feeding rate, ', italic(f[i]))))
  points(datr$N[datr$Year==1969], datr$f[datr$Year==1969],
         pch = 21,
         bg = 'grey30')
  points(datr$N[datr$Year==2004], datr$f[datr$Year==2004],
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