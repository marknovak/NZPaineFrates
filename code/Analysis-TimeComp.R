#########################################################################
#########################################################################
#########################################################################
rm(list = ls())
options(stringsAsFactors = F)
library(dplyr)
library(sfsmisc) # for eaxis
library(binom) # for confidence intervals
library(Hmisc) # for LaTeX table export
  options(xdvicmd='open')
library(xtable) # for LaTeX summary table export
############################

# All data at hand
focalSites <- read.csv("../data/derived/NZ-1969_2004-Site-focal.csv")$x

sizes <- 
  read.csv("../data/derived/NZ-1969_2004-tab_Sizes.csv")

h_all_69 <-
  as.matrix(read.csv("../data/derived/NZ-1969_2004-tab_Htime_all_69.csv",
                     row.names = 1))
n_all_69 <-
  as.matrix(read.csv("../data/derived/NZ-1969_2004-tab_Fobs_all_69.csv",
                     row.names = 1))
N_all_69 <-
  as.matrix(read.csv("../data/derived/NZ-1969_2004-tab_Abund_all_69.csv",
                     row.names = 1))
h_all_04 <-
  as.matrix(read.csv("../data/derived/NZ-1969_2004-tab_Htime_all_04.csv",
                     row.names = 1))
n_all_04 <-
  as.matrix(read.csv("../data/derived/NZ-1969_2004-tab_Fobs_all_04.csv",
                     row.names = 1))
N_all_04 <-
  as.matrix(read.csv("../data/derived/NZ-1969_2004-tab_Abund_all_04.csv",
                     row.names = 1))

# Sites with only feeding surveys
h_compf_69 <-
  as.matrix(read.csv("../data/derived/NZ-1969_2004-tab_Htime_compf_69.csv", 
                     row.names = 1))
n_compf_69 <-
  as.matrix(read.csv("../data/derived/NZ-1969_2004-tab_Fobs_compf_69.csv",
                     row.names = 1))

h_compf_04 <-
  as.matrix(read.csv("../data/derived/NZ-1969_2004-tab_Htime_compf_04.csv", 
                     row.names = 1))
n_compf_04 <-
  as.matrix(read.csv("../data/derived/NZ-1969_2004-tab_Fobs_compf_04.csv", 
                     row.names = 1))

# Sites with both feeding and abundance surveys
h_comp_69 <-
  as.matrix(read.csv("../data/derived/NZ-1969_2004-tab_Htime_comp_69.csv", 
                     row.names = 1))
n_comp_69 <-
  as.matrix(read.csv("../data/derived/NZ-1969_2004-tab_Fobs_comp_69.csv", 
                     row.names = 1))
N_comp_69 <-
  as.matrix(read.csv("../data/derived/NZ-1969_2004-tab_Abund_comp_69.csv", 
                     row.names = 1))

h_comp_04 <-
  as.matrix(read.csv("../data/derived/NZ-1969_2004-tab_Htime_comp_04.csv", 
                     row.names = 1))
n_comp_04 <-
  as.matrix(read.csv("../data/derived/NZ-1969_2004-tab_Fobs_comp_04.csv", 
                     row.names = 1))
N_comp_04 <-
  as.matrix(read.csv("../data/derived/NZ-1969_2004-tab_Abund_comp_04.csv", 
                     row.names = 1))

############################
# Restrict abundance data to those species that occurred in diets
N_comp_69[is.na(h_comp_69)] <- NA
N_comp_04[is.na(h_comp_04)] <- NA

############################
# Feeding proportions
fp_compf_69 <-
  (t(t(n_compf_69) / colSums(n_compf_69, na.rm = TRUE)))[-1, ]
fp_compf_04 <-
  (t(t(n_compf_04) / colSums(n_compf_04, na.rm = TRUE)))[-1, ]

fp_comp_69 <- (t(t(n_comp_69) / colSums(n_comp_69, na.rm = TRUE)))[-1, ]
fp_comp_04 <- (t(t(n_comp_04) / colSums(n_comp_04, na.rm = TRUE)))[-1, ]

# feeding rates
# (able to calculate for all species)
# f_i=n_i/n * 1/h_i
f_compf_69 <-
  (t(t(n_compf_69) / colSums(n_compf_69, na.rm = TRUE)) / h_compf_69)[-1, ]
f_compf_04 <-
  (t(t(n_compf_04) / colSums(n_compf_04, na.rm = TRUE)) / h_compf_04)[-1, ]

f_comp_69 <-
  (t(t(n_comp_69) / colSums(n_comp_69, na.rm = TRUE)) / h_comp_69)[-1, ]
f_comp_04 <-
  (t(t(n_comp_04) / colSums(n_comp_04, na.rm = TRUE)) / h_comp_04)[-1, ]

# attack rates
# (only able to calculate for species that also showed up in abundance surveys)
# a_i= n_i/n_0 * 1/(h_i*N_i)
a_comp_69 <-
  (t(t(n_comp_69) / n_comp_69[1, ]) / (h_comp_69 * N_comp_69))[-1, ]
a_comp_04 <-
  (t(t(n_comp_04) / n_comp_04[1, ]) / (h_comp_04 * N_comp_04))[-1, ]

a_comp_04[is.infinite(a_comp_04)] <- NA


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

# feeding rate - feeding survey only sites
xylim <- range(c(f_compf_69, f_compf_04), na.rm = TRUE)
f.rho.all <-
  cor.test(as.vector(f_compf_69), as.vector(f_compf_04),
      use = 'complete.obs')
flog.rho.all <-
  cor.test(log10(as.vector(f_compf_69)), log10(as.vector(f_compf_04)), 
      use = 'complete.obs')
f.rhos.all <-
  cor.test(as.vector(f_compf_69), as.vector(f_compf_04),
      use = 'complete.obs', method = 'spearman')
matplot(
  f_compf_69,
  f_compf_04,
  pch = 21:25,
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
  legend = rev(c('Leigh (ER)', 'Leigh (TR)')),
  pch = rev(24:25),
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

# Diet proportion - feeding survey only sites
xylim <- range(c(fp_compf_69, fp_compf_04), na.rm = TRUE)
fp.rho.all <-
  cor.test(as.vector(fp_compf_69), as.vector(fp_compf_04), 
      use = 'complete.obs')
fplog.rho.all <-
  cor.test(log10(as.vector(fp_compf_69)), log10(as.vector(fp_compf_04)),
      use = 'complete.obs')
fp.rhos.all <-
  cor.test(as.vector(fp_compf_69), as.vector(fp_compf_04), 
      use = 'complete.obs', method = 'spearman')
matplot(
  fp_compf_69,
  fp_compf_04,
  pch = 21:25,
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
  pch = rev(21:23),
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

# detection/handling time - feeding survey only sites
xylim <- range(c(h_compf_69, h_compf_04), na.rm = TRUE)
h.rho.all <-
  cor.test(as.vector(h_compf_69), as.vector(h_compf_04), 
      use = 'complete.obs')
hlog.rho.all <-
  cor.test(log10(as.vector(h_compf_69)), log10(as.vector(h_compf_04)), 
      use = 'complete.obs')
h.rhos.all <-
  cor.test(as.vector(h_compf_69), as.vector(h_compf_04), 
      use = 'complete.obs', method = 'spearman')
matplot(
  h_compf_69,
  h_compf_04,
  pch = 21:25,
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

# feeding rate - feeding and abundance survey sites
xylim <- range(c(f_comp_69, f_comp_04), na.rm = TRUE)
f.rho <-
  cor.test(as.vector(f_comp_69), as.vector(f_comp_04), 
      use = 'complete.obs')
flog.rho <-
  cor.test(log10(as.vector(f_comp_69)), log10(as.vector(f_comp_04)), 
      use = 'complete.obs')
f.rhos <-
  cor.test(as.vector(f_comp_69), as.vector(f_comp_04), 
      use = 'complete.obs', method = 'spearman')
matplot(
  f_comp_69,
  f_comp_04,
  pch = 21:23,
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

# abundance
xylim <- range(c(N_comp_69, N_comp_04), na.rm = TRUE) + 1E-1
N.rho <-
  cor.test(as.vector(N_comp_69), as.vector(N_comp_04), 
      use = 'complete.obs')
Nlog.rho <-
  cor.test(log10(as.vector(N_comp_69)), log10(as.vector(N_comp_04)), 
      use = 'complete.obs')
N.rhos <-
  cor.test(as.vector(N_comp_69), as.vector(N_comp_04), 
      use = 'complete.obs', method = 'spearman')
matplot(
  N_comp_69,
  N_comp_04,
  pch = 21:23,
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

# attack rate
xylim <- range(c(a_comp_69, a_comp_04), na.rm = TRUE)
a.rho <-
  cor.test(as.vector(a_comp_69), as.vector(a_comp_04), 
      use = 'complete.obs')
alog.rho <-
  cor.test(log10(as.vector(a_comp_69)), log10(as.vector(a_comp_04)), 
      use = 'complete.obs')
a.rhos <-
  cor.test(as.vector(a_comp_69), as.vector(a_comp_04), 
      use = 'complete.obs', method = 'spearman')
matplot(
  a_comp_69,
  a_comp_04,
  pch = 21:23,
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

#########################
# Compare relationships between feeding rate and abundance
#########################

# Plaxiphora's abundance is 0, so can't log
N_comp_04[N_comp_04 == 0] <- NA

# Separate models
reg69 <- lm(log10(as.vector(f_comp_69)) ~ log10(as.vector(N_comp_69[-1, ])))
reg04 <- lm(log10(as.vector(f_comp_04)) ~ log10(as.vector(N_comp_04[-1, ])))

# Single model
dat <- data.frame(
  f = c(as.vector(f_comp_69), as.vector(f_comp_04)),
  N = c(as.vector(N_comp_69[-1, ]), as.vector(N_comp_04[-1, ])),
  Year = rep.int(c('1969', '2004'), c(length(
    as.vector(f_comp_69)
  ), length(
    as.vector(f_comp_04)
  ))),
  Species = c(rownames(f_comp_69), rownames(f_comp_04))
)


regFNYearInt <- lm(log10(f) ~ log10(N) * Year, data = dat)
summary(regFNYearInt) # Year & intxn not signif
regFNYear <- lm(log10(f) ~ log10(N) + Year, data = dat)
summary(regFNYear) # Year not signif
regFN <- lm(log10(f) ~ log10(N), data = dat)
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

xlim <- range(dat$N[!is.na(dat$f)], na.rm = TRUE)
ylim <- range(dat$f, na.rm = TRUE)

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

  plot(dat$N, dat$f,
       log = 'xy',
       type = 'n',
       axes = FALSE,
       xlab = expression(paste('Abundance, ', italic(N[i]))),
       ylab = expression(paste('Feeding rate, ', italic(f[i]))))
  points(dat$N[dat$Year==1969], dat$f[dat$Year==1969],
         pch = 21,
         bg = 'grey30')
  points(dat$N[dat$Year==2004], dat$f[dat$Year==2004],
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

##############################
# Pred and prey sizes
##############################
# Size distributions
sizes <- sizes[sizes$Site%in%focalSites,]

preysizes <- sizes[sizes$PreySize > 0 ,]
breaks = seq(0, max(sizes$PreySize) + 1, 1)
xlims <- range(breaks)
Prey <- hist(preysizes$PreySize, breaks = breaks, 
             plot = FALSE)
Prey04 <- hist(preysizes$PreySize[preysizes$Year==2004], breaks = breaks,
               plot = FALSE)

breaks = seq(0, max(sizes$PredSize) + 1, 1)
ylims <- range(breaks)
Pred <- hist(sizes$PredSize, breaks = breaks, plot = FALSE)
Pred04 <- hist(sizes$PredSize[sizes$Year==2004], breaks = breaks, 
               plot = FALSE)

top <- c(max(Prey$counts), max(Pred$counts))

# KS tests on distributions
ks.prey <- ks.test(preysizes$PreySize[preysizes$Year==1969],
                   preysizes$PreySize[preysizes$Year==2004])
lapply(list(preysizes$PreySize[preysizes$Year==1969],
            preysizes$PreySize[preysizes$Year==2004]), mean)

ks.pred <- ks.test(sizes$PredSize[sizes$Year==1969],
                   sizes$PredSize[sizes$Year==2004])
lapply(list(sizes$PredSize[sizes$Year==1969],
            sizes$PredSize[sizes$Year==2004]), mean)

# Pred-prey size selectivity
regPredPreyYearInt <- lm(log(PredSize) ~ log(PreySize) * Year, data = preysizes)
summary(regPredPreyYearInt)
regPredPreyYear <- lm(log(PredSize) ~ log(PreySize) + Year, data = preysizes)
summary(regPredPreyYear)
regPredPrey <- lm(log(PredSize) ~ log(PreySize), data = preysizes)
summary(regPredPrey)


latex(
  xtable(summary(regPredPreyYearInt)),
  file='../tables/Paine-PredPreySizeYearInt.tex',
  # rowname = NULL,
  rowlabel = '',
  label = 'tab:SizeYint',
  # center = 'centering',
  first.hline.double = FALSE,
  caption="Summary table for the regression of predator size on prey size,
  time period (\\emph{Year}), and their interaction.",
  digits = 3,
  where = "!htbp"
)
latex(
  xtable(summary(regPredPreyYear)),
  file='../tables/Paine-PredPreySizeYear.tex',
  # rowname = NULL,
  rowlabel = '',
  label = 'tab:SizeY',
  # center = 'centering',
  first.hline.double = FALSE,
  caption="Summary table for the regression of predator size on prey size and 
  time period (\\emph{Year}).",
  digits = 3,
  where = "!htbp"
)
latex(
  xtable(summary(regPredPrey)),
  file='../tables/Paine-PredPreySize.tex',
  # rowname = NULL,
  rowlabel = '',
  label = 'tab:Size',
  # center = 'centering',
  first.hline.double = FALSE,
  caption="Summary table for the regression of predator size 
  on prey size.",
  digits = 3,
  where = "!htbp"
)


regPredPrey69 <- lm(log(PredSize) ~ log(PreySize), 
                  data = preysizes[preysizes$Year==1969,])
summary(regPredPrey69)
rng69 <- range(preysizes[preysizes$Year==1969,]$PreySize)

regPredPrey04 <- lm(log(PredSize) ~ log(PreySize), 
                  data = preysizes[preysizes$Year==2004,])
summary(regPredPrey04)
rng04 <- range(preysizes[preysizes$Year==2004,]$PreySize)

pdf('../figs/Paine-Sizes.pdf',
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

nf <- layout(matrix(c(2,0,1,3),
                    2,2,
                    byrow = TRUE), 
             c(7,1), c(1,5), TRUE)
# layout.show(nf)

par(mar = c(3,3,0,0))
  plot(sizes$PreySize, sizes$PredSize,
       type = 'n',
       xlab = 'Prey size (mm)',
       ylab = 'Predator size (mm)',
       xlim = xlims, 
       ylim = ylims
  )
  abline(0,1, col = 'grey70', lty = 2)
  lat = rng69[1]:rng69[2]
  lines(lat,
        exp(predict(regPredPrey69, data.frame('PreySize'=lat))),
        col = 'grey40',
        lwd = 3)
  lines(lat,
        exp(predict(regPredPrey69, data.frame('PreySize'=lat))),
        col = 'white',
        lwd = 1)
  lat = rng04[1]:rng04[2]
  lines(lat,
        exp(predict(regPredPrey04, data.frame('PreySize'=lat))),
        col= 'grey60',
        lwd = 3)
  lines(lat,
        exp(predict(regPredPrey04, data.frame('PreySize'=lat))),
        col= 'white',
        lwd = 1)
  points(sizes$PreySize[sizes$Year==2004], 
         sizes$PredSize[sizes$Year==2004],
         pch = 21,
         bg = 'grey70',
         cex = 0.7)
  points(sizes$PreySize[sizes$Year==1969], 
         sizes$PredSize[sizes$Year==1969],
         pch = 23,
         bg = 'grey30',
         cex = 0.7)
  legend('bottomright',
         legend = c('1968-9','2004'),
         pch = c(21, 23),
         pt.bg = c('grey30','grey70'),
         inset = 0.01)
  box(lwd = 1)

  
par(mar = c(0,3,1,0))
  barplot(Prey$counts, axes = FALSE, ylim = c(0, top[1]), space = 0,
          col = 'grey30')
  axis(2, cex.axis = 0.5)
  barplot(Prey04$counts, axes = FALSE, ylim = c(0, top[1]), space = 0,
          col = 'grey70', add = TRUE)
  
par(mar = c(3,0,0,1))
  barplot(Pred$counts, axes = FALSE, xlim = c(0, top[2]), space = 0, 
          horiz = TRUE,
          col = 'grey30')
  axis(1, cex.axis = 0.5, mgp = c(0, 0, 0))
  barplot(Pred04$counts, axes = FALSE, xlim = c(0, top[2]), space = 0, 
          horiz = TRUE,
          col = 'grey70', add = TRUE)
  
dev.off()

#################################
# Summary table
###############
# Total count of observations (feeding + not feeding)
nTot_69 <- apply(n_compf_69, 2, sum, na.rm = TRUE)
nTot_04 <- apply(n_compf_04, 2, sum, na.rm = TRUE)

# Total count of observations (feeding)
nTotF_69 <- apply(n_compf_69[-1,], 2, sum, na.rm = TRUE)
nTotF_04 <- apply(n_compf_04[-1,], 2, sum, na.rm = TRUE)

# Fraction feeding overall
fF_69 <- binom.confint(nTot_69-n_compf_69['Not Feeding',], nTot_69, 
                       methods = 'wilson')[,-1]
fF_04 <- binom.confint(nTot_04-n_compf_04['Not Feeding',], nTot_04, 
                       methods = 'wilson')[,-1]
rownames(fF_69) <- rownames(fF_04) <- names(nTot_69)

# Fraction feeding on H. scobina of those observed feeding
fHs_69 <- binom.confint(n_compf_69['Haustrum scobina',], nTotF_69, 
                       methods = 'wilson')[,-1]
n_compf_04['Haustrum scobina',][is.na(n_compf_04['Haustrum scobina',])] <- 0
fHs_04 <- binom.confint(n_compf_04['Haustrum scobina',], nTotF_04, 
                       methods = 'wilson')[,-1]
rownames(fHs_69) <- rownames(fHs_04) <- names(nTot_69)

fF_69 <- fF_69[order(rownames(fF_69)),]
fF_04 <- fF_04[order(rownames(fF_04)),]
fHs_69 <- fHs_69[order(rownames(fHs_69)),]
fHs_04 <- fHs_04[order(rownames(fHs_04)),]

size.stats <- sizes %>% group_by(Site, Year) %>%
  summarise(mean = round(mean(PredSize, na.rm = TRUE),1),
            max = max(PredSize, na.rm = TRUE))
size.stats <- as.data.frame(size.stats)
size.stats <- size.stats[size.stats$Site%in%focalSites,]

s69 <- data.frame(nT = nTot_69[ord],
                fF = 100*round(fF_69[,-c(1,2)],3),
                fHs = 100*round(fHs_69[,-c(1,2)],3),
                size.stats[size.stats$Year==1969,-c(1,2)])
s04 <- data.frame(nT = nTot_04[ord],
                     fF = 100*round(fF_04[,-c(1,2)],3),
                     fHs = 100*round(fHs_04[,-c(1,2)],3),
                     size.stats[size.stats$Year==2004,-c(1,2)])

summ <- data.frame(
   Site = rownames(s04),
   s69$nT,
   s04$nT,
   fF69 = paste0(s69$fF.mean,' (',s69$fF.lower,'-',s69$fF.upper,')'),
   fF04 = paste0(s04$fF.mean,' (',s04$fF.lower,'-',s04$fF.upper,')'),
   fHs69 = paste0(s69$fHs.mean,' (',s69$fHs.lower,'-',s69$fHs.upper,')'),
   fHs04 = paste0(s04$fHs.mean,' (',s04$fHs.lower,'-',s04$fHs.upper,')')
)

sum_fF_69 <- apply(s69, 2, sum)
sum_fF_04 <- apply(s04, 2, sum)
mean_fF_69 <- apply(s69, 2, mean)
mean_fF_04 <- apply(s04, 2, mean)

summ <- rbind(summ, c('Sum/Average',
                      sum_fF_69['nT'], sum_fF_04['nT'],
                      mean_fF_69['fF.mean'], mean_fF_04['fF.mean'],
                      mean_fF_69['fHs.mean'], mean_fF_04['fHs.mean'],
                      mean_fF_69['mean'], mean_fF_04['mean'],
                      mean_fF_69['max'], mean_fF_04['max']))

summ$Site[which(summ$Site == "Leigh...Waterfall.Rocks")] <- "Leigh - Waterfall Rocks"
summ$Site[which(summ$Site == "LeighTabletopRocksandBoulders")] <- "Leigh - Tabletop Rocks"
summ$Site[which(summ$Site == "LeighEchinodermReef")] <- "Leigh - Echinoderm Reef"
summ$Site[which(summ$Site == "Rangitoto.Island...Whites.Beach")] <- "Rangitoto Island - Whites Beach"
summ$Site[which(summ$Site == "Red.Beach...Whangaparaoa")] <- "Red Beach - Whangaparaoa"

latex(
  summ,
  file='../tables/Paine-SiteSumm.tex',
  cgroup = c('', 'Observations', '\\% feeding', 
             '\\% feeding on \\emph{H. scobina}'),
  n.cgroup = c(1, 2, 2, 2),
  colheads = c('Site',
               '1968-9', '2004',
               '1968-9', '2004',
               '1968-9', '2004'),
  n.rgroup = c(5,1),
  rowname = NULL,
  label = 'tab:summ',
  center = 'centering',
  first.hline.double = FALSE,
  caption="Summary of Paine's 1968-9 and my 2004 feeding observations.  Observerations refers to the total number of whelks inspected. \\% feeding refers to the proportion of observed whelks that were feeding. \\% feeding on \\emph{H. scobina} refers to the proportion of feeding whelks that were feeding on \\emph{Haustrum scobina}.  Parentheticals are the biomial confidence interval (95\\% coverage probability) calculated using the Wilson method.",
)

#########################################################################
#########################################################################
#########################################################################