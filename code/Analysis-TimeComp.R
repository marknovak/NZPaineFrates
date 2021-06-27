###############################################################################
###############################################################################
###############################################################################
rm(list = ls())
options(stringsAsFactors = F)
library(sfsmisc) # for eaxis
library(binom) # for confidence intervals
library(Hmisc) # for LaTeX table export
  options(xdvicmd='open')
############################

# All data at hand
siteInfo <- read.csv("../data/derived/NZ-1969_2004-Site-LatLon.csv")

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
fR2all <-
  cor(as.vector(f_compf_69), as.vector(f_compf_04),
      use = 'complete.obs') ^ 2
flogR2all <-
  cor(log10(as.vector(f_compf_69)), log10(as.vector(f_compf_04)), 
      use = 'complete.obs') ^ 2
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
  legend = c(as.expression(bquote(R ^ 2 == .(
    round(fR2all, 2)
  ))), as.expression(bquote(R[10] ^ 2 == .(
    round(flogR2all, 2)
  )))),
  bty = 'n',
  inset = 0.0,
  y.intersp = 1,
  cex = 0.7
)

# Diet proportion - feeding survey only sites
xylim <- range(c(fp_compf_69, fp_compf_04), na.rm = TRUE)
fpR2all <-
  cor(as.vector(fp_compf_69), as.vector(fp_compf_04), 
      use = 'complete.obs') ^ 2
fplogR2all <-
  cor(log10(as.vector(fp_compf_69)), log10(as.vector(fp_compf_04)),
      use = 'complete.obs') ^ 2
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
  legend = c(as.expression(bquote(R ^ 2 == .(round(fpR2all, 3)))), 
             as.expression(bquote(R[10] ^ 2 == .(round(fplogR2all, 2))))),
  bty = 'n',
  inset = 0.0,
  y.intersp = 1,
  cex = 0.7
)

# handling time - feeding survey only sites
xylim <- range(c(h_compf_69, h_compf_04), na.rm = TRUE)
hR2all <-
  cor(as.vector(h_compf_69), as.vector(h_compf_04), 
      use = 'complete.obs') ^ 2
hlogR2all <-
  cor(log10(as.vector(h_compf_69)), log10(as.vector(h_compf_04)), 
      use = 'complete.obs') ^ 2
matplot(
  h_compf_69,
  h_compf_04,
  pch = 21:25,
  col = 'black',
  bg = 'grey',
  log = 'xy',
  xlab = '',
  ylab = '2004 handling time',
  axes = FALSE,
  xlim = xylim,
  ylim = xylim
)
title(xlab = '1968-9 handling time', line = 1.4)
abline(0, 1, lty = 2, col = 'grey70')
ats <- 10 ^ seq(-1, 4)
eaxis(1, at = ats)
eaxis(2, at = ats)
box(lwd = 1)
legend(
  'bottomright',
  legend = c(as.expression(bquote(R ^ 2 == .(round(hR2all, 2)))), 
             as.expression(bquote(R[10] ^ 2 == .(round(hlogR2all, 2))))),
  bty = 'n',
  inset = 0.0,
  y.intersp = 1,
  cex = 0.7
)

# feeding rate - feeding and abundance survey sites
xylim <- range(c(f_comp_69, f_comp_04), na.rm = TRUE)
fR2 <-
  cor(as.vector(f_comp_69), as.vector(f_comp_04), 
      use = 'complete.obs') ^ 2
flogR2 <-
  cor(log10(as.vector(f_comp_69)), log10(as.vector(f_comp_04)), 
      use = 'complete.obs') ^ 2
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
  legend = c(as.expression(bquote(R ^ 2 == .(round(fR2, 2)))), 
             as.expression(bquote(R[10] ^ 2 == .(round(flogR2, 2))))),
  bty = 'n',
  inset = 0.0,
  y.intersp = 1,
  cex = 0.7
)

# abundance
xylim <- range(c(N_comp_69, N_comp_04), na.rm = TRUE) + 1E-1
NR2 <-
  cor(as.vector(N_comp_69), as.vector(N_comp_04), 
      use = 'complete.obs') ^ 2
NlogR2 <-
  cor(log10(as.vector(N_comp_69)), log10(as.vector(N_comp_04)), 
      use = 'complete.obs') ^ 2
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
  legend = c(as.expression(bquote(R ^ 2 < .(format(round(NR2, 4), 
                                                    scientific = FALSE)
  ))), as.expression(bquote(R[10] ^ 2 == .(round(NlogR2, 2))))),
  bty = 'n',
  inset = 0.0,
  y.intersp = 1,
  cex = 0.7
)

# attack rate
xylim <- range(c(a_comp_69, a_comp_04), na.rm = TRUE)
aR2 <-
  cor(as.vector(a_comp_69), as.vector(a_comp_04), 
      use = 'complete.obs') ^ 2
alogR2 <-
  cor(log10(as.vector(a_comp_69)), log10(as.vector(a_comp_04)), 
      use = 'complete.obs') ^ 2
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
  legend = c(as.expression(bquote(R ^ 2 == .(round(aR2, 3)
  ))),
  as.expression(bquote(R[10] ^ 2 < .(format(round(NR2, 4), 
                                            scientific = FALSE))
  # as.expression(bquote(R[10] ^ 2 == .(round(alogR2, 2)
  ))),
  bty = 'n',
  inset = 0.0,
  y.intersp = 1,
  cex = 0.7
)

dev.off()

#########################
# Compare relationships between feeding rate and abundance

NfR2 <-
  cor(c(as.vector(N_comp_69[-1, ]), as.vector(N_comp_04[-1, ])), c(as.vector(f_comp_69), as.vector(f_comp_04)), use =
        'complete.obs') ^ 2
NflogR2 <-
  cor(log10(c(
    as.vector(N_comp_69[-1, ]), as.vector(N_comp_04[-1, ])
  )), log10(c(
    as.vector(f_comp_69), as.vector(f_comp_04)
  )), use = 'complete.obs') ^ 2

# Plaxiphora's abundance is 0, so can't log
N_comp_04[N_comp_04 == 0] <- NA

# Seperate models
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


summary(lm(log10(f) ~ log10(N) * Year, data = dat)) # Year & intxn not signif
summary(lm(log10(f) ~ log10(N) + Year, data = dat)) # Year not signif
reg3 <- lm(log10(f) ~ log10(N), data = dat) # Intercept and slope signif
summary(reg3)

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
       xlab = expression(paste('Abundance ', (m^{-2}))),
       ylab = expression(paste('Feeding rate ', (day^{-1}))))
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
  fNlogR2 <- round(summary(reg3)$adj.r.squared,2)
  legend(
    'bottomright',
    legend = as.expression(
      bquote(R[10] ^ 2 == .(round(fNlogR2, 2) ))),
    bty = 'n',
    inset = 0.0,
    y.intersp = 1,
    cex = 0.7
  )
  clip(xlim[1], xlim[2], ylim[1], ylim[2])
  abline(reg3, untf = F)
dev.off()

##############################
# Pred and prey sizes
sizes <- sizes[sizes$Site%in%siteInfo$Site,]

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


regPred <- lm(log(PredSize) ~ Year, data = preysizes)
summary(regPred)
regPrey <- lm(log(PreySize) ~ Year, data = preysizes)
summary(regPrey)
regPredPrey <- lm(log(PredSize) ~ log(PreySize) * Year, data = preysizes)
summary(regPredPrey)
regPredPrey <- lm(log(PredSize) ~ log(PreySize) + Year, data = preysizes)
summary(regPredPrey)


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

ord <- order(names(nTot_69))


# Total count of observations (feeding)
nTotF_69 <- apply(n_compf_69[-1,], 2, sum, na.rm = TRUE)
nTotF_04 <- apply(n_compf_04[-1,], 2, sum, na.rm = TRUE)

# Fraction feeding overall
fF_69 <- binom.confint(nTot_69-n_compf_69['Not Feeding',], nTot_69, 
                       methods = 'wilson')[,-1]
fF_04 <- binom.confint(nTot_04-n_compf_04['Not Feeding',], nTot_04, 
                       methods = 'wilson')[,-1]

# Fraction feeding on H. scobina of those observed feeding
fHs_69 <- binom.confint(n_compf_69['Haustrum scobina',], nTotF_69, 
                       methods = 'wilson')[,-1]
n_compf_04['Haustrum scobina',][is.na(n_compf_04['Haustrum scobina',])] <- 0
fHs_04 <- binom.confint(n_compf_04['Haustrum scobina',], nTotF_04, 
                       methods = 'wilson')[,-1]




size.stats <- sizes %>% group_by(Site, Year) %>%
  summarise(mean = round(mean(PredSize, na.rm = TRUE),1),
            max = max(PredSize, na.rm = TRUE))
size.stats <- as.data.frame(size.stats)
size.stats <- size.stats[size.stats$Site%in%siteInfo$Site,]


s69 <- data.frame(nT = nTot_69[ord],
                fF = 100*round(fF_69[,-c(1,2)],3)[ord,],
                fHs = 100*round(fHs_69[,-c(1,2)],3)[ord,],
                size.stats[size.stats$Year==1969,-c(1,2)])
s04 <- data.frame(nT = nTot_04[ord],
                     fF = 100*round(fF_04[,-c(1,2)],3)[ord,],
                     fHs = 100*round(fHs_04[,-c(1,2)],3)[ord,],
                     size.stats[size.stats$Year==2004,-c(1,2)])

mean_fF_69 <- apply(s69, 2, mean)
mean_fF_04 <- apply(s04, 2, mean)

siteInfo <- siteInfo[order(siteInfo$Lat, decreasing = TRUE),]
siteInfo[,c(2,3)] <- round(siteInfo[,c(2,3)], 4)

summ <- data.frame(
   Site = siteInfo[,1],
   s69$nT,
   s04$nT,
   fF69 = paste0(s69$fF.mean,' (',s69$fF.lower,'-',s69$fF.upper,')'),
   fF04 = paste0(s04$fF.mean,' (',s04$fF.lower,'-',s04$fF.upper,')'),
   fHs69 = paste0(s69$fHs.mean,' (',s69$fHs.lower,'-',s69$fHs.upper,')'),
   fHs04 = paste0(s04$fHs.mean,' (',s04$fHs.lower,'-',s04$fHs.upper,')')
)

summ <- rbind(summ, c('Average',
                      mean_fF_69['nT'], mean_fF_04['nT'],
                      mean_fF_69['fF.mean'], mean_fF_04['fF.mean'],
                      mean_fF_69['fHs.mean'], mean_fF_04['fHs.mean'],
                      mean_fF_69['mean'], mean_fF_04['mean'],
                      mean_fF_69['max'], mean_fF_04['max']))

summ$Site[which(summ$Site == "LeighTabletopRocksandBoulders")] <- "Leigh - Tabletop Rocks"
summ$Site[which(summ$Site == "LeighEchinodermReef")] <- "Leigh - Echinoderm Reef"
summ$Site[which(summ$Site == "Rangitoto Island - Whites Beach")] <- "Rangitoto Island"
summ$Site[which(summ$Site == "Red Beach - Whangaparaoa")] <- "Whangaparaoa"

# summ <- as.matrix(t(summ))
latex(
  summ,
  file='../tables/Paine-SiteSumm.tex',
  cgroup = c('Site', 'Observations', '\\% feeding', 
             '\\% feeding on \\emph{H. scobina}'),
  n.cgroup = c(1, 2, 2, 2),
  colheads = c('','1968-9', '2004','1968-9', '2004','1968-9', '2004'),
  n.rgroup = c(5,1),
  rowname = NULL,
  label = 'tab:summ',
  center = 'centering',
  first.hline.double = FALSE,
  caption="Summary of Paine's 1968-9 and my 2004 feeding observations.  Observerations refers to the total number of whelks inspected. \\% feeding refers to the proportion of observed whelks that were feeding. \\% feeding on \\emph{H. scobina} refers to the proportion of feeding whelks that were feeding on \\emph{Haustrum scobina}.  Parentheticals are the biomial confidence interval (95\\% coverage probability) calculated using the Wilson method.",
)

##########################################################################
##########################################################################
##########################################################################

