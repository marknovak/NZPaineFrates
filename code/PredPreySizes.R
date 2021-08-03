#####################
# Pred and prey sizes
#####################
#########################################################################
#########################################################################
#########################################################################
rm(list = ls())
options(stringsAsFactors = F)

library(sfsmisc) # for eaxis
library(Hmisc) # for LaTeX table export
options(xdvicmd='open')
library(xtable) # for LaTeX summary table export
############################
focalSites <- read.csv("../data/derived/NZ-1969_2004-Site-focal.csv")$x
sizes <- read.csv("../data/derived/NZ-1969_2004-tab_Sizes.csv")

############################
sizes <- sizes[sizes$Site%in%focalSites,]
sizes <- sizes[!is.na(sizes$PredSize),]

preysizes <- sizes[sizes$PreySize > 0 & !is.na(sizes$PreySize) ,]
breaks = seq(0, max(preysizes$PreySize) + 1, 1)
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
            preysizes$PreySize[preysizes$Year==2004]), 
       function(x){c(mean = mean(x),
                     sd = sd(x))})

ks.pred <- ks.test(sizes$PredSize[sizes$Year==1969],
                   sizes$PredSize[sizes$Year==2004])
lapply(list(sizes$PredSize[sizes$Year==1969],
            sizes$PredSize[sizes$Year==2004]), 
       function(x){c(mean = mean(x),
                     sd = sd(x))})

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
  digits = 3,
  where = "!htbp",
  caption="Summary table for the regression of predator 
  size on prey size,
  time period (\\textit{Year}), 
  and their interaction."
)
latex(
  xtable(summary(regPredPreyYear)),
  file='../tables/Paine-PredPreySizeYear.tex',
  # rowname = NULL,
  rowlabel = '',
  label = 'tab:SizeY',
  # center = 'centering',
  first.hline.double = FALSE,
  where = "!htbp",
  digits = 3,
  caption="Summary table for the regression of predator size 
  on prey size and 
  time period (\\textit{Year})."
)
latex(
  xtable(summary(regPredPrey)),
  file='../tables/Paine-PredPreySize.tex',
  # rowname = NULL,
  rowlabel = '',
  label = 'tab:Size',
  # center = 'centering',
  first.hline.double = FALSE,
  digits = 3,
  where = "!htbp",
  caption="Summary table for the regression of predator size 
  on prey size."
)


regPredPrey69 <- lm(log(PredSize) ~ log(PreySize), 
                    data = preysizes[preysizes$Year==1969,])
summary(regPredPrey69)
rng69 <- range(preysizes[preysizes$Year==1969,]$PreySize)

regPredPrey04 <- lm(log(PredSize) ~ log(PreySize), 
                    data = preysizes[preysizes$Year==2004,])
summary(regPredPrey04)
rng04 <- range(preysizes[preysizes$Year==2004,]$PreySize)


# To visualize non-feeding predators by year
sizes$PreySize[sizes$Prey=="Not Feeding" & sizes$Year==2004] <- 0.25
sizes$PreySize[sizes$Prey=="Not Feeding" & sizes$Year==1969] <- -0.25


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
     xlab = expression(paste('Prey size (',italic(mm),')')),
     ylab = expression(paste('Predator size (',italic(mm),')')),
     xlim = xlims, 
     ylim = ylims
)
abline(0,1, col = 'grey70', lty = 2)
lat = rng69[1]:rng69[2]
lines(lat,
      exp(predict(regPredPrey69, data.frame('PreySize'=lat))),
      col = 'grey20',
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
       pt.bg = c('grey20','grey70'),
       inset = 0.01)
box(lwd = 1)


par(mar = c(0,3,1,0))
barplot(Prey$counts, axes = FALSE, ylim = c(0, top[1]), space = 0,
        col = 'grey20')
axis(2, cex.axis = 0.5)
barplot(Prey04$counts, axes = FALSE, ylim = c(0, top[1]), space = 0,
        col = 'grey70', add = TRUE)

par(mar = c(3,0,0,1))
barplot(Pred$counts, axes = FALSE, xlim = c(0, top[2]), space = 0, 
        horiz = TRUE,
        col = 'grey20')
axis(1, cex.axis = 0.5, mgp = c(0, 0, 0))
barplot(Pred04$counts, axes = FALSE, xlim = c(0, top[2]), space = 0, 
        horiz = TRUE,
        col = 'grey70', add = TRUE)

dev.off()

#######################################################################
#######################################################################
#######################################################################