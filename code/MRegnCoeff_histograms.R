# Summary of detection-time sensitivities (regression coefficients)
# from Novak 2013 as requested by reviewer

htimes <-
  read.csv("../data/orig/NZ-HandlingTimes-MRegnCoeff-MeasuredOnly.csv",
           skip = 2)

# As done in 'DataPrep.R'
htimes <- subset(htimes, ConLevel == 0.1 & Type != 'Unweighted')

# Here also remove 'mean' estimates (i.e. estimates based on only 1 observation)
# but keep H.scobina since they are used too
htimes <- subset(htimes, Type != 'Mean')

# Calculate mean effect size
summ <- round(apply(cbind(htimes$logPredSizeC, 
                    htimes$logPreySizeC, 
                    htimes$logTempC), 
              2, 
              function(x){c(mean = mean(x), 
                            sd = sd(x))}),
              2)

pdf('../figs/MRegnCoeff_histograms.pdf', height = 5, width = 3)
par(
  cex = 0.8,
  cex.axis = 0.9,
  cex.lab = 1,
  tcl = -0.2,
  mar = c(3, 4, 1, 1),
  mgp = c(1.5, 0.3, 0),
  las = 1,
  yaxs = 'i',
  mfrow = c(3, 1)
)
  bks <- 10
  
  hist(htimes$logPredSizeC, 
       breaks = bks, 
       main = '',
       xlab = 'Effect size')
    abline(v = summ['mean', 1], 
           col = 'grey50', 
           lwd = 3)
    box(lwd = 1,
        bty = 'l')
    mtext('(a)', 
          3, 
          adj = -0.2, 
          line = 0, 
          cex = 0.8)
    mtext(bquote(bar(italic(x)) == .(summ['mean', 1])), 
          3,
          adj = 1,
          line = -1,
          cex = 0.7)
    mtext(bquote(italic(sd) == .(summ['sd', 1])), 
          3,
          adj = 1,
          line = -2, 
          cex = 0.7)
    
  hist(htimes$logPreySizeC, 
       breaks = bks, 
       main = '',
       xlab = 'Effect size')
    abline(v = summ['mean', 2], 
         col = 'grey50', 
         lwd = 3)
    box(lwd = 1,
        bty = 'l')
    mtext('(b)', 
          3, 
          adj = -0.2, 
          line = 0, 
          cex = 0.8)
    mtext(bquote(bar(italic(x)) == .(summ['mean', 2])), 
          3,
          adj = 1,
          line = -1,
          cex = 0.7)
    mtext(bquote(italic(sd) == .(summ['sd', 2])), 
          3,
          adj = 1,
          line = -2, 
          cex = 0.7)
    
  hist(htimes$logTempC, 
       breaks = bks, 
       main = '',
       xlab = 'Effect size')
    abline(v = summ['mean', 3], 
         col = 'grey50', 
         lwd = 3)
    box(lwd = 1,
        bty = 'l')
    mtext('(c)', 
          3, 
          adj = -0.2, 
          line = 0, 
          cex = 0.8)
    mtext(bquote(bar(italic(x)) == .(summ['mean', 3])), 
          3,
          adj = 1,
          line = -1,
          cex = 0.7)
    mtext(bquote(italic(sd) == .(summ['sd', 3])), 
          3,
          adj = 1,
          line = -2, 
          cex = 0.7)
dev.off()

