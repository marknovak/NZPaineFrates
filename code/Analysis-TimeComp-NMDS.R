library(vegan)
library(lubridate)

abund.raw <- read.csv("../data/orig/NZ-1969_2004-Abunds_orig.csv")

abund <- subset(
  abund.raw,
  Site == 'Leigh - Waterfall Rocks' & grepl('Transect#2', Zone) |
  Site == 'Leigh - Waterfall Rocks' & year(dmy(Date)) == 2004 |
  Site == 'Rangitoto Island - Whites Beach' & Zone == 'Oyster Zone' |
  Site == 'Rangitoto Island - Whites Beach' & year(dmy(Date)) == 2004 |
  Site == 'Red Beach - Whangaparaoa' & grepl('Oyster', Zone)
)


dat <- abund %>% select(-one_of(c('Date','Site','Zone','Quad_Size','Quad'))) 
dat <- dat[, -which(apply(dat, 2, sum, na.rm = TRUE)==0)]
dat[is.na(dat)] <- 0
mat <- as.matrix(dat)

row.names(mat) <- paste0(abbreviate(abund$Site),'-',year(dmy(abund$Date)))

              
# Non-metric Multi-dimensional Scaling
nmds <- metaMDS(mat, 
                distance = "bray",
                maxit = 999, 
                trymax = 500)

site.pch <- c(23,22,21)
year.bg <- c('grey30','grey90')
grps <- table(rownames(mat))
grp.cols <- c('grey20','grey80')

pdf('../figs/Paine-NMDS.pdf',
    width = 8,
    height = 4)
par(
  mfrow = c(1,2),
  pty = 's',
  cex = 0.8,
  cex.axis = 0.9,
  cex.lab = 1,
  tcl = -0.2,
  mar = c(3, 4, 1, 1),
  mgp = c(2, 0.3, 0),
  las = 1
)
  stressplot(nmds)
  
  plot(nmds, type = "n")
  legend('bottomleft',
         legend = bquote(Stress == .(round(nmds$stress,3))),
         cex = 0.9,
         bty = 'n')
  legend('bottomright',
         legend = c('Leigh (WR)','Rangitoto','Whangaporoa'),
         pch = site.pch,
         pt.bg = NA,
         cex = 0.9)
  legend('topright',
         legend = levels(factor(year(dmy(abund$Date)))),
         pch = 22,
         pt.bg = year.bg,
         cex = 0.9)
  ordihull(nmds, 
           groups = rownames(mat), 
           draw = "polygon", 
           lty = 1,
           col = grp.cols,
           alpha = 0.2)
  points(nmds, 
           display = "sites", 
           pch = site.pch[as.numeric(factor(abund$Site))], 
           bg = year.bg[as.numeric(factor(year(dmy(abund$Date))))],
           cex = 1)
dev.off()
####################################################################
####################################################################
####################################################################