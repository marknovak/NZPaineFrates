###############################################################################
###############################################################################
###############################################################################
rm(list = ls())
options(stringsAsFactors = F)

library(gdata) # for xls import
library(tidyr)
library(dplyr)
library(lubridate)

##############################
##############################
qArea <- 0.33 # quadrat area (m^2)
assu.Temp <- 14 # assumed temperature

##############################
fobs <- read.csv("../data/orig/NZ-1969_2004-FeedObs.csv")
abund.raw <- read.xls("../data/orig/NZ-1969_2004-Abunds_orig.xls")
htimes <-
  read.csv("../data/orig/NZ-HandlingTimes-MRegnCoeff-MeasuredOnly.csv",
           skip = 2) # File copied over from NZ Rscripts folder

fobs <- subset(fobs, ObservationType == 'Systematic Census')
###########################

abund.raw <-
  abund.raw[, -which(apply(abund.raw, 2, 
                           function(x) { all(is.na(x)) }))]
abund <-
  gather(
    abund.raw,
    Species,
    Count,
    6:ncol(abund.raw),
    -Date,
    -Site,
    -Zone,
    -Quad_Size,
    -Quad,
    na.rm = FALSE
  )
abund$Year <- year(dmy(abund$Date))
abund$Species <- gsub('([[:punct:]])|\\s+', ' ', abund$Species)
write.csv(abund, "../data/derived/NZ-1969_2004-Abunds.csv", row.names = FALSE)
# ---------------
abund$Dens <- abund$Count / qArea
#####################################
abund.stats <-
  abund %>% group_by(Site, Zone, Year, Species) %>% 
    summarise(N = mean(Dens, na.rm = TRUE))
write.csv(abund.stats,
          '../data/derived/NZ-1969_2004-AbundZoneMeans.csv',
          row.names = FALSE)

###################################
# Align abundance and feeding obs site names
fobs$Site[which(fobs$Site == "LeighWaterfallandPenny'sRocks")] <-
  unique(abund$Site[grep('Waterfall', abund$Site)])
fobs$Site[which(fobs$Site == "RangitotoIslandWhiteBeach")] <-
  unique(abund$Site[grep('Rangitoto', abund$Site)])
fobs$Site[which(fobs$Site == "RedBeach")] <-
  unique(abund$Site[grep('Red Beach - Whangaparaoa', abund$Site)])
fobs$Site[which(fobs$Site == "WhangareiBaptistCamp")] <-
  unique(abund$Site[grep('Whangarei', abund$Site)])

tab.fobs <- table(unique(fobs[, c('Site', 'Year')]))
tab.abund <- table(unique(abund[, c('Site', 'Year')]))
tab.fobs
tab.abund

sites.fcomp <- names(which(apply(tab.fobs, 1, sum) > 1 & tab.fobs[, 3] == 1))

# Pick sites that have both feeding obs and abundance surveys
# and select abundance transects to compare
x <-
  unique(data.frame(abund.stats$Site, abund.stats$Zone, abund.stats$Year))
x <- x[order(x[, 3]), ]
x

y <- abund
Y <-
  subset(
    y,
    Site == 'Leigh - Waterfall Rocks' & grepl('Transect#2', y$Zone) |
      Site == 'Leigh - Waterfall Rocks' & Year == 2004 |
      Site == 'Rangitoto Island - Whites Beach' &
      Zone == 'Oyster Zone' |
      Site == 'Rangitoto Island - Whites Beach' &
      Year == 2004 |
      Site == 'Red Beach - Whangaparaoa' &
      grepl('Oyster', y$Zone)
  )

abund_comp <-
  Y %>% group_by(Site, Year, Species) %>% summarise(N = mean(Dens, na.rm = TRUE))

################################
fobs$Temp <-
  assu.Temp  # Should probably break up by site (latitude and island side)
fobs$Year[which(fobs$Year == 1968)] <- 1969 # For convenience

fobs$Counter[!grepl("Novak", fobs$Counter)] <- 'Paine'
table(fobs$Site, fobs$Counter)
# nPaine<-length(grep("Novak",fobs$Counter,invert=TRUE))
# nMe<-length(grep("Novak",fobs$Counter))
# nPaineF<-length(grep("Novak",fobs$Counter[which(fobs$Prey!='Not Feeding')],invert=TRUE))
# nMeF<-length(grep("Novak",fobs$Counter[which(fobs$Prey!='Not Feeding')]))
# c(Paine=nPaine,Me=nMe,Total=nrow(fobs))
# c(PaineFracFeed=nPaineF/nPaine,MeFracFeed=nMeF/nMe,Total=nrow(subset(fobs,Prey!='Not Feeding'))/nrow(fobs))

prey <- sort(unique(fobs$Prey[fobs$Prey != "Not Feeding"]))

htimes <- subset(htimes, ConLevel == 0.1 & Type != 'Unweighted')

meas <- unique(data.frame(Pred = htimes$Pred, Prey = htimes$Prey))
match <- merge(data.frame(Prey = prey), meas, all.x = TRUE)[, c(2, 1)]
match <- data.frame(match, MatchToPred = NA, MatchToPrey = NA)
write.csv(match,
          '../data/derived/NZ-HandlingTimes-Paine-Prey2match.csv',
          row.names = FALSE)
Match <- read.csv('../data/orig/NZ-HandlingTimes-Paine-PreyMatch.csv')

# drop unidentified prey
fobs <- subset(fobs, Prey != 'UNID')
# drop observations with unmeasured pred sizes
fobs <- subset(fobs, PredSize > 0)

fobs <- merge(fobs, Match, all.x = TRUE)
fobs <-
  merge(
    fobs,
    htimes,
    all.x = TRUE,
    by.x = c('MatchToPred', 'MatchToPrey'),
    by.y = c('Pred', 'Prey')
  )
fobs$htime <-
  exp(
    fobs$logIntC 
    + fobs$logPredSizeC * log(fobs$PredSize) 
    + fobs$logPreySizeC * log(fobs$PreySize) 
    + fobs$logTempC * log(fobs$Temp)
  )

# Look for htime problems
aggregate(list(htime = fobs$htime), by = list(Prey = fobs$Prey), mean)
subset(fobs[, c('Site',
                'CensusID',
                'Pred',
                'Prey',
                'PredSize',
                'PreySize',
                'htime')], htime > 100)
subset(fobs, is.na(htime) & Prey != 'Not Feeding')
fobs$htime[is.infinite(fobs$htime) == TRUE] <- NA
mHtime <-
  aggregate(list(htime = fobs$htime),
            by = list(Prey = fobs$Prey),
            mean,
            na.rm = TRUE)

# Replace problems with species means
fobs$htime[fobs$htime == 0] <- NA
for (i in 1:length(prey)) {
  mn <- mean(fobs$htime[fobs$Prey == prey[i]], na.rm = TRUE)
  fobs$htime[which(fobs$Prey == prey[i] &
                     is.na(fobs$htime) == TRUE)] <- mn
}

###########################
h <-
  aggregate(list(htime = fobs$htime),
            by = list(
              Site = fobs$Site,
              Year = fobs$Year,
              Prey = fobs$Prey
            ),
            mean)
n <-
  aggregate(list(n = fobs$Prey),
            by = list(
              Site = fobs$Site,
              Year = fobs$Year,
              Prey = fobs$Prey
            ),
            length)

nSites <- sort(unique(fobs$Site))

abund.stats$SiteZone <- paste(abund.stats$Site, '-', abund.stats$Zone)
abund.stats <-
  data.frame(abund.stats)
abund.stats <- abund.stats[, -c(1, 2)]

############################
h_all <- spread(h, Site, htime)
n_all <- spread(n, Site, n)
N_all <- spread(abund.stats, SiteZone, N)

###########################################
# Revisited sites with just feeding surveys
h_compf <- subset(h, Site %in% sites.fcomp)
h_compf <- spread(h_compf, Prey, htime)
h_compf <- gather(h_compf, Prey, htime, -Site, -Year)
h_compf <- spread(h_compf, Site, htime)

n_compf <- subset(n, Site %in% sites.fcomp)
n_compf <- spread(n_compf, Prey, n)
n_compf <- gather(n_compf, Prey, n, -Site, -Year)
n_compf <- spread(n_compf, Site, n)

#########################################################
# Revisited sites with both abundance and feeding surveys
h_comp <-
  subset(
    h,
    Site == 'Leigh - Waterfall Rocks' |
      Site == 'Rangitoto Island - Whites Beach' |
      Site == 'Red Beach - Whangaparaoa'
  )
h_comp <- spread(h_comp, Prey, htime)
h_comp <- gather(h_comp, Prey, htime, -Site, -Year)
h_comp <- spread(h_comp, Site, htime)

n_comp <-
  subset(
    n,
    Site == 'Leigh - Waterfall Rocks' |
      Site == 'Rangitoto Island - Whites Beach' |
      Site == 'Red Beach - Whangaparaoa'
  )
n_comp <- spread(n_comp, Prey, n)
n_comp <- gather(n_comp, Prey, n, -Site, -Year)
n_comp <- spread(n_comp, Site, n)

N_comp <-
  merge(
    data.frame(Prey = unique(n_comp$Prey)),
    data.frame(abund_comp),
    all.x = TRUE,
    by.x = 'Prey',
    by.y = 'Species'
  )
N_comp$Site[is.na(N_comp$Site) == TRUE] <- 'Leigh - Waterfall Rocks'
N_comp$Year[is.na(N_comp$Year) == TRUE] <- '2004'
N_comp <- spread(N_comp, Year, N)
N_comp <- gather(N_comp, Year, N, -Site, -Prey)
N_comp <- spread(N_comp, Prey, N)
N_comp <- gather(N_comp, Prey, N, -Site, -Year)
N_comp <- spread(N_comp, Site, N)

###############################################
h_all_69 <- subset(h_all, Year == 1969)[, -1]
h_all_04 <- subset(h_all, Year == 2004)[, -1]
n_all_69 <- subset(n_all, Year == 1969)[, -1]
n_all_04 <- subset(n_all, Year == 2004)[, -1]
N_all_69 <- subset(N_all, Year == 1969)[, -1]
N_all_04 <- subset(N_all, Year == 2004)[, -1]

# Reorder sites to align with feeding and abundance sites
h_compf <- h_compf[, c(1:3, 6, 7, 4, 5)]
n_compf <- n_compf[, c(1:3, 6, 7, 4, 5)]

h_compf_69 <- subset(h_compf, Year == 1969)[, -1]
h_compf_04 <- subset(h_compf, Year == 2004)[, -1]
n_compf_69 <- subset(n_compf, Year == 1969)[, -1]
n_compf_04 <- subset(n_compf, Year == 2004)[, -1]

h_comp_69 <- subset(h_comp, Year == 1969)[, -1]
h_comp_04 <- subset(h_comp, Year == 2004)[, -1]
n_comp_69 <- subset(n_comp, Year == 1969)[, -1]
n_comp_04 <- subset(n_comp, Year == 2004)[, -1]
N_comp_69 <- subset(N_comp, Year == 1969)[, -1]
N_comp_04 <- subset(N_comp, Year == 2004)[, -1]

nf <- grep('Not Feeding', h_compf_69$Prey)
h_compf_69 <- h_compf_69[c(nf, seq(1, nrow(h_compf_69))[-nf]), ]
nf <- grep('Not Feeding', n_compf_69$Prey)
n_compf_69 <- n_compf_69[c(nf, seq(1, nrow(n_compf_69))[-nf]), ]

nf <- grep('Not Feeding', h_compf_04$Prey)
h_compf_04 <- h_compf_04[c(nf, seq(1, nrow(h_compf_04))[-nf]), ]
nf <- grep('Not Feeding', n_compf_04$Prey)
n_compf_04 <- n_compf_04[c(nf, seq(1, nrow(n_compf_04))[-nf]), ]

nf <- grep('Not Feeding', h_comp_69$Prey)
h_comp_69 <- h_comp_69[c(nf, seq(1, nrow(h_comp_69))[-nf]), ]
nf <- grep('Not Feeding', n_comp_69$Prey)
n_comp_69 <- n_comp_69[c(nf, seq(1, nrow(n_comp_69))[-nf]), ]
nf <- grep('Not Feeding', N_comp_69$Prey)
N_comp_69 <- N_comp_69[c(nf, seq(1, nrow(N_comp_69))[-nf]), ]

nf <- grep('Not Feeding', h_comp_04$Prey)
h_comp_04 <- h_comp_04[c(nf, seq(1, nrow(h_comp_04))[-nf]), ]
nf <- grep('Not Feeding', n_comp_04$Prey)
n_comp_04 <- n_comp_04[c(nf, seq(1, nrow(n_comp_04))[-nf]), ]
nf <- grep('Not Feeding', N_comp_04$Prey)
N_comp_04 <- N_comp_04[c(nf, seq(1, nrow(N_comp_04))[-nf]), ]

############################
# Predator sizes
sizes <- fobs[,c('Pred','Year','Site','Area','PredSize','Prey','PreySize')]

############################
write.csv(h_all_69,
          "../data/derived/NZ-1969_2004-tab_Htime_all_69.csv",
          row.names = FALSE)
write.csv(n_all_69,
          "../data/derived/NZ-1969_2004-tab_Fobs_all_69.csv",
          row.names = FALSE)
write.csv(N_all_69,
          '../data/derived/NZ-1969_2004-tab_Abund_all_69.csv',
          row.names = FALSE)

write.csv(h_all_04,
          "../data/derived/NZ-1969_2004-tab_Htime_all_04.csv",
          row.names = FALSE)
write.csv(n_all_04,
          "../data/derived/NZ-1969_2004-tab_Fobs_all_04.csv",
          row.names = FALSE)
write.csv(N_all_04,
          '../data/derived/NZ-1969_2004-tab_Abund_all_04.csv',
          row.names = FALSE)

write.csv(h_compf_69,
          "../data/derived/NZ-1969_2004-tab_Htime_compf_69.csv",
          row.names = FALSE)
write.csv(n_compf_69,
          "../data/derived/NZ-1969_2004-tab_Fobs_compf_69.csv",
          row.names = FALSE)

write.csv(h_compf_04,
          "../data/derived/NZ-1969_2004-tab_Htime_compf_04.csv",
          row.names = FALSE)
write.csv(n_compf_04,
          "../data/derived/NZ-1969_2004-tab_Fobs_compf_04.csv",
          row.names = FALSE)

write.csv(h_comp_69,
          "../data/derived/NZ-1969_2004-tab_Htime_comp_69.csv",
          row.names = FALSE)
write.csv(n_comp_69,
          "../data/derived/NZ-1969_2004-tab_Fobs_comp_69.csv",
          row.names = FALSE)
write.csv(N_comp_69,
          '../data/derived/NZ-1969_2004-tab_Abund_comp_69.csv',
          row.names = FALSE)

write.csv(h_comp_04,
          "../data/derived/NZ-1969_2004-tab_Htime_comp_04.csv",
          row.names = FALSE)
write.csv(n_comp_04,
          "../data/derived/NZ-1969_2004-tab_Fobs_comp_04.csv",
          row.names = FALSE)
write.csv(N_comp_04,
          '../data/derived/NZ-1969_2004-tab_Abund_comp_04.csv',
          row.names = FALSE)

write.csv(sizes,
          '../data/derived/NZ-1969_2004-tab_Sizes.csv',
          row.names = FALSE)

#################################################################################
#################################################################################
#################################################################################

