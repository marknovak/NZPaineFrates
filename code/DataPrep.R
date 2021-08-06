###################################################################
###################################################################
###################################################################
rm(list = ls())
options(stringsAsFactors = F)

library(tidyr)
library(dplyr)
library(lubridate)
library(rgdal) # for lat-long conversion
library(Hmisc) # for LaTeX table export
  options(xdvicmd='open')

##############################
# Files copied over from 
# https://github.com/NovakLabOSU/NZ_Intertidal_Novak-PhD
fobs <- read.csv("../data/orig/NZ-1969_2004-FeedObs.csv")
abund.raw <- read.csv("../data/orig/NZ-1969_2004-Abunds_orig.csv")
htimes <-
  read.csv("../data/orig/NZ-HandlingTimes-MRegnCoeff-MeasuredOnly.csv",
           skip = 2)
sitecoord <- read.csv("../data/orig/NZ-Access-NZwide-SiteInfo.csv")
leightemps <- read.csv("../data/orig/LeighTemps/data/1-data/SEATEMPW7-monthly_means.csv")[,1:3]

##############################
# Quad-size in database is incorrect!?!
qArea <- 0.3 * 0.3 # quadrat area (m^2) (30 by 30 cm)
abund.raw$Quad_Size <- qArea

#############################
# Abundance surveys
#############################
# Reorganize abundance data
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
abund$Density <- abund$Count / qArea

write.csv(abund, 
          "../data/derived/NZ-1969_2004-Abunds.csv", 
          row.names = FALSE)

colnames(abund)[colnames(abund)=='Species'] <- 'Prey'

#######################
# Feeding surveys
#######################

# Remove anecdotal observations
fobs <- subset(fobs, ObservationType == 'Systematic Census')

############################
# Convert site locations (NZMG) to WGS84
############################
# https://stackoverflow.com/questions/58585146/convert-nzmg-coordinates-to-lat-long

dat <- data.frame(id = sitecoord$Site, 
                  x = sitecoord$NZ.E.coordinate, 
                  y = sitecoord$NZ.N.coordinate)

dat <- dat[dat$id%in%fobs$Site,]
dat <- subset(dat, !is.na(x))

sp::coordinates(dat) = ~x+y

proj4string <- "+proj=nzmg +lat_0=-41 +lon_0=173 +x_0=2510000 +y_0=6023150 +ellps=intl +datum=nzgd49 +units=m +towgs84=59.47,-5.04,187.44,0.47,-0.1,1.024,-4.5993 +nadgrids=nzgd2kgrid0005.gsb +no_defs"

sp::proj4string(dat) = sp::CRS(proj4string) 
data_wgs84 <- sp::spTransform(dat, sp::CRS('+init=epsg:4326'))
siteInfo <- data.frame(Site = data_wgs84$id,
                       Lat = round(coordinates(data_wgs84)[,2], 4),
                       Lon = round(coordinates(data_wgs84)[,1], 4))

#############################
siteInfo$Site <- gsub("([[:lower:]])([[:upper:]])", "\\1 \\2", siteInfo$Site)
fobs$Site <- gsub("([[:lower:]])([[:upper:]])", "\\1 \\2", fobs$Site)

# Align feeding-obs site names with Abund survey names
fobs$Site[which(fobs$Site == "Leigh Waterfalland Penny's Rocks")] <-
  unique(abund$Site[grep('Waterfall', abund$Site)])
fobs$Site[which(fobs$Site == "Leigh Echinoderm Reef")] <-
  'Leigh - Echinoderm Reef'
fobs$Site[which(fobs$Site == "Leigh Tabletop Rocksand Boulders")] <-
  'Leigh - Tabletop Rocks and Boulders'
fobs$Site[which(fobs$Site == "Rangitoto Island White Beach")] <-
  unique(abund$Site[grep('Rangitoto', abund$Site)])
fobs$Site[which(fobs$Site == "Red Beach")] <-
  unique(abund$Site[grep('Red Beach - Whangaparaoa', abund$Site)])
fobs$Site[which(fobs$Site == "Whangarei Baptist Camp")] <-
  'Whangarei'

siteInfo$Site[which(siteInfo$Site == "Leigh Waterfalland Penny's Rocks")] <-
  unique(abund$Site[grep('Waterfall', abund$Site)])
siteInfo$Site[which(siteInfo$Site == "Leigh Echinoderm Reef")] <-
  'Leigh - Echinoderm Reef'
siteInfo$Site[which(siteInfo$Site == "Leigh Tabletop Rocksand Boulders")] <-
  'Leigh - Tabletop Rocks and Boulders'
siteInfo$Site[which(siteInfo$Site == "Rangitoto Island White Beach")] <-
  unique(abund$Site[grep('Rangitoto', abund$Site)])
siteInfo$Site[which(siteInfo$Site == "Red Beach")] <-
  unique(abund$Site[grep('Red Beach - Whangaparaoa', abund$Site)])
siteInfo$Site[which(siteInfo$Site == "Whangarei Baptist Camp")] <-
  'Whangarei'

#############################
# Which sites can be compared
tab.fobs <- table(unique(fobs[, c('Site', 'Year')]))
  tab.fobs
tab.abund <- table(unique(abund[, c('Site', 'Year')]))
  tab.abund
sites.fcomp <- names(which(apply(tab.fobs, 1, sum) > 1 
                           & tab.fobs[, 3] == 1))

# Pick sites that have both feeding obs and abundance surveys
# and select abundance transects to compare
x <- unique(data.frame(abund$Site, abund$Zone, abund$Year))
x <- x[order(x[, 3]), ]
x

abund <- subset(
            abund,
            Site == 'Leigh - Waterfall Rocks' 
              & grepl('Transect#2', abund$Zone) |
            Site == 'Leigh - Waterfall Rocks' 
              & Year == 2004 |
            Site == 'Rangitoto Island - Whites Beach' 
              & Zone == 'Oyster Zone' |
            Site == 'Rangitoto Island - Whites Beach' 
              & Year == 2004 |
            Site == 'Red Beach - Whangaparaoa' 
              & grepl('Oyster', abund$Zone)
            )

abund <-
  abund %>% group_by(Site, Year, Prey) %>%
  summarise(N.mean = mean(Density, na.rm = TRUE),
            N.sd = sd(Density, na.rm = TRUE),
            N.n = n())
abund$N.mean[is.na(abund$N.mean)] <- 0

################################
# Feeding observations
################################
# Add Year-Month mean temperatures
fobs$Month <- month(mdy(fobs$Date))
colnames(leightemps)[1] <- 'Temp'
fobs <- merge(fobs, leightemps,
              by.x = c('Year','Month'),
              by.y = c('year','month'),
              all.x = TRUE)

fobs$Year[which(fobs$Year == 1968)] <- 1969 # For convenience
Sites <- sort(unique(fobs$Site))

sort(unique(fobs$Counter))

# Assign all observations not by me to Paine
fobs$Counter[!grepl("Novak", fobs$Counter)] <- 'Paine'
table(fobs$Site, fobs$Counter)

# Match prey to the lab-based detection (handling) time coefficients
htimes <- subset(htimes, ConLevel == 0.1 & Type != 'Unweighted')
prey <- unique(fobs$Prey)

meas <- unique(data.frame(Pred = htimes$Pred, Prey = htimes$Prey))
tomatch <- merge(data.frame(Prey = prey), meas, all.x = TRUE)[, c(2, 1)]
tomatch <- data.frame(tomatch, MatchToPred = NA, MatchToPrey = NA)
write.csv(tomatch,
          '../data/derived/NZ-HandlingTimes-Paine-Prey2match.csv',
          row.names = FALSE)

Match <- read.csv('../data/orig/NZ-HandlingTimes-Paine-PreyMatch.csv')

rem <- which(Match$Pred == Match$MatchToPred 
             & Match$Prey == Match$MatchToPrey)
latex(
  Match[-rem,],
  file='../tables/Paine-PredPreyMatches.tex',
  cgroup = c('Unmeasured','Matched to measured'),
  n.cgroup = c(2, 2),
  colheads = c('Predator','Prey','Predator','Prey'),
  rowname = NULL,
  label = 'tab:matches',
  first.hline.double = FALSE,
  where = "!htbp",
  caption="Prey for which detection times had not been measured in the laboratory experiments of \\citep{Novak:2013qg} were assigned the regression coefficients of prey species for which they had been measured."
)

# Observations with unmeasured pred or prey sizes get NA
fobs$PredSize[fobs$PredSize <= 0] <- NA
fobs$PreySize[fobs$PreySize <= 0] <- NA

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

# Divide by 24 hours to convert to days
fobs$htime <- fobs$htime/24

# Look for problems
subset(fobs, is.na(htime) & Prey != 'Not Feeding')

# Replace problems with species means
mHtime <-
  aggregate(list(htime = fobs$htime),
            by = list(Prey = fobs$Prey),
            mean,
            na.rm = TRUE)

prey <- unique(fobs$Prey)
for (i in 1:length(prey)) {
  mn <- mean(fobs$htime[fobs$Prey == prey[i]], 
             na.rm = TRUE)
  fobs$htime[which(fobs$Prey == prey[i] &
                     is.na(fobs$htime) == TRUE)] <- mn
}

###########################

nh <- 
  fobs %>% group_by(Site, Prey, Year) %>%
  summarise(n.obs = n(),
            h.mean = mean(htime) )

N <- data.frame(abund[, c('Site','Year','Prey','N.mean')])

# Calculate 
# proportion feeding (p_i = n_i / sum_n)
#         (for feeding rate calculation)
# and
# ratio of feeding to not feeding
#         (for attack rate calculation)
p <- 
  nh %>% group_by(Site, Year) %>%
  mutate(pi = n.obs / sum(n.obs),
         pi0 = n.obs / n.obs[Prey == "Not Feeding"])

# Merge all into single dataframe
all.summ <- Reduce(function(x, y) merge(x, y, all=TRUE), list(nh, p, N))
all.summ$N.mean[all.summ$N.mean==0] <- NA

############################
# Predator sizes
sizes <- fobs[,c('Pred','Year','Site','Area','PredSize','Prey','PreySize')]

############################
# Site information
fobs.SiteYr <- as.matrix(table(fobs$Site, fobs$Year))
fobs.SiteYr[fobs.SiteYr > 0] <- 'x'
fobs.SiteYr[fobs.SiteYr == 0] <- ''
fobs.SiteYr <- data.frame(Site = rownames(fobs.SiteYr),
                     'F1968.9' = fobs.SiteYr[,1],
                     'F2004' = fobs.SiteYr[,2]  )

abund.SiteYr <- as.matrix(table(abund$Site, abund$Year))
abund.SiteYr[abund.SiteYr > 0] <- 'x'
abund.SiteYr[abund.SiteYr == 0] <- ''
abund.SiteYr <- data.frame(Site = rownames(abund.SiteYr),
                     'A1968-9' = abund.SiteYr[,1],
                     'A2004' = abund.SiteYr[,2]  )

SiteYr <- merge(fobs.SiteYr, abund.SiteYr, all = TRUE)
SiteYr[is.na(SiteYr)] <- ''

siteInfo <- merge(siteInfo, SiteYr, all.y = TRUE)
# siteInfo[is.na(siteInfo)] <- ''
siteInfo <- siteInfo[order(-xtfrm(siteInfo$F1968.9),-siteInfo$Lat),]

latex(
  siteInfo,
  file='../tables/Paine-Sites.tex',
  cgroup = c('', '', '','Feeding', 'Abundance'),
  n.cgroup = c(1, 1, 1, 2, 2),
  colheads = c('Site','Latitude','Longitude',
               '1968-9', '2004','1968-9', '2004'),
  rowname = NULL,
  label = 'tab:sites',
  first.hline.double = FALSE,
  where = "!htbp",
  caption="The locations where Paine and I surveyed \\emph{Haustrum haustorium}'s diet and the abundances of its prey in 1968-9 and 2004 for which data are posted to the public repositories indicated in the main text.  Missing coordinates are unknown."
)

############################
write.csv(all.summ,
          '../data/derived/NZ-1969_2004-tab_Summarized.csv',
          row.names = FALSE)

write.csv(sizes,
          '../data/derived/NZ-1969_2004-tab_Sizes.csv',
          row.names = FALSE)

write.csv(siteInfo,
          '../data/derived/NZ-1969_2004-Site-LatLon.csv',
          row.names = FALSE)

write.csv(sites.fcomp,
          '../data/derived/NZ-1969_2004-Site-focal.csv',
          row.names = FALSE)


###################################################################
###################################################################
###################################################################