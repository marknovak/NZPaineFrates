#########################################################################
#########################################################################
#########################################################################
###############
# Summary table
###############
rm(list = ls())
options(stringsAsFactors = F)

library(dplyr)
library(binom) # for confidence intervals
library(Hmisc) # for LaTeX table export
options(xdvicmd='open')
library(xtable) # for LaTeX summary table export

###########################
dat <- read.csv('../data/derived/NZ-1969_2004-tab_Summarized.csv')
focalSites <- read.csv("../data/derived/NZ-1969_2004-Site-focal.csv")$x
sizes <- read.csv("../data/derived/NZ-1969_2004-tab_Sizes.csv")
###########################

# Total count of observations (feeding + not feeding)
dat$n.obs[is.na(dat$n.obs) & dat$Site == 'Red Beach - Whangaparaoa'
          & dat$Year == 2004 & dat$Prey =='Haustrum scobina'] <- 0

dat <- dat %>% 
  filter(Site%in%focalSites)

obs <- dat %>%
  group_by(Prey, Year) %>%
  summarise(tot.n.obs = sum(n.obs, na.rm = TRUE)) %>%
  pivot_wider(id_cols = Prey,
                       names_from = Year,
                       values_from = tot.n.obs) %>%
  data.frame() %>%
  replace_na(list(X1969=0, X2004=0))
  
obs <- obs[order(apply(obs[,-1], 1, sum), 
                 decreasing = TRUE),]
obs <- obs[-which(obs$X1969 == 0 & obs$X2004 == 0),]
apply(obs[,-1], 2, sum)

# Prey observed by Paine but not me
obs[which(obs$X1969 > 0 & obs$X2004 == 0),]
# Prey observed by me but not Paine
obs[which(obs$X1969 == 0 & obs$X2004 > 0),]
sum(obs[which(obs$X1969 == 0 & obs$X2004 > 0),]$X2004)

summ <- dat %>%
  group_by(Site, Year) %>%
  summarise(tot.obs = sum(n.obs, na.rm = TRUE),
            feed.obs = sum(n.obs[Prey != "Not Feeding"], na.rm = TRUE),
            Hs.obs = sum(n.obs[Prey == 'Haustrum scobina'], na.rm = TRUE)) %>%
  mutate(fF = binom.confint(feed.obs, tot.obs, # percent feeding
                            method = 'wilson')[,-c(1:3)],
         fF.Hs = binom.confint(Hs.obs, feed.obs, # percent of feeding that fed on Hs
                               method = 'wilson')[,-c(1:3)],
         fF = round(100*fF, 1),
         fF.Hs = round(100*fF.Hs, 1)) %>%
  select(-feed.obs, -Hs.obs)


summ.w1 <- pivot_wider(summ, 
                      id_cols = 'Site',
                      names_from = 'Year',
                      values_from = c('tot.obs','fF','fF.Hs')) 

cor.test(summ.w1$tot.obs_1969, summ.w1$tot.obs_2004)
cor.test(log10(summ.w1$tot.obs_1969), log10(summ.w1$tot.obs_2004))
cor.test(summ.w1$tot.obs_1969, summ.w1$tot.obs_2004, method = 'spearman')

summ.tot <- c('Sum/Average',
              c(sum(summ.w1$tot.obs_1969), sum(summ.w1$tot.obs_2004)),
              format(
                round(
                  c(mean(summ.w1$fF_1969$mean), mean(summ.w1$fF_2004$mean),
                    mean(summ.w1$fF.Hs_1969$mean), mean(summ.w1$fF.Hs_2004$mean)), 1),
                nsmall = 1)  )

summ$fF = paste0(summ$fF$mean,' (',summ$fF$lower,'-',summ$fF$upper,')')
summ$fF.Hs = paste0(summ$fF.Hs$mean,' (',summ$fF.Hs$lower,'-',summ$fF.Hs$upper,')')

summ.w2 <- pivot_wider(summ, 
                     id_cols = 'Site',
                     names_from = 'Year',
                     values_from = c('tot.obs','fF','fF.Hs')) 

summ <- rbind(data.frame(summ.w2), summ.tot)

summ <- summ[,1:5] # Drop Hs columns afterall

latex(
  summ,
  file='../tables/Paine-SiteSumm.tex',
  cgroup = c('', 'Observations', '\\% feeding'),
  n.cgroup = c(1, 2, 2),
  colheads = c('Site',
               '1968-9', '2004',
               '1968-9', '2004'),
  n.rgroup = c(5,1),
  rowname = NULL,
  label = 'tab:summ',
  center = 'centering',
  first.hline.double = FALSE,
  where = "!htbp",
  caption="Summary of Paine's 1968-9 and my 2004 feeding observations.  Observations refers to the total number of whelks inspected. \\% feeding refers to the proportion of observed whelks that were feeding. Parentheticals are the binomial confidence interval (95\\% coverage probability) calculated using the Wilson method."
)
#########################################################################
#########################################################################
#########################################################################
