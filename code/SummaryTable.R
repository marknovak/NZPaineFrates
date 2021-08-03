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
summ <-
  dat %>% 
  filter(Site%in%focalSites) %>%
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
  where = "!htbp",
  caption="Summary of Paine's 1968-9 and my 2004 feeding observations.  Observerations refers to the total number of whelks inspected. \\% feeding refers to the proportion of observed whelks that were feeding. \\% feeding on \\emph{H. scobina} refers to the proportion of feeding whelks that were feeding on \\emph{Haustrum scobina}.  Parentheticals are the biomial confidence interval (95\\% coverage probability) calculated using the Wilson method."
  
)
#########################################################################
#########################################################################
#########################################################################
