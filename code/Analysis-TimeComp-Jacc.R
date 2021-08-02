#########################################################################
#########################################################################
#########################################################################
# Calculate and plot between-period Jaccard similarities
# of feeding observations and prey abundances
#########################################################################
#########################################################################
rm(list = ls())
options(stringsAsFactors = F)
source('JaccardIndices.R')

library(dplyr)
library(tidyr)
library(Hmisc) # for LaTeX table export
options(xdvicmd='open')

############################
dat <- read.csv('../data/derived/NZ-1969_2004-tab_Summarized.csv')
focalSites <- read.csv("../data/derived/NZ-1969_2004-Site-focal.csv")$x

############################
dat <- filter(dat, Site%in%focalSites & Prey!='Not Feeding' )
############################
# Time-periods side-by-side
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Feeding survey sites
datw <- pivot_wider(data = dat, 
                    id_cols = c('Site','Prey'),
                    names_from = 'Year',
                    values_from = c('n.obs','N.mean')) %>%
  mutate(n.obs_1969 = replace_na(n.obs_1969, 0),
         n.obs_2004 = replace_na(n.obs_2004, 0),
         N.mean_1969 = replace_na(N.mean_1969, 0),
         N.mean_2004 = replace_na(N.mean_2004, 0))


#############################
J.nobs.abund <- 
  datw %>% group_by(Site) %>%
  summarise(Jclass.nobs = round(J.class(n.obs_1969, n.obs_2004), 2),
            Jabd.nobs = round(J.abd(n.obs_1969, n.obs_2004), 2),
            Jest.nobs = round(J.est(n.obs_1969, n.obs_2004), 2),
            Jclass.N = round(J.class(N.mean_1969, N.mean_2004), 2),
            Jabd.N = round(J.abd(N.mean_1969, N.mean_2004), 2),
            Jest.N = round(J.est(N.mean_1969, N.mean_2004), 2)) %>%
  replace_na(list(Jclass.N = "-", Jabd.N = "-", Jest.N = "-") )

J.nobs.abund



latex(
  J.nobs.abund,
  file='../tables/Paine-Jacc-nobsabund.tex',
  cgroup = c('', 'Feeding observations', 'Prey abundances'),
  n.cgroup = c(1, 3, 3),
  colheads = c('Site','$J_{class}$', '$J_{abd}$','$J_{est}$',
                      '$J_{class}$', '$J_{abd}$','$J_{est}$'),
  rowname = NULL,
  label = 'tab:Jnobsabund',
  center = 'centering',
  first.hline.double = FALSE,
  caption="The between time period similarity of
  \\textit{Haustrum haustorium}'s apparent diet 
  -- at the sites where Paine and I performed either feeding surveys only
  or both feeding and abundance surveys --
  as quantified by 
  the incidence-based Jaccard index ($J_{class}$), as well as
  the abundance-based Jaccard index ($J_{abd}$) and
  the estimator for the abundace-based Jaccard index ($J_{est}$) 
  of \\citet{Chao:2005tf}.
  Unlike the correlation and distance-based comparisons of the main text which
  included only prey species which both Paine and I observed 
  \\textit{H. haustorium} feeding on at a given site,
  the feeding observation comparisons include prey species which only one of us observed
  and the abundance surveys include prey species which only one of us observed as well 
  as additional (non-prey) mobile species which Paine or I observed."
)


##############################################################################
##############################################################################
##############################################################################


