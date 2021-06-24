###############################################################################
###############################################################################
###############################################################################
rm(list = ls())
options(stringsAsFactors = FALSE)

library(sfsmisc)
library(corrplot)
############################

# All available 1968-69 and 2004 data
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

############################
# Remove non-feeding
n_all_69 <- n_all_69[!c(rownames(n_all_69)%in%'Not Feeding'),]
n_all_04 <- n_all_04[!c(rownames(n_all_04)%in%'Not Feeding'),]

# Remove surveys with less than nMin (feeding + non-feeding) observations 
nMin <- 20
n_all_69 <- n_all_69[, apply(n_all_69, 2, sum, na.rm=TRUE) >= nMin]
n_all_04 <- n_all_04[, apply(n_all_04, 2, sum, na.rm=TRUE) >= nMin]

# Remove relic species
n_all_69 <- n_all_69[apply(n_all_69, 1, sum, na.rm=TRUE) > 0 ,]
n_all_04 <- n_all_04[apply(n_all_04, 1, sum, na.rm=TRUE) > 0 ,]

# Convert counts to site-specific proportions
p_all_69 <- apply(n_all_69, 2, function(x) {x/sum(x, na.rm=TRUE)} )
p_all_04 <- apply(n_all_04, 2, function(x) {x/sum(x, na.rm=TRUE)} )

# Pairwise correlations
R2.69 <- cor(p_all_69, use = 'pairwise.complete.obs') ^ 2
R2.log69 <- cor(log10(p_all_69), use = 'pairwise.complete.obs') ^ 2
R2.04 <- cor(p_all_04, use = 'pairwise.complete.obs') ^ 2
R2.log04 <- cor(log10(p_all_04), use = 'pairwise.complete.obs') ^ 2

corrplot(R2.69, 
         method = 'ellipse',
         type = 'lower', 
         is.cor = FALSE, 
         tl.col = 'black',
         tl.pos = 'l',
         cl.lim = c(0,1))
corrplot(R2.log69, 
         add = TRUE,
         method = 'ellipse', 
         type = 'upper', 
         diag = FALSE, 
         tl.pos = 'n',
         is.cor = FALSE, 
         cl.pos = 'n',
         cl.lim = c(0,1))


corrplot(R2.04, 
         method = 'ellipse',
         type = 'lower', 
         is.cor = FALSE, 
         tl.col = 'black',
         tl.pos = 'l',
         cl.lim = c(0,1))
corrplot(R2.log04, 
         add = TRUE,
         method = 'ellipse', 
         type = 'upper', 
         diag = FALSE, 
         tl.pos = 'n',
         is.cor = FALSE, 
         cl.pos = 'n',
         cl.lim = c(0,1))

# High values come from there being only 2 or 3 prey species in common.
# Concl: There are just not enough observations per site to get good estimates
# of spatial correlations.
###############################################################################
###############################################################################
###############################################################################