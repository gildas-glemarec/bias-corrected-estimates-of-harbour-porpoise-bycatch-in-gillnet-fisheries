##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Porpoise bycatch mortality estimates ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Set to TRUE to discretize spatial coordinates to ices square
## and speed up the computations for testing.
lowSpatialResolution = FALSE

## NOTE: This code uses a custom branch of glmmTMB called "aggregate":
## https://github.com/glmmTMB/glmmTMB/tree/aggregate
##
##
## Install like this:
### https://stackoverflow.com/questions/52447227/how-do-i-install-an-r-package-under-another-name
# devtools::install_github("glmmTMB/glmmTMB/glmmTMB", ref="aggregate", 
#                          lib = 'H:/DTU Projects/Articles/Fidelity/lib1')
##
##
## It has new arguments for glmmTMB:::predict.glmmTMB : "aggregate" and "do.bias.correct"
## It lets us get confidence intervals and bias-correction for SUM of the predictions
## "aggregate" is a factor that predictions will be summed across.
## Note, use the (bias.correct) version of the mean

## LOAD LIBRARIES AND FUNCTIONS ####
gc()
par(mfrow = c(1,1))
Sys.setenv(LANG = "en") # force error messages to English
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Data management
library(tidyverse)
library(data.table)
library(lubridate)
## Statistical modelling and model validation
library(glmmTMB)
library(TMB)
library(DHARMa)
library(bbmle)
library(parameters)
library(performance) ## NB: remotes::install_github("easystats/parameters")
library(mapplots)
library(see)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## LOAD FUNCTIONS ----
`%notin%` <- Negate(`%in%`)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## LOAD DATA ----
work.data <- fread("work.data.csv",
                   header = TRUE,
                   sep = ';', dec = ',')
work.data.FD2 <- work.data #work.data[,!c('tot.landings')]
work.data.FD2$fid <- factor(work.data.FD2$fid)
work.data.FD2$oal <- factor(work.data.FD2$oal)
work.data.FD2$ices.area <- factor(work.data.FD2$ices.area)
work.data.FD2$icesrect <- factor(work.data.FD2$icesrect)
work.data.FD2$y <- factor(work.data.FD2$y)
work.data.FD2$m <- factor(work.data.FD2$m)
work.data.FD2$quarter <- factor(work.data.FD2$quarter)
work.data.FD2$f.mesh <- factor(work.data.FD2$f.mesh)
work.data.FD2$sp.group <- factor(work.data.FD2$sp.group)
work.data.FD2[, tot.landings := NULL]
work.data.FD2 <- work.data.FD2[complete.cases(work.data.FD2)]
## f.length: no sampling of <8m in the EM dataset -> merge with 8-10m
work.data.FD2 <- work.data.FD2 %>%
  dplyr::mutate(f.length = 
                  case_when(vessel.length < 10 ~ "<10m",
                            between(vessel.length, 10, 12) ~ "10-12m",
                            between(vessel.length, 12, 15) ~ "12-15m",
                            vessel.length > 15 ~ ">15m",
                            TRUE ~ as.character(vessel.length)))
work.data.FD2$f.length <- factor(work.data.FD2$f.length)

work.data.FD2[, population := fifelse(ices.area %notin% c("4.a","4.b","3.a.20"),
                                      'IDW_pop', 'NS_pop')]
work.data.FD2$population <- factor(work.data.FD2$population)

if(lowSpatialResolution){
  lonlat = sapply(work.data.FD2$icesrect,
                  DATRAS::icesSquare2coord,format="midpoint")
  work.data.FD2$lon = unlist(lonlat[1,])
  work.data.FD2$lat = unlist(lonlat[2,])
  work.data.FD2$pos <- numFactor(work.data.FD2$lon, work.data.FD2$lat)
  
} else {
  work.data.FD2$pos <- factor(work.data.FD2$pos)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## MODEL FIT ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bc.model <- glmmTMB(porpoise.p.FD.sq ~
                      f.length + f.mesh + log(netlength.p.FD.sq) + log(soak.p.FD.sq) +
                      (exp(pos + 0 | quarter)) + (exp(pos + 0 | y)) +
                      f.mesh:log(soak.p.FD.sq)
                    ,
                    family = nbinom2,
                    data = work.data.FD2)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Summarise best-fitting model
### See https://www.rdocumentation.org/packages/parameters/versions/0.19.0/topics/model_parameters
model_parameters(bc.model,
                 bootstrap = F, effects = "fixed",
                 # exponentiate = TRUE,
                 wb_component = TRUE)
# Fixed Effects

# Parameter                               | Coefficient |   SE |           95% CI |     z     |   df |      p
# --------------------------------------------------------------------------------------------------------
# (Intercept)                             |      -17.28 | 2.22 | [-21.64, -12.93] |     -7.78 | 6029 | < .001
# f.length [>15m]                         |        1.59 | 0.54 | [  0.53,   2.65] |      2.93 | 6029 | 0.003 
# f.length [10-12m]                       |        1.29 | 0.38 | [  0.55,   2.03] |      3.42 | 6029 | < .001
# f.length [12-15m]                       |        1.00 | 0.47 | [  0.08,   1.93] |      2.12 | 6029 | 0.034 
# f.mesh [>200mm]                         |        1.22 | 2.42 | [ -3.52,   5.95] |      0.50 | 6029 | 0.615 
# f.mesh [120-200mm]                      |        5.20 | 1.87 | [  1.53,   8.87] |      2.77 | 6029 | 0.006 
# netlength.p.FD.sq [log]                 |        0.74 | 0.13 | [  0.49,   1.00] |      5.71 | 6029 | < .001
# soak.p.FD.sq [log]                      |        1.50 | 0.50 | [  0.53,   2.47] |      3.03 | 6029 | 0.002 
# f.mesh [>200mm] * soak.p.FD.sq [log]    |   -4.46e-03 | 0.59 | [ -1.16,   1.15] | -7.59e-03 | 6029 | 0.994 
# f.mesh [120-200mm] * soak.p.FD.sq [log] |       -0.95 | 0.50 | [ -1.94,   0.03] |     -1.90 | 6029 | 0.057  
# 
# Uncertainty intervals (equal-tailed) and p-values (two-tailed) computed using a Wald z-distribution approximation.

# bc.model %>% 
#   model_parameters() %>% 
#   # model_parameters(exponentiate = TRUE) %>% 
#   plot()

# View(broom.mixed::tidy(bc.model))
# car::Anova(bc.model)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~## Check Model ##~~~~~~##
## Check for multicollinearity
# check_collinearity(bc.model)

# plot(check_collinearity(bc.model))

# check_model(bc.model)

## Check residuals(by simulations)
SimOut <- simulateResiduals(bc.model, 1000, plot = TRUE)
testDispersion(SimOut, type = "PearsonChisq", alternative = 'greater', plot = T)

plotResiduals(SimOut, form = work.data.FD2$d2shore.p.FD.sq)
plotResiduals(SimOut, form = work.data.FD2$depth.p.FD.sq)
plotResiduals(SimOut, form = work.data.FD2$depth)
plotResiduals(SimOut, form = work.data.FD2$y)
plotResiduals(SimOut, form = work.data.FD2$quarter)
plotResiduals(SimOut, form = work.data.FD2$f.length)
plotResiduals(SimOut, form = work.data.FD2$ices.area)
plotResiduals(SimOut, form = work.data.FD2$m)
plotResiduals(SimOut, form = work.data.FD2$soak.p.FD.sq)
plotResiduals(SimOut, form = work.data.FD2$mesh.p.FD.sq)
plotResiduals(SimOut, form = work.data.FD2$f.mesh)
