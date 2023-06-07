##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Porpoise bycatch mortality model selection with AIC ####
## author: Gildas Glemarec (ggle@aqua.dtu.dk)
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
library(bbmle)
library(mapplots)
library(MuMIn)
library(performance)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## LOAD FUNCTIONS ----
`%notin%` <- Negate(`%in%`)

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

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
## Spatial autocorrelation? ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### Spatial correlograms  examine patterns of spatial autocorrelation 
##### in your data or model residuals. They show how correlated are pairs of 
##### spatial observations when you increase the distance (lag) between them

## Mantel cross-correlograms 
library(ncf)
fit1 <- correlog(x = work.data.FD2$lon.haul, 
                 y = work.data.FD2$lat.haul,
                 z = work.data.FD2$porpoise.p.FD.sq,
                 w = as.numeric(work.data.FD2$quarter),
                 increment = 2,
                 resamp = 0,
                 latlon = T)
plot(fit1)
fit2 <- correlog(x = work.data.FD2$lon.haul, 
                 y = work.data.FD2$lat.haul,
                 z = work.data.FD2$porpoise.p.FD.sq,
                 w = as.numeric(work.data.FD2$y),
                 increment = 2,
                 resamp = 0,
                 latlon = T)
plot(fit2)

## Variograms
# https://rdrr.io/cran/DHARMa/man/testSpatialAutocorrelation.html

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## MODEL SELECTION ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bc.model.no.var <- glmmTMB(porpoise.p.FD.sq ~ 1,
                           family = nbinom2,
                           data = work.data.FD2,
)

bc.model.full.NoRE <- glmmTMB(porpoise.p.FD.sq ~
                                # 0 +
                                ## Discrete
                                f.length +
                                f.mesh +
                                # population + ## DO NOT USE !!!
                                ## Continuous
                                log(netlength.p.FD.sq) +
                                log(soak.p.FD.sq) +
                                log(depth.p.FD.sq) #+
                              ## Interactions
                              # log(netlength.p.FD.sq) : quarter +
                              # log(soak.p.FD.sq) : quarter +
                              # log(depth.p.FD.sq): quarter  +
                              # ## Random effects and spatial autocorrelation
                              # (1|fid) +
                              # exp(pos + 0 | sp.group) +
                              # exp(pos + 0 | quarter) +
                              # exp(pos + 0 | y) +
                              # exp(pos + 0 | y:quarter)
                              ,
                              family = nbinom2,
                              data = work.data.FD2)
bc.model1 <- update(bc.model.full.NoRE, . ~ . 
                    + (1 | y)
                    + exp(pos + 0 | sp.group)
                    #  + exp(pos + 0 | quarter)
                    #  + exp(pos + 0 | y)
                    #  + exp(pos + 0 | y:quarter)
)
bc.model1.no.year <- update(bc.model1, . ~ . 
                            - (1 | y)
)
bc.model2 <- update(bc.model.full.NoRE, . ~ . 
                    + (1 | y)
                    + exp(pos + 0 | sp.group) 
                    + exp(pos + 0 | quarter) 
                    # + exp(pos + 0 | y) 
                    # + exp(pos + 0 | y:quarter)
)
bc.model2.no.year <- update(bc.model2, . ~ . 
                            - (1 | y)
)
bc.model3 <- update(bc.model.full.NoRE, . ~ .  
                    + (1 | y)
                    + exp(pos + 0 | sp.group) 
                    + exp(pos + 0 | quarter) 
                    + exp(pos + 0 | y) 
                    # + exp(pos + 0 | y:quarter)
)
bc.model3.no.year <- update(bc.model3, . ~ . 
                            - (1 | y)
)
bc.model4 <- update(bc.model.full.NoRE, . ~ . 
                    + (1 | y)
                    + exp(pos + 0 | sp.group) 
                    + exp(pos + 0 | quarter) 
                    + exp(pos + 0 | y) 
                    + exp(pos + 0 | y:quarter)
)
bc.model4.no.year <- update(bc.model4, . ~ . 
                            - (1 | y)
)
bc.model1.RE <- update(bc.model.full.NoRE, . ~ . 
                       + (1 | fid) 
                       + (1 | y)
                       + exp(pos + 0 | sp.group)
                       #  + exp(pos + 0 | quarter)
                       #  + exp(pos + 0 | y)
                       #  + exp(pos + 0 | y:quarter)
)
bc.model2.RE <- update(bc.model.full.NoRE, . ~ . 
                       + (1 | fid) 
                       + (1 | y)
                       + exp(pos + 0 | sp.group) 
                       + exp(pos + 0 | quarter) 
                       # + exp(pos + 0 | y) 
                       # + exp(pos + 0 | y:quarter)
)
bc.model3.RE <- update(bc.model.full.NoRE, . ~ . 
                       + (1 | fid) 
                       + (1 | y)
                       + exp(pos + 0 | sp.group) 
                       + exp(pos + 0 | quarter) 
                       + exp(pos + 0 | y) 
                       # + exp(pos + 0 | y:quarter)
)
bc.model4.RE <- update(bc.model.full.NoRE, . ~ . 
                       + (1 | fid) 
                       + (1 | y)
                       + exp(pos + 0 | sp.group) 
                       + exp(pos + 0 | quarter) 
                       + exp(pos + 0 | y) 
                       + exp(pos + 0 | y:quarter)
)
bc.model1.RE.no.year <- update(bc.model.full.NoRE, . ~ . + (1 | fid)
                               + exp(pos + 0 | sp.group)
                               #  + exp(pos + 0 | quarter)
                               #  + exp(pos + 0 | y)
                               #  + exp(pos + 0 | y:quarter)
)
bc.model2.RE.no.year <- update(bc.model.full.NoRE, . ~ . + (1 | fid)
                               + exp(pos + 0 | sp.group) 
                               + exp(pos + 0 | quarter) 
                               # + exp(pos + 0 | y) 
                               # + exp(pos + 0 | y:quarter)
)
bc.model3.RE.no.year <- update(bc.model.full.NoRE, . ~ . + (1 | fid)
                               + exp(pos + 0 | sp.group) 
                               + exp(pos + 0 | quarter) 
                               + exp(pos + 0 | y) 
                               # + exp(pos + 0 | y:quarter)
)
bc.model4.RE.no.year <- update(bc.model.full.NoRE, . ~ . + (1 | fid)
                               + exp(pos + 0 | sp.group) 
                               + exp(pos + 0 | quarter) 
                               + exp(pos + 0 | y) 
                               + exp(pos + 0 | y:quarter)
)

AICctab(bc.model.no.var, bc.model.full.NoRE,
        bc.model1, bc.model1.RE, bc.model1.no.year, bc.model1.RE.no.year,
        bc.model2, bc.model2.RE, bc.model2.no.year, bc.model2.RE.no.year,
        bc.model3, bc.model3.RE, bc.model3.no.year, bc.model3.RE.no.year,
        bc.model4, bc.model4.RE, bc.model4.no.year, bc.model4.RE.no.year
)
#                    dAICc df
# bc.model3            0.0 16
# bc.model3.RE         0.7 17
# bc.model2           18.5 14
# bc.model2.RE        19.0 15
# bc.model1           47.8 12
# bc.model1.RE        48.5 13
# bc.model.full.NoRE  95.2 10
# bc.model.no.var    424.7 2 
# bc.model4             NA 18
# bc.model4.RE          NA 19

# bc.model3.no.year      0.0 16
# bc.model3.RE.no.year   1.0 17
# bc.model3              2.0 17
# bc.model3.RE           3.0 18
# bc.model2             18.3 15
# bc.model2.RE          18.8 16
# bc.model2.no.year     19.5 14
# bc.model2.RE.no.year  20.0 15
# bc.model1             48.5 13
# bc.model1.RE          49.7 14
# bc.model1.no.year     50.2 12
# bc.model1.RE.no.year  50.9 13
# bc.model.full.NoRE   121.6 10
# bc.model.no.var      518.3 2 
# bc.model4               NA 19
# bc.model4.RE            NA 20
# bc.model4.no.year       NA 18
# bc.model4.RE.no.year    NA 19


## For all models with sp. autocorrelation, adding vessel ID as a random factor
## only marginally affect the fit (delta AIC <2). This suggests that the variation
## between vessels are already accounted for when taking into account variations
## in effort (fixed effect part) and in porpoise densities (spatial correlation
## part).

# summary(bc.model3)
## (...)
##   Overdispersion parameter for nbinom2 family (): 1.12 
## 
##   Conditional model:
##   Estimate Std. Error z value Pr(>|z|)    
##   (Intercept)            -14.2377     1.5577  -9.140  < 2e-16 ***
##   f.length>15m             1.6271     0.6054   2.688 0.007197 ** 
##   f.length10-12m           1.5094     0.4154   3.633 0.000280 ***
##   f.length12-15m           1.1018     0.4992   2.207 0.027303 *  
##   f.mesh>200mm             2.2515     0.6273   3.589 0.000331 ***
##   f.mesh120-200mm          1.6611     0.5613   2.959 0.003083 ** 
##   log(netlength.p.FD.sq)   0.6231     0.1372   4.540 5.62e-06 ***
##   log(soak.p.FD.sq)        0.6512     0.1122   5.802 6.54e-09 ***
##   log(depth.p.FD.sq)       0.3403     0.2035   1.673 0.094386 . 
##   
## All FE are significant but depth is borderline, so let's try without 
##   

bc.model5.no.year <- update(bc.model3.no.year , . ~ . - log(depth.p.FD.sq))
bc.model5.RE.no.year <- update(bc.model3.RE.no.year, . ~ . - log(depth.p.FD.sq))
AICctab(bc.model3.no.year, bc.model5.no.year, bc.model5.RE.no.year)
##              dAICc df
## bc.model3     0.0  16
## bc.model5.RE  0.7  16
## bc.model5     0.8  15
##   
## Delta AIC < 2 between all models so, using the parcimony principle, we 
## select the simplest model (i.e., bc.model5)
##   


## Can we reduce the model even more?
bc.model6 <- update(bc.model5.no.year, . ~ .  
                    - exp(pos + 0 | sp.group)
)
bc.model7 <- update(bc.model5.no.year, . ~ .  
                    - exp(pos + 0 | quarter)
)
bc.model8 <- update(bc.model5.no.year, . ~ .  
                    - exp(pos + 0 | y)
)
bc.model9 <- update(bc.model5.no.year, . ~ .
                    - exp(pos + 0 | sp.group)
                    - exp(pos + 0 | quarter)
)
bc.model10 <- update(bc.model5.no.year, . ~ .
                     - exp(pos + 0 | sp.group)
                     - exp(pos + 0 | y)
)

AICctab(bc.model5.no.year, bc.model6, bc.model7,
        bc.model8, bc.model9, bc.model10
)
##             dAICc df
## bc.model6   0.0  13
## bc.model3   1.8  16
## bc.model5   2.6  15
## bc.model10 20.5  11
## bc.model8  21.0  13
## bc.model7  26.8  13
## bc.model9  31.5  11

summary(bc.model6)
## Overdispersion parameter for nbinom2 family (): 1.12 
## 
## Conditional model:
##   Estimate Std. Error z value Pr(>|z|)    
##   (Intercept)            -13.3561     1.3445  -9.934  < 2e-16 ***
##   f.length>15m             1.8563     0.5752   3.228 0.001249 ** 
##   f.length10-12m           1.5224     0.3974   3.831 0.000127 ***
##   f.length12-15m           1.2731     0.4881   2.608 0.009106 ** 
##   f.mesh>200mm             2.3204     0.6258   3.708 0.000209 ***
##   f.mesh120-200mm          1.8069     0.5546   3.258 0.001123 ** 
##   log(netlength.p.FD.sq)   0.6172     0.1351   4.568 4.92e-06 ***
##   log(soak.p.FD.sq)        0.6545     0.1122   5.831 5.52e-09 ***
##

## Out of curiosity:
bc.model6.1 <- update(bc.model6, . ~ . + log(depth.p.FD.sq))
AICctab(bc.model6, bc.model6.1)
##             dAICc df
## bc.model6.1  0.0  14
## bc.model6    1.8  13
## Here again, adding depth as a FE only marginally improves the fit, so we
## decide to drop it for good

## Let's now also have a look at interactions
## Because depth is correlated to position (which is itself aggregated), I'd
## rather ignore the cross effects of depth, and focus on the crossed effects of
## mesh with soak, net length. I also need to test crossed soak * net length
bc.model11 <- update(bc.model6, . ~ .  
                     + f.mesh : log(soak.p.FD.sq)
                     + (1 | y)
)
bc.model11.no.year <- update(bc.model6, . ~ .  
                     + f.mesh : log(soak.p.FD.sq)
)
bc.model12 <- update(bc.model6, . ~ .  
                     + f.mesh : log(netlength.p.FD.sq)
)
bc.model13 <- update(bc.model6, . ~ .  
                     + log(netlength.p.FD.sq) : log(soak.p.FD.sq)
)

AICctab(bc.model6, bc.model11, bc.model12, bc.model13)
##            dAICc df
## bc.model11  0.0  15
## bc.model6   9.5  13
## bc.model13 10.0  14
## bc.model12 10.1  15

bc.model14 <- update(bc.model6, . ~ . 
                     + f.mesh : log(soak.p.FD.sq)
                     + log(netlength.p.FD.sq) : log(soak.p.FD.sq)
)
AICctab(bc.model11, bc.model14)
##            dAICc df
## bc.model11  0.0  15
## bc.model14  0.9  16

## OK, so now we have a winner, but what about fid, is it still useless?
bc.model11.RE <- update(bc.model11, . ~ . + (1 | fid) + (1 | y))
bc.model11.RE.no.year <- update(bc.model11.RE, . ~ . - (1 | y))
AICctab(bc.model11, bc.model11.RE, bc.model11.no.year, bc.model11.RE.no.year)
##               dAICc df
## bc.model11.RE  0.0  16
## bc.model11     0.8  15

## I think we can stop here: bc.model11 it is!!!
summary(bc.model11)
## Overdispersion parameter for nbinom2 family (): 1.13 
## 
## Conditional model:
##   Estimate Std. Error z value Pr(>|z|)    
##   (Intercept)                       -16.730631   2.295360  -7.289 3.13e-13 ***
##   f.length>15m                        1.839569   0.551931   3.333 0.000859 ***
##   f.length10-12m                      1.407714   0.390820   3.602 0.000316 ***
##   f.length12-15m                      1.045419   0.484746   2.157 0.031034 *  
##   f.mesh>200mm                        1.091791   2.505684   0.436 0.663036    
##   f.mesh120-200mm                     5.464620   1.946756   2.807 0.005000 ** 
##   log(netlength.p.FD.sq)              0.656541   0.135473   4.846 1.26e-06 ***
##   log(soak.p.FD.sq)                   1.567225   0.525472   2.983 0.002859 ** 
##   f.mesh>200mm:log(soak.p.FD.sq)      0.006736   0.618667   0.011 0.991312    
##   f.mesh120-200mm:log(soak.p.FD.sq)  -1.086322   0.533172  -2.037 0.041603 *
formula(bc.model11)
# porpoise.p.FD.sq ~ 
#   f.length + f.mesh + log(netlength.p.FD.sq) + log(soak.p.FD.sq) + 
#   (exp(pos + 0 | quarter)) + (exp(pos + 0 | y)) +
#   f.mesh:log(soak.p.FD.sq)
bc.model11.1 <- update(bc.model11, . ~ . - (exp(pos + 0 | quarter)))
bc.model11.2 <- update(bc.model11, . ~ . - (exp(pos + 0 | y)))
bc.model11.3 <- update(bc.model11, . ~ . 
                       - (exp(pos + 0 | quarter))
                       - (exp(pos + 0 | y)))
AICctab(bc.model11, bc.model11.1, bc.model11.2, bc.model11.3)
##              dAICc df
## bc.model11     0.0 15
## bc.model11.2  15.6 13
## bc.model11.1  34.0 13
## bc.model11.3 110.8 11

## 
## LIST ALL THE HYPOTHESES THAT WERE JUST TESTED HERE
## 
## 

## Create an AIC model selection table with model formulas included ####
modelTABLE <- model.sel(
  bc.model.no.var, bc.model.full.NoRE,
  bc.model1, bc.model1.RE, bc.model1.no.year, bc.model1.RE.no.year,
  bc.model2, bc.model2.RE, bc.model2.no.year, bc.model2.RE.no.year,
  bc.model3, bc.model3.RE, bc.model3.no.year, bc.model3.RE.no.year,
  bc.model4, bc.model4.RE, bc.model4.no.year, bc.model4.RE.no.year,
  bc.model6, bc.model5.no.year, bc.model5.RE.no.year,
  bc.model7,
  bc.model8,
  bc.model9,
  bc.model10,
  bc.model11, bc.model11.RE, bc.model11.no.year, bc.model11.RE.no.year,
  # bc.model11.1, bc.model11.2, bc.model11.3,
  bc.model12,
  bc.model13,
  bc.model14
)
View(modelTABLE)  ##No AIC values, just AICc, no R-squared, and model name (i.e, a~b) not present
# modelTABLE <- model.sel(modelTABLE, rank = AIC)
write.csv2(modelTABLE, 'modelTABLE.csv')


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Examples of bias-corrected predictions ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Load glmmTMB incl. aggregate function
# detach("package:glmmTMB", unload = TRUE)
# library(glmmTMB, lib.loc = 'path to the branched version of glmmTMB')

## Data.for.prediction is a dataset merging logbook/salesnotes data from Sweden 
## and Denmark and aggregating the effort to ICES statistical square level

#### 1. Aggregate by year #### 

sumpred.y <- predict(bc.model,
                     newdata=Data.for.prediction,
                     allow.new.levels=TRUE,
                     aggregate=factor(Data.for.prediction$y),
                     se.fit=TRUE,
                     do.bias.correct=TRUE,
                     type="response")
## CV = sd / mean (bias correct version), which is approx SD of log mean.
sdlogpred <- sumpred.y[,4] / sumpred.y[,3]

## Gather the results in a data frame:

results.year = data.frame(year = levels(Data.for.prediction$y), 
                                    est = sumpred.y[,3], 
                                    cilow = exp( log(sumpred.y[,3]) - 2*sdlogpred), 
                                    cihigh = exp( log( sumpred.y[,3]) + 2*sdlogpred))
fwrite(results.year,
       'results.year.bias.corr.csv',
       sep = ';', dec = ',')

#### 2. Aggregate by quarter/area ####
## Mean (and 95% CI) quarterly fleet-wide harbour porpoise bycatch estimates 
gc()
sumpred.q <- lapply(split(Data.for.prediction,
                          by = 'ices.area',
                          flatten = TRUE),
                    function(x){predict(bc.model,
                                        newdata=x,
                                        allow.new.levels=TRUE,
                                        aggregate=factor(x$quarter),
                                        se.fit=TRUE,
                                        do.bias.correct=TRUE,
                                        type="response")}
)
# ## There are 11 years of data aggregated here, so:
sumpred.q <- lapply(sumpred.q,function(x){x / n_distinct(Data.for.prediction$y)}) 

## CV = sd / mean (bias correct version), which is approx SD of log mean.
rm(sdlogpred)
sdlogpred <- lapply(sumpred.q,
                    function(x){as.data.table(x[,4] / x[,3])})
sumpred.q <- 
  map2(sumpred.q, sdlogpred, ~cbind(.x, setNames(.y, names(.x))))

results.quarter <- lapply(sumpred.q,
                                    function(x){
                                      data.frame(quarter = levels(Data.for.prediction$quarter), 
                                                 est = x[,3], 
                                                 cilow = exp( log(x[,3]) - 2*x[,5]), 
                                                 cihigh = exp( log( x[,3]) + 2*x[,5]))
                                    }
) 
results.quarter <- rbindlist(results.quarter, idcol = TRUE)
setnames(results.quarter,
         old = c('.id','Est...bias.correct.',
                 'Est...bias.correct..1','Est...bias.correct..2'),
         new = c('ices.area','est.bc',
                 'lwr.bc','upr.bc'))

fwrite(results.quarter,
       'results.quarter.bias.corr.csv',
       sep = ';', dec = ',')