### LB Updated MRSea Model Script for TAST JAE Paper
### 18 November 2023

## read in data
padat <- read.csv("pa_data.csv")

## load libraries
library(tidyverse)
require(car)
library(lawstat)
library(MRSea)
library(mgcv)
library(splines)
library(geepack)

## format data
df <- padat %>%
  select(pa, treatment, obs, x, y, session_id, id, jul_day, start_hr, obs, fish)
df$pa <- as.numeric(df$pa) 
df$treatment <- as.factor(df$treatment)
df$obs <- as.factor(df$obs)

names(df)[4] <- "x.pos"
names(df)[5] <- "y.pos"

#  STEP 1: check for collinearity between covariates in glm model ---------
fullmod_linear <- glm(pa ~ treatment + jul_day + start_hr + obs + fish + x.pos + y.pos + x.pos:y.pos,
                      na.action = "na.fail", family = "binomial", data = df)

vif(fullmod_linear)
AIC(fullmod_linear)
summary(fullmod_linear)


# STEP 2: make splineParams object with knots at the mean of each  --------
splineParams<-makesplineParams(data=df, varlist=c("x.pos", "y.pos", "fish","jul_day", "start_hr"))

# fit new full mod with splines binomial dist
fullmod <- glm(pa ~ treatment + obs + bs(x.pos, knots = splineParams[[2]]$knots) 
               + bs(y.pos, knots = splineParams[[3]]$knots)
               + bs(fish,knots = splineParams[[4]]$knots)
               + bs(jul_day,knots = splineParams[[5]]$knots)
               + bs(start_hr, knots =splineParams[[6]]$knots),
               na.action = "na.fail", family = "binomial", data = df)


# STEP 3: check for correlated residuals in discrete covariates --------

varList <- c("y.pos", "x.pos", "fish", "jul_day", "start_hr")

# ACF plot with blocks as survey sesh, day hour
df$day_hour_as_block <- as.factor(paste(df$jul_day, df$start_hr,sep = ""))
ps <- par(ask = FALSE, mfrow = c(1,1))
runACF(df$day_hour_as_block, fullmod, store = F) #some positive AC, use this as our cor block


# STEP 4: plot cumulative residuals for model to check over or und --------

par(ask = F, mfrow = c(3,3))
plotCumRes(fullmod, varlist= varList, splineParams) 

# compare with model that just has linear terms
par(ask = F, mfrow = c(3,3))
plotCumRes(fullmod_linear, varlist= c(varList), splineParams) 

# STEP 5: Run SALSA1D -----------------------------------------------------

# Create the variable "response" needed by SALSA
df$response <- df$pa

## Set initial model without spline-based terms
initialModel <- glm(response ~ treatment, family=binomial, data=df)

## Set SALSA arguments
factorList <- c("obs", "treatment")
varList <- c("x.pos", "y.pos", "fish", "jul_day", "start_hr")
#create salsa list 
salsa1DList <- list(fitnessMeasure="AIC", minKnots_1d=rep(1, 5),
                    maxKnots_1d=rep(5, 5), startKnots_1d=rep(1, 5),
                    degree=rep(2, 5), maxIterations=100,
                    gaps=rep(0, 5))

## Run SALSA with removal
set.seed(53195)
salsa1D_RT <- MRSea::runSALSA1D(initialModel, salsa1DList, varList,
                                factorList, varlist_cyclicSplines=NULL,
                                splineParams=NULL, datain=df,
                                suppress.printout=TRUE, removal=T,
                                panelid=NULL)
bestModel1D<- salsa1D_RT$bestModel


# STEP 6: Select initial knot locations using a space-filling desin --------

# Next we add a two dimensional smooth of geographic coordinates (s(x.pos,
# y.pos))# To test for a redistribution of animals with a treatment effect we fit an
# interaction term between the smooth of coordinates and impact.
# SALSA is used to determine spatially adaptive knot locations for this
# smooth term.

## create knotgrid
knotgrid <- MRSea::getKnotgrid(coordData=cbind(df$x.pos, df$y.pos),
                               numKnots=300, plot=FALSE)


# STEP 7: Create distance matrix ------------------------------------------
# i.e. distance between each data points and knots)

## make distance matrices for datatoknots and knottoknots
distMats <- makeDists(cbind(df$x.pos, df$y.pos), na.omit(knotgrid))

## create sequence of radii
r_seq <- getRadiiChoices(numberofradii= 8, distMats$dataDist, basis = "gaussian")

# ### STEP 8: RUN SALSA 2D ------------------------------------------------
salsa2DList6 <- list(fitnessMeasure="AIC", knotgrid=knotgrid,
                     startKnots=6, minKnots=2, maxKnots=20,
                     gap=1, r_seq = r_seq, interactionTerm = "treatment")


# splineParams must be an object in workspace
## update splineParams with the SALSA1D results
splineParams <- salsa1D_RT$splineParams


set.seed(53195)
salsa2D_6 <- MRSea::runSALSA2D(bestModel1D, salsa2DList6, d2k=distMats$dataDist,
                               k2k=distMats$knotDist, splineParams=NULL,
                               tol=0, chooserad=F, panels=NULL,
                               suppress.printout=TRUE)
## Store "best" model (tested several startign knots, 6 was best)
baseModel <- salsa2D_6$bestModel

## Summary of results
summary(baseModel)
aov(baseModel, test="Chisq")

## calculate Deviance
library(modEvA)
Dsquared(model = baseModel, adjust = TRUE)
Dsquared(model = baseModel, adjust = F)

## update spline parameter object
splineParams <- baseModel$splineParams

# check for residual correlation
runs.test(residuals(baseModel, type = "pearson")) #still neg stat, still tiny p


# STEP 9: fit as GEE--------
# significant positive resid correlation therefore refit as GEE
# the data must be ordered by block (which this is)
# and the blockid must be numeric

## specify parameters for local radial:
radiusIndices <- splineParams[[1]]$radiusIndices
dists <- splineParams[[1]]$dist
radii <- splineParams[[1]]$radii
aR <- splineParams[[1]]$knotPos ## Index of knot locations. The index contains numbers selected by SALSA from 1 to the number of legal knot locations  

## update model in workspace with parameters for spatial smooth (above)
baseModel <- update(baseModel, . ~ .)
df$day_hour_as_block <- as.character(df$day_hour_as_block)

## update spline parameter object
splineParams <- baseModel$splineParams

## Re-fit the chosen model as a GEE (based on SALSA knot placement) 
geeModel <- geepack::geeglm(formula(baseModel), family = binomial, data = df, id = day_hour_as_block)
summary(geeModel)

## GEE p-valuesd
anova.gamMRSea(geeModel) 
summary(geeModel)

# Remove interaction to compare
geeModel$formula 

# grab model from above ^ & remove interaction
noint.model <- response ~  treatment + bs(x.pos, knots = splineParams[[2]]$knots,degree = splineParams[[2]]$degree, Boundary.knots = splineParams[[2]]$bd) +
  bs(y.pos, knots = splineParams[[3]]$knots, degree = splineParams[[3]]$degree,Boundary.knots = splineParams[[3]]$bd) + 
  bs(fish, knots = splineParams[[4]]$knots, degree = splineParams[[4]]$degree, Boundary.knots = splineParams[[4]]$bd) + 
  bs(jul_day, knots = splineParams[[5]]$knots, degree = splineParams[[5]]$degree,Boundary.knots = splineParams[[5]]$bd) + 
  bs(start_hr,knots = splineParams[[6]]$knots, degree = splineParams[[6]]$degree, Boundary.knots = splineParams[[6]]$bd) + 
  LRF.g(radiusIndices, dists, radii, aR)

nointModel <- geepack::geeglm(formula(noint.model), data = df, family = "binomial",id = day_hour_as_block)

## compare p-values
anova.gamMRSea(nointModel)
anova.gamMRSea(geeModel) 

# The interaction is significant in geeglm, but the treatment term in the model without
# the interaction is not significant then there has been a re-distribution of
# animals but no overall change!!! 

# Partial residual plots 
par(mfrow = c(2, 3))
runPartialPlots(model = geeModel, data = df, factorlist.in = c("treatment"), save =T) #, factorlist.in = factorList, varlist.in = varList, data = df, showKnots = T)
runPartialPlots(model = baseModel, data = df, varlist.in = c("fish"), save =T)
runPartialPlots(model = baseModel, data  = df, varlist.in = c("start_hr"), save =T)
runPartialPlots(model = baseModel, data  = df, varlist.in = c("x.pos"), save =T)
runPartialPlots(model = baseModel, data  = df, varlist.in = c("y.pos"), save =T)
runPartialPlots(baseModel, data = df, varlist.in = c("jul_day"), save =T)


# Step 10: Diagnostics ----------------------------------------------------

# create observed vs fitted and fitted vs residual plots
par(mfrow = c(1, 1))
runDiagnostics(geeModel, save = T) 

# How to interpret these with binary output
install.packages("arm")
arm::binnedplot(fitted(geeModel), 
                residuals(geeModel, type = "response"), 
                nclass = NULL, 
                xlab = "Expected Values", 
                ylab = "Average residual", 
                main = "Binned residual plot", 
                cex.pts = 0.8, 
                col.pts = 1, 
                col.int = "gray")

# Cumulative Residuals and runs profile 
plotCumRes(geeModel, varlist= c("x.pos", "y.pos"))
plotRunsProfile(geeModel, varlist= c("x.pos", "y.pos"), save = T) 
plotCumRes(geeModel, varlist= c( "start_hr"))
plotRunsProfile(geeModel, varlist= c("jul_day", "start_hr"), save = T) 

## COVRATIO and PRESS Statistic --- This will take 276 minutes or 4.6 hours... run over night
timeInfluenceCheck(geeModel, df$day_hour_as_block, dists, splineParams)
# influence plots (covratio and press statistics)
influence <- runInfluence(geeModel, df$day_hour_as_block,  dists, splineParams)

## RAW residual plots
resids <- fitted(geeModel) - df$pa
dims <- getPlotdimensions(df$x.pos, df$y.pos, 10, 10)
#par(mfrow = c(1, 2), mar = c(5, 7, 5, 5))
quilt.plot(df$x.pos[df$treatment == "off"], df$y.pos[df$treatment == "off"], 
           resids[df$treatment == "off"], asp = 1, ncol = dims[2], nrow = dims[1],
           zlim = c(-2.2, 2.2), main = "off", xlab  = "x position (m)", ylab = "y position (m)", add.legend = F)
par(mar = c(5, 5, 5, 7))
quilt.plot(df$x.pos[df$treatment == "on"], df$y.pos[df$treatment == "on"], 
           resids[df$treatment == "on"], asp = 1, ncol = dims[2], nrow = dims[1],
           zlim = c(-2.2, 2.2), main = "on", xlab  = "x position (m)", 
           ylab = "y position (m)", legend.lab = "Magnitude of Residuals", legend.shrink = 0.7, legend.args=list(cex=1.5, side=1))

## Influence plots
timeInfluenceCheck(geeModel, df$block_as_session, dists, splineParams)
# influence plots (covratio and press statistics) this will take 6hrs
influence <- runInfluence(geeModel, df$block_as_session, dists, splineParams)


# Step 11: Predictions --------------------------------------------------------------

# loading the prediction grid data
df$offset <- rep(100, length(df$pa)) #add area column
predictData <- df

# create the distance matrix for predictions
g2k <- makeDists(cbind(predictData$x.pos, predictData$y.pos),na.omit(knotgrid), knotmat = FALSE)$dataDist

# use baseModel to make predictions 
predslink<- predict.gamMRSea(object = baseModel, predictData, type= "link")

# reversing the logit-link to convert predictions back to the response scale
preds <- plogis(predslink)
length(preds)
predictData$preds <- preds

max_on <- max(predictData$preds[predictData$treatment == "ON"])
mean_on <- mean(predictData$preds[predictData$treatment =="ON"])
max_ofF <- max(predictData$preds[predictData$treatment =="OFF"])
mean_off <- mean(predictData$preds[predictData$treatment =="OFF"])

### Step 19 visualising predictions
# plotting the predictions for before and after impact
# get the plot dimensions. We know each cell is 10x10m
dims<-getPlotdimensions(x.pos=predictData$x.pos, predictData$y.pos,
                        segmentWidth=10, segmentLength=10)
par(mfrow=c(1,2))#, mar=c(5,5,3,5))
quilt.plot(predictData$x.pos[predictData$treatment=="OFF"],
           predictData$y.pos[predictData$treatment=="OFF"],
           preds[predictData$treatment=="OFF"], asp=1, nrow=dims[1], ncol=dims[2], 
           zlim=c(0, maxlim = 0.3), main = "TAST OFF", add.legend = F,ylab = "y position", xlab = "x position")

quilt.plot(predictData$x.pos[predictData$treatment=="ON"],
           predictData$y.pos[predictData$treatment=="ON"], preds[predictData$treatment=="ON"],
           asp=1,nrow=dims[1], ncol=dims[2], zlim=c(0, maxlim = 0.3), legend.lab = "Probability of Occurance", horizontal = F, 
           main = "TAST ON",  ylab = "y position", xlab = "x position")


## Step 20: bootstrap CI ------------------------------------------------------
# do the bootstrap at 999 when you get the chance

bootPreds <- do.bootstrap.cress.robust(model.obj = baseModel, 
                                       predictionGrid = predictData,
                                       g2k = g2k,
                                       B=1000,
                                       robust = T) # try this at 999

cis <- makeBootCIs(bootPreds)
predictData$LCL <- cis[,1]
predictData$UCL <- cis[,2]
predictData$mean <- rowMeans(bootPreds)
predictData$median <- apply(bootPreds,1,median)
head(predictData)

# bootstrap plots showing upper and lower limits

# make percentile confidence intervals
cison <- makeBootCIs(bootPreds[predictData$treatment=="ON",])

cisoff <- makeBootCIs(bootPreds[predictData$treatment=="OFF",])

## Step 21: visualising boot CI
par(mfrow=c(2,2), mar=c(3,4,3,3)) 

#low OFFs
quilt.plot(predictData$x.pos[predictData$treatment=="OFF"],
           predictData$y.pos[predictData$treatment=="OFF"],
           cisoff[,1], asp=1, nrow=dims[1], ncol=dims[2], 
           zlim=c(0, maxlim = 0.4), main = " Lower CI  - TAST OFF", add.legend = F)
par(mar=c(3,1,3,6))
# up OFF
quilt.plot(predictData$x.pos[predictData$treatment=="OFF"],
           predictData$y.pos[predictData$treatment=="OFF"], cisoff[,2],
           asp=1,nrow=dims[1], ncol=dims[2], zlim=c(0, maxlim = 0.4), 
           main = "Upper CI - TAST OFF")
par(mar=c(3,4,3,3))
#low ONs
quilt.plot(predictData$x.pos[predictData$treatment=="ON"],
           predictData$y.pos[predictData$treatment=="ON"],
           cison[,1], asp=1, nrow=dims[1], ncol=dims[2], 
           zlim=c(0, maxlim = 0.4), main = " Lower CI  - TAST ON", add.legend = F)
par(mar=c(3,1,3,6))
# up ON
quilt.plot(predictData$x.pos[predictData$treatment=="ON"],
           predictData$y.pos[predictData$treatment=="ON"], cison[,2],
           asp=1,nrow=dims[1], ncol=dims[2], zlim=c(0, maxlim = 0.4), 
           main = "Upper CI - TAST ON")

## Step 22: Identifying Differences 
differences <- getDifferences(beforePreds = bootPreds[predictData$treatment=="OFF",],
                              afterPreds = bootPreds[predictData$treatment=="ON",])

# Step 23: Visualising differences
# The median for each after - before difference
mediandiff <- differences$mediandiff

# The marker for each after - before difference:
# positive ('1') and negative ('-1') significant differences
marker <- differences$significanceMarker

#plot positive differences
par(mfrow = c(1, 1), mar = c(3, 6, 3, 3))
quilt.plot(predictData$x.pos[predictData$treatment=="OFF"],
           predictData$y.pos[predictData$treatment=="OFF"],
           mediandiff, asp = 1, nrow = 7, ncol = 9, add.legend = T, xlab = "x", 
           y.lab = "y", legend.lab = "Pr(occurance|Off) - Pr(occurance|On) ")

# add + depending on significance of cells. Just
# requires one significance out of all to be allocated
points(predictData$x.pos[predictData$treatment=="OFF"][marker==1],
       predictData$y.pos[predictData$treatment=="OFF"][marker==1], pch="+",
       col="black", cex=1)

# location of observer
points(0,0,pch= "*", col="black", cex=3)

#select pred dat off
pred.datoff <- predictData[predictData$treatment=="OFF",]
#randomly sample from each level of id pred off and select one from each val
predoff.samp <- pred.datoff[tapply(1:nrow(pred.datoff), pred.datoff$id, sample, 1),]
#check there are 235
nrow(predoff.samp)

#select pred dat on
pred.daton <- predictData[predictData$treatment=="ON",]
#randomly sample from each level of id pred off and select one from each val
predon.samp <- pred.daton[tapply(1:nrow(pred.daton), pred.daton$id, sample, 1),]
#check
nrow(predon.samp)

#now add points to graph
points(predictData$x.pos[predictData$treatment=="OFF"][marker==(-1)],
       predictData$y.pos[predictData$treatment=="OFF"][marker==(-1)], col="white",
       cex=1)

points(predictData$x.pos[predictData$treatment=="OFF"][marker==(0)],
       predictData$y.pos[predictData$treatment=="OFF"][marker==(0)], col="black",
       cex=.5) 
# location of observer
points(0,0,pch= "*", col="green", cex=3)

# Find linear distance between tast and predicted peak density for ON/OFF--------
## calculate distance from 0 for each cell
predictData$lin_d2_TAST <- sqrt((predictData$x.pos^2)+(predictData$y.pos^2))

## OFF
off_preds <- predictData[predictData$treatment == "OFF",]

ggplot(off_preds) +
  theme_classic()+
  stat_smooth(aes(lin_d2_TAST, mean), se =F, color = "black")+
  stat_smooth(aes(lin_d2_TAST, LCL), se =F, color = "red", size = 0.5, linetype = 4) +
  stat_smooth(aes(lin_d2_TAST, UCL), se =F, color = "red", size = 0.5,linetype = 4) +
  labs(x = "Distance to TAST (m)", y = "Pr(Presence)")
  
off_gam <- gam(mean ~ s(lin_d2_TAST, bs = "cs"), data = predictData[predictData$treatment == "OFF",] )
off_gam_lcl <- gam(LCL ~ s(lin_d2_TAST, bs = "cs"), data = predictData[predictData$treatment == "OFF",] )
off_gam_ucl <- gam(UCL ~ s(lin_d2_TAST, bs = "cs"), data = predictData[predictData$treatment == "OFF",] )


## ON
on_preds <- predictData[predictData$treatment == "ON",]

ggplot(on_preds) +
  theme_classic()+
  stat_smooth(aes(lin_d2_TAST, mean), se =F, level = 0.95, color = "black")+
  stat_smooth(aes(lin_d2_TAST, LCL), se =F, color = "red", size = 0.5, linetype = 4) +
  stat_smooth(aes(lin_d2_TAST, UCL), se =F, color = "red", size = 0.5,linetype = 4) +
  labs(x = "Distance to TAST (m)", y = "Pr(Presence)")

on_gam <- gam(mean ~ s(lin_d2_TAST, bs = "cs"), data = predictData[predictData$treatment == "ON",] )
on_gam_lcl <- gam(LCL ~ s(lin_d2_TAST, bs = "cs"), data = predictData[predictData$treatment == "ON",] )
on_gam_ucl <- gam(UCL ~ s(lin_d2_TAST, bs = "cs"), data = predictData[predictData$treatment == "ON",] )


# Make predictions using the model
x <- c(25, 50, 100, 150)

predictions_on <- predict.gam(on_gam, newdata = data.frame(lin_d2_TAST = x), type = "response")
lcl_preds_on <- predict.gam(on_gam_lcl, newdata = data.frame(lin_d2_TAST = x), type = "response")
ucl_preds_on <- predict.gam(on_gam_ucl, newdata = data.frame(lin_d2_TAST = x), type = "response")

predictions_off <- predict.gam(off_gam, newdata = data.frame(lin_d2_TAST = x), type = "response")
lcl_preds_off <- predict.gam(off_gam_lcl, newdata = data.frame(lin_d2_TAST = x), type = "response")
ucl_preds_off <- predict.gam(off_gam_ucl, newdata = data.frame(lin_d2_TAST = x), type = "response")

# make tables
onpred_tab <- data.frame(distance = x, pr_occurance = predictions_on, LCL = lcl_preds_on, UCL = ucl_preds_on)
offpred_tab <- data.frame(distance = x, pr_occurance = predictions_off, LCL = lcl_preds_off, UCL = ucl_preds_off)
library("gt")
gt(data.frame(onpred_tab))
gt(data.frame(offpred_tab))

## Double plot
ggplot(predictData) +
  theme_classic()+
  stat_smooth(aes(lin_d2_TAST, mean), se =F, color = "black")+
  stat_smooth(aes(lin_d2_TAST, LCL), se =F, color = "red", size = 0.5, linetype = 4) +
  stat_smooth(aes(lin_d2_TAST, UCL), se =F, color = "red", size = 0.5,linetype = 4) +
  labs(x = "Distance to TAST (m)", y = "Pr(Presence)") +
  facet_wrap(~treatment)


#### try to calculate difference between two gams
new_data <- runif(n = 1000, min = 0, max = 250)

Dpredictions_on <- predict.gam(on_gam, newdata = data.frame(lin_d2_TAST = new_data, treatment = "TAST ON"), type = "response")
Dlcl_preds_on <- predict.gam(on_gam_lcl, newdata = data.frame(lin_d2_TAST = new_data, treatment = "TAST ON"), type = "response")
Ducl_preds_on <- predict.gam(on_gam_ucl, newdata = data.frame(lin_d2_TAST = new_data, treatment = "TAST ON"), type = "response")

Dpredictions_off <- predict.gam(off_gam, newdata = data.frame(lin_d2_TAST = new_data, treatment = "TAST OFF"), type = "response")
Dlcl_preds_off <- predict.gam(off_gam_lcl, newdata = data.frame(lin_d2_TAST = new_data, treatment = "TAST OFF"), type = "response")
Ducl_preds_off <- predict.gam(off_gam_ucl, newdata = data.frame(lin_d2_TAST = new_data, treatment = "TAST OFF"), type = "response")

predictions_diff <- Dpredictions_off  - Dpredictions_on
lcl_preds_diff <- 
ucl_preds_difff <- 

# Specify the two treatment conditions
treatment_condition1 <- "OFF"
treatment_condition2 <- "ON"

# Sample 1000 rows from each treatment condition
sampled_data <- predictData %>%
  filter(treatment %in% c(treatment_condition1, treatment_condition2)) %>%
  group_by(treatment) %>%
  sample_n(size = 2000, replace = FALSE)

# Reset row names
sampled_data <- sampled_data %>%
  ungroup() %>%
  rownames_to_column()

# Print the first few rows of the sampled data
## Double plot
ggplot(sampled_data) +
  theme_classic()+
  stat_smooth(aes(lin_d2_TAST, mean), se =F, color = "black")+
  stat_smooth(aes(lin_d2_TAST, LCL), se =F, color = "red", size = 0.5, linetype = 4) +
  stat_smooth(aes(lin_d2_TAST, UCL), se =F, color = "red", size = 0.5,linetype = 4) +
  stat_smooth(aes(lin_d2_TAST, pred_diff), se =F, color = "blue", size = 0.5,linetype = 4)+
  labs(x = "Distance to TAST (m)", y = "Pr(Presence)") 

  
 pred_diff <- sampled_data$mean[sampled_data$treatment == "OFF"] - sampled_data$mean[sampled_data$treatment == "ON"]
