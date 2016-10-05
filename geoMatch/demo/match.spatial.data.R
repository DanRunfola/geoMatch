library(devtools)
library(sp)
library(MatchIt)
detach("package:geoMatch", unload=TRUE)
load_all("~/Desktop/Github/geo.match/geoMatch/R")

###
### An Example Script for Obtaining Matched Data when you have
### Spatial information
###
data(lalonde)

##Traditional, non-spatial MatchIt
match.out1 <- matchit(treat ~ age + educ + black + hispan + nodegree + married + re74 + re75, 
                  method = "nearest", data = lalonde)

##Example model performed after matching, including both Control and Treatment groups
lm.out1 <- lm(re78 ~ treat + age + educ + black + hispan + nodegree + married + re74 + re75 +
                  distance, data = match.data(match.out1))

summary(lm.out1)

##Simulate Latitude and Longtiude information for each point
set.seed(424)
coords = cbind(runif(nrow(lalonde),37.1708,37.3708), runif(nrow(lalonde),76.6069,76.8069))

##Create a spatial points data frame
spatial_lalonde <- SpatialPointsDataFrame(coords, lalonde)

##Optionally plot the new spatial data
spplot(spatial_lalonde, z=c("age"))

##Spatially-weighted 
match.out2 <- geoMatch(treat ~ age + educ + black + hispan + nodegree + married + re74 + re75, 
                      method = "nearest", data = lalonde)#, outcome.variable="re78_adjusted")

##Example model performed after spatial adjustment, including both Control and Treatment groups
#lm.out2 <- lm(re78_adjusted ~ treat + age + educ + black + hispan + nodegree + married + re74 + re75 +
#                distance, data = match.data(match.out2))

#summary(lm.out2)