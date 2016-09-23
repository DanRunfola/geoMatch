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

##Simulate Latitude and Longtiude information for each point
set.seed(424)
coords = cbind(runif(nrow(lalonde),37.1708,37.3708), runif(nrow(lalonde),76.6069,76.8069))

##Create a spatial points data frame
spdf_LL <- SpatialPointsDataFrame(coords, lalonde)

##Traditional, non-spatially weighted matching
m.out1 <- matchit(treat ~ re74 + re75 + age + educ, data=lalonde,
                  method="nearest", distance="logit", caliper=0.25)

