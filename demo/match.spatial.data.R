library(devtools)
install_github("itpir/geoMatch")
library(geoMatch)

###
### An Example Script for Obtaining Matched Data when you have
### Spatial information
###
data(lalonde)

##Traditional, non-spatial MatchIt
match.out1 <- matchit(treat ~ age + educ + black + hispan + 
                        nodegree + married + re74 + re75, 
                  method = "nearest", data = lalonde)

##Example model performed after matching, including both Control and Treatment groups
lm.out1 <- lm(re78 ~ treat + age + educ + black + hispan + 
                nodegree + married + re74 + re75 + distance, 
              data = match.data(match.out1))
summary(lm.out1)

##Simulate Latitude and Longtiude information for each point, 
##with enforced spatial correlation.
set.seed(500)
coords = cbind(runif(nrow(lalonde),37.1708,37.3708), 
               runif(nrow(lalonde),76.6069,76.8069))

##Create a spatial points data frame
spatial_lalonde <- SpatialPointsDataFrame(coords, lalonde)

##Matching and adjusting for spillover effects
##See ?geoMatch for more parameters specific to spatial data.
##See ?MatchIt for more options for matching methods.
match.out2 <- geoMatch(treat ~ age + educ + black + hispan + nodegree + 
                         married + re74 + re75, 
                      method = "nearest", 
                      caliper=0.25, 
                      data = spatial_lalonde, 
                      outcome.variable="re78", 
                      outcome.suffix="_adjusted")


##Example maps 
spplot(match.out2$spdf, z="matched", col.regions=c("red","green"), 
       main="Map of Matched Pairs")
spplot(match.out2$spdf, z="distance", main="Propensity Scores")
spplot(match.out2$spdf, z="est_spillovers", main="Estimated Spillovers")
spplot(match.out2$spdf, z="re78_adjusted", main="Adjusted Outcome")

#Percent of outcomes attributable to spillovers
match.out2$spdf@data["spill_percent"] <- 100 * 
  (match.out2$spdf@data["est_spillovers"] / match.out2$spdf@data["re78"])

spplot(match.out2$spdf[!is.infinite(match.out2$spdf@data$spill_percent) & 
                         match.out2$spdf@data$treat == 0,], 
       z="spill_percent", 
       main="% Outcome Attributable to Spillover",
       pretty=TRUE,
       cuts=5)

##Example model performed after spatial spillover adjustment, using matched data
lm.out2 <- lm(re78_adjusted ~ treat + age + educ + black + hispan + nodegree + 
                married + re74 + re75 + distance, 
              data = match.data(match.out2))
summary(lm.out2)

##Example model with spatial lag after spatial spillover adjustment
library(spdep)
coords <- coordinates(match.out2$spdf[match.out2$spdf@data$matched == 1,])
#Weights matrix based on 15 nearest neighbors
k1 <- knn2nb(knearneigh(coords, k=15))
all.linked <- max(unlist(nbdists(k1, coords)))
neighbor.list <- dnearneigh(coords, 0, all.linked)
sl.out3 <- lagsarlm(re78_adjusted ~ treat + age + educ + black + hispan + 
                      nodegree + married + re74 + re75 + distance, 
                    data=match.out2$spdf[match.out2$spdf@data$matched == 1,],
                        nb2listw(neighbor.list, style="W"), method="eigen", quiet=FALSE)
summary(sl.out3)
