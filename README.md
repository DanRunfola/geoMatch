
# geoMatch

> geoMatch

A package designed to provide matched data for propensity-score based causal analyses using spatially-explicit data.  Fully compatible with matching strategies available in the MatchIt package, geoMatch provides adjusted outcome measures which mitigate bias associated with spillover from areas which received interventions to proximate areas that did not.

## Installation
```r
library(devtools)
install.github("itpir/geoMatch")
```
geoMatch imports both the sp and MatchIt packages.

## Usage

```r
#A full demo for usage can be found in /demo/match.spatial.data.R
library(geoMatch)
##Simulate Latitude and Longtiude information for each point, 
##with enforced spatial correlation.
library(sp)
set.seed(500)
coords = cbind(runif(nrow(lalonde),37.1708,37.3708), 
               runif(nrow(lalonde),76.6069,76.8069))

##Create a spatial points data frame
spatial_lalonde <- SpatialPointsDataFrame(coords, lalonde)

match.out2 <- geoMatch(treat ~ age + educ + black + hispan + nodegree + 
                         married + re74 + re75, 
                      method = "nearest", 
                      caliper=0.25, 
                      data = spatial_lalonde, 
                      outcome.variable="re78", 
                      outcome.suffix="_adjusted")
lm.out2 <- lm(re78_adjusted ~ treat + age + educ + black + hispan + nodegree + 
                married + re74 + re75 + distance, 
              data = match.data(match.out2))
summary(lm.out2)
```

## License

Copyright 2016 Daniel Miller Runfola (dan@danrunfola.com)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

