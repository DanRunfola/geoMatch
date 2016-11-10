#' geoMatch
#'
#' geoMatch improves models using spatial data for the purposes of causal inference by selecting matched subsets of the original treated and control groups.
#' It provides an extension of the R package MatchIt (see Ho, Imai, King, and Stuart (2004)), enabling the use 
#' of spatial data frames as well as providing an adjustment factor to mitigate potential spillover between 
#' treated and control units (i.e., in cases where stable unit treatment value assumptions (SUTVA) may not be met
#' due to spatial spillovers).  geoMatch maintains the full functionality of MatchIt, including all 
#' matching strategies and caliper functions. Full documentation on MatchIt is available online at http://gking.harvard.edu/matchit, 
#' and help for specific commands is available through help.matchit.
#'
#'
#' @docType package
#' @name geoMatch
#' @export
#' @import sp
#' @import MatchIt
#' @details
#' geoMatch overcomes two challenges.  First, it seeks to first remove spillover from control outcome measures,
#' for example a case where a clinic may improve health outcomes in both the geographic neighborhood
#' it is located in, as well as nearby neighborhoods.  Failing to adjust for this spillover can result in
#' erroneous estimates of impact. Second, it allows for the use of spatial points data frames in conjunction
#' with the MatchIt matching framework.
#' geoMatch returns an adjusted version of the outcome variable for each control case as specified by the user, 
#' defined as Y*.  Y* can be interpreted as the estimated outcome if the spatial spillovers from treated
#' cases are netted out. This adjusted Y* can be used in the second stage model for more accurate estimations of
#' treatment effects.  Y* is calculated through a multiple step process.
#' First, a distance matrix (Dct) is constructed which provides the Euclidean distances between each treated and control case:
#' \preformatted{
#' [Dc1t1  Dc1t2   Dc1t3  .   Dc1tn]
#' [  .      .       .    .     .  ]
#' [  .      .       .    .     .  ]
#' [Dcnt1  Dcnt2   Dcnt3  .   Dcntn]
#' }
#' Second, two vector (Yt and Yc) are constructed which contain the known, observed outcomes in both Yt and Yc.
#' Third, a spherical distance decay function is fit for each T simultanesouly, in which the parameter Ut is 
#' solved for across all units T in the spatial function sf(distance, U):
#' \preformatted{
#' (1 - [[3/2] * (Dct / Ut) - [1/2] * (Dct/Ut)^3])
#' }
#' where we seek to optimize the absolute difference between the observed outcome at each control location (Yc) and the estimated
#' outcome (Yc_e) as a product of neighboring units Yt:
#' \preformatted{
#' Yc_e =  sf(Dct, Ut) * Yt 
#' minimize(abs(Yc_e-Yc))
#' }
#' The vector of distances Ut is used to estimate the 
#' adjusted Yc*, which - for each control case - provides an estimate of the outcome
#' if spillover is removed:
#' \preformatted{
#' Yc* = Yc - (sf[Dct, Ut]*Yt).
#' }
#' geoMatch then proceeds as usual with MatchIt, returning a matched set of spatial locations alongside Yc*.
#'@param ... The first parameters provided to geoMatch should be a traditional MatchIt specification.  Full documentation on MatchIt is available online at http://gking.harvard.edu/matchit. Dataframe must be a spatial points dataframe.
#'@param outcome.variable The name of the outcome variable that will be modeled to establish causal effect.  This must be an existing attribute in the spatial dataframe passed to geoMatch.
#'@param outcome.suffix Suffix for the returned column name with spillover-adjusted outcome data.
#'@return This function will return a MatchIt object, with a spatial data frame accesible in $spdf.
#'@example /demo/match.spatial.data.R
#'@source \url{http://geo.aiddata.org/}
#'@section Authors:
#'Dan Miller Runfola dan@danrunfola.com; Ariel BenYishay abenyishay@wm.edu; Seth Goodman sgoodman@aiddata.org
#'@section Citation:
#' \preformatted{
#' If you find this package useful, please cite:
#'Runfola, D., BenYishay, A., Goodman, S., 2016. geoMatch: An R Package for Spatial Propensity Score Matching. R package. http://geo.aiddata.org.
#'Daniel Ho, Kosuke Imai, Gary King, and Elizabeth Stuart (2007). Matching as Nonparametric Preprocessing for Reducing Model Dependence in Parametric Causal Inference. Political Analysis 15(3): 199-236. http://gking.harvard.edu/files/abs/matchp-abs.shtml
#'}

geoMatch <- function (..., outcome.variable,outcome.suffix="_adjusted", max.it = 10000, verbose=FALSE)
{
  spatial_match = -1
  a <- list(...)
  if (missing(outcome.variable) & isS4(a[['data']]))
  {
    warning("You are using a spatial dataframe, but did not assign a outcome.variable. This is required for spatial adjustments.  MatchIt will be performed ignoring spatial options.")
    spatial_match = FALSE
    if("data" %in% names(a))
    {
      #Pass the data name forward for a cleaner summary
      dName <- as.character(match.call()$data)
    }
    dfA <- as.data.frame(a[['data']]@data)
    a[['data']] <- dfA
  }
  
  if(!missing(outcome.variable) & !isS4(a[['data']]) & spatial_match==-1)
  {
    warning("You assigned an outcome variable, but your data is not a spatial data frame. A spatial data frame is required for spatial adjustments.  MatchIt will be performed ignoring spatial options.")
    spatial_match = FALSE  
    if("data" %in% names(a))
    {
      #Pass the data name forward for a cleaner summary
      dName <- as.character(match.call()$data)
    }
  }
  if(missing(outcome.variable) & !isS4(a[['data']]) & spatial_match==-1)
  {
    warning("To use the spatial functions of geoMatch, you must specify a spatial data frame and outcome.variable.  MatchIt is being run ignoring spatial options.")
    spatial_match = FALSE
    if("data" %in% names(a))
    {
      #Pass the data name forward for a cleaner summary
      dName <- as.character(match.call()$data)
    }
    
  }
  
  if(spatial_match == FALSE)
  {
    m.out <- do.call("matchit", a)
    m.out$call$data <- dName
    return(m.out)
  }
  
  if(!missing(outcome.variable) & isS4(a[['data']]) & spatial_match != FALSE)
  {
    x.coord <- coordinates(a[['data']])[,2]
    y.coord <- coordinates(a[['data']])[,1]
    
    if(max(x.coord) > 180.0 | min(x.coord) < -180.0 | max(y.coord) > 90.0 | min(y.coord) < -90.0)
    {
      warning("geoMatch currently only supports data projected with latitude and longitude information (WGS84).  Please reproject your data.")
      return(0)
    }
    
    if(class(a[['data']])[1] == "SpatialPointsDataFrame")
    {
      o_var <- outcome.variable
      geoMatch.Core(..., outcome.variable = o_var, m.it=max.it, verb = verbose) 
    }
    else
    {
      warning("Currently only spatial points data frames are supported.  Please convert your spatial data to a Spatial Points Data Frame.")
      return(0)
    }
  }
  else
  {
    warning("geoMatch encountered an error. Ensure that you have specified a valid spatial dataframe and outcome.variable.")
    return(0)
  }
  
}



