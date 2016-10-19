#' geoMatch
#'
#' A package designed to provide matched data for propensity-score based causal analyses using spatially-explicit data
#'
#' @docType package
#' @name geoMatch
#' @export
#' @import sp
#' @import MatchIt
require(MatchIt)
require(sp)

geoMatch <- function (..., outcome.variable)
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
    
    if(class(spatial_lalonde)[1] == "SpatialPointsDataFrame")
    {
      o_var <- outcome.variable
      geoMatch.Core(..., outcome.variable = o_var)
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
  


