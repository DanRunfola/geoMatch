#' geoMatch
#'
#' A package designed to provide matched data for propensity-score based causal analyses using spatially-explicit data
#'
#' @docType package
#' @name geoMatch
require(MatchIt)
geoMatch <- function (..., outcome.variable)
{
  a <- list(...)
  if (missing(outcome.variable) & isS4(a[['data']]))
    {
      warning("You are using a spatial dataframe, but did not assign a outcome.variable.  \nThis is required for spatial adjustments.  \nMatchIt will be performed ignoring spatial options.")
      spatial_match = FALSE
      if("data" %in% names(a))
      {
        #Pass the data name forward for a cleaner summary
        dName <- as.character(match.call()$data)
      }
      dfA <- as.data.frame(a[['data']]@data)
      a[['data']] <- dfA
  }
  
  if(!missing(outcome.variable) & !isS4(a[['data']]))
    {
    warning("You assigned an outcome variable, but your data is not a spatial data frame.  \nA spatial data frame is required for spatial adjustments.  \nMatchIt will be performed ignoring spatial options.")
    spatial_match = FALSE  
    if("data" %in% names(a))
      {
        #Pass the data name forward for a cleaner summary
        dName <- as.character(match.call()$data)
      }
  }
  if(missing(outcome.variable) & !isS4(a[['data']]))
  {
    warning("To use the spatial functions of geoMatch, you must specify a spatial \ndata frame and outcome.variable.\nMatchIt is being run ignoring spatial options.")
    spatial_match == FALSE
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
}
  
  #First, calculate the OVB function (third order poly for now).
  
  #Second, for every Treated unit
    #use MatchIt as specified by the user to identify the best match
    #(or set of matches).
    #Record the matches found, and calculate the distances between Cs and Ts (matrix for all Ts)
  
  #Parameterize a model which best predicts Cs as a function of the distance to the matched Ts
  #Across all matrices.  For each individual T, this function defines it's spillover.
  
  #Third, run a normal, full MatchIt model according to the users specified MatchIt.
  
  #Fourth, Calculate the Distances between all control cases included in the model and all Ts.
  
  #Fifth, adjust each matched C's outcome measure (Y*) according to the estimated spillover, given
  #the distance to Ts as well as the OVB function.
  
  #Sixth, this function returns the new dataframe with the distance-adjusted outcomes,
  #which do not include spillover effects, as well as the MatchIt pairs.
  
  #Outside of this function, the user runs their model with the returned, adjusted, 
  #matched and transformed dataset.

