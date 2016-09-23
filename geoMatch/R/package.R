#' geoMatch
#'
#' A package designed to provide matched data for propensity-score based causal analyses using spatially-explicit data
#'
#' @docType package
#' @name geoMatch

geoMatch <- function (match,mod,dta)
{
  #First, calculate the OVB function (third order poly for now).
  
  #Second, for every Treated unit
    #use MatchIt as specified by the user to identify the best match
    #(or set of matches).
    #Take the difference in outcomes from each T to the matched
    #C(s). If multiple Cs are matched, use MatchIt weights.
  
      #Repeat this over multiple distance thresholds (for every D,
      #Remove all Cs which are within the D threshold of ANY T, not just the
      #One being analyzed).
      #Record the distance (or, average dist) of the best found match from the T unit.
  
    #After all distance thresholds have been run
    #Parameterize a distance-decay model (3rd order poly for now)
    #that summarizes the estimates spillover for that TE (Si)
  
  #Third, run a normal, full MatchIt model according to the users specified MatchIt.
  
  #Fourth, Calculate the Distances between all control cases included in the model and the nearest T.
  
  #Fifth, adjust each C's outcome measure according to the estimated spillover, given
  #the distance to the closest T (assuming non-additive spillover) as well as the OVB function.
  
  #Sixth, this function returns the new dataframe with the distance-adjusted outcomes,
  #which do not include spillover effects, as well as the MatchIt pairs.
  
  #Outside of this function, the user runs their model with the returned, adjusted, 
  #matched and transformed dataset.

  
  
}
