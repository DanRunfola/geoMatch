#Core functionality for the geoMatch MatchIt wrapper.
geoMatch.Core <- function (..., outcome.variable,outcome.suffix="_adjusted", m.it = 10000, v = FALSE){
  a <- list(...)
  #Remove spillover from all C outcomes, then matches.
  #Mitigates OVB through PSM matching, but does not take
  #advantage of OVB back-extrapolation.
  #To correct for potential OVB in spatially-correlated systems,
  #see geoOVB.
  
  #For each matched unit, this outputs a new variable (Y*), which is
  #the outcome with spatial spillovers from T netted out.
  #This should be used in your second stage model.
  
  #Dct - Distance between each control case and treatment case (matrix)
  #[Dc1t1  Dc1t2   Dc1t3  .   Dc1tn]
  #[  .      .       .    .     .  ]
  #[  .      .       .    .     .  ]
  #[Dcnt1  Dcnt2   Dcnt3  .   Dcntn]


  
  t.name <- strsplit(as.character(a[[1]]),"~")[[2]]
  t.coords <- a[['data']][as.vector(a[['data']]@data[t.name]==1),]
  c.coords <- a[['data']][as.vector(a[['data']]@data[t.name]==0),]
  Dct <- matrix(, nrow = length(c.coords), ncol = length(t.coords))
  rownames(Dct) <- rownames(c.coords@data)
  colnames(Dct) <- rownames(t.coords@data)
  for(i in 1:(length(c.coords)))
  {
    ct <- spDistsN1(t.coords, c.coords[i,], longlat = TRUE)
    Dct[i,] <- ct
  }
  
  #Yt - Vector. Known outcome at each T.
  Yt.gen <- a[['data']][as.vector(a[['data']]@data[t.name]==1),]
  Yt <- Yt.gen@data[outcome.variable]
  
  #Yc - Vector. Known outcome at each C.
  Yc.gen <- a[['data']][as.vector(a[['data']]@data[t.name]==0),]
  Yc <- Yc.gen@data[outcome.variable]
  
  #e - residual variance for units c and/or t (depending on subscript)
  #sf(Dct, Ut) - A spatial decay function, variably defined.
  #This iteration of the script assumes spherical:
  #(1 - [[3/2] * (Dct / Ut) - [1/2] * (Dct/Ut)^3])
  
  #Ut - Vector. Parameter "U" in spherical 
  #distance decay relationships, solved for 
  #each T simultaneously
  
  #Solve for Ut in Yc =  sf(Dct, Ut) * Yt + ec
  

  sf.opt <- function(Ut, ...)
  {

    S <- tail(Ut, length(Ut)/2) 
    D <- Ut[1:length(S)] * ltyp.scale
    Yc.spill.est.genA = S * ((3/2) * (Dct / D) - (1/2) * (Dct/D)^3)
    Yc.spill.est.genA[Yc.spill.est.genA < 0.0] <- 0
    Yc.spill.est.genB <- sweep(Yc.spill.est.genA,MARGIN=2,Yt[[1]],'*')
    Yc.spill.est <- rowSums(Yc.spill.est.genB)
    #print(mean(abs(Yc - Yc.spill.est)))
    print(abs(Yc - Yc.spill.est)[1])
    Yc.err = mean(as.numeric(abs(Yc - Yc.spill.est)[[1]]))
    return(Yc.err)
  }
  
  sf <- function(...)
  {
    S <- tail(Ut, length(Ut)/2) 
    D <- Ut[1:length(S)] * ltyp.scale
    Yc.spill.est.genA = S * ((3/2) * (Dct / D) - (1/2) * (Dct/D)^3)
    Yc.spill.est.genA[Yc.spill.est.genA < 0.0] <- 0
    Yc.spill.est.genB <- sweep(Yc.spill.est.genA,MARGIN=2,Yt[[1]],'*')
    Yc.spill.est <- rowSums(Yc.spill.est.genB)
    return(Yc.spill.est)
  }
  

  #Optimize Ut vector, limiting distance threshold to be
  #greater than 0, and less than the circumfrence of the earth 
  #(Approx. 40,100 km)
  #Use random starting points between the minimum and maximum observed distances
  #between C and T.
  Ut <- runif(nrow(Yt)*2,.1,1)
  high_init <- runif(nrow(Yt)*2,4.0,4.0) #Maximum decay relative to largest distance observed.
  low_init <- runif(nrow(Yt)*2,.00001,.00001) 
  
  #Parameter scaling
  ltyp.scale <- max(Dct)

  if(v == TRUE)
  {
    q <- FALSE
    t <- TRUE
  }
  else
  {
    q <- TRUE
    t <- FALSE
  }
  
  Ut.optim <- 
    spg(par = Ut, 
        fn=sf.opt, 
        gr=NULL,
        lower = .00000000001,
        ftol = .0000001,
        upper= 10,
        control=list(trace=t, maxit = m.it, checkGrad = FALSE, M=15),
        Dct = Dct,
        quiet = q,
        ltyp.scale = ltyp.scale)
  
  if(Ut.optim$convergence != 0)
  {
    warning("No optimal spatial decay function was found, which can indicate a lack of spatial autocorrelation or a highly complex system. Consider using the unadjusted estimate.")
    warning(Ut.optim$message)
  }
  

    Ut <- Ut.optim$par
    
  

  #Calculate adjusted Yc*, which - for each C - removes spatial spillover.
  #Yc* = Yc - (sf[Dct, Ut]*Yt) [Note: Yt multiplier is applied in the function
  #to make this code easier to read].
  spillovers_c <- sf(Ut, Dct, p.scale)
  
  #Adjsut spillovers to account for spillover that could be explained by
  #covariate information.
  #This is a very conservative approach to overcoming potential bias due to
  #y = 1*B + 0 * C, where B is highly spatially autocorrelated.
  #Further, if omitted variables are correlated with ancillary data,
  #this additionally mitigates potential OVB (when taken in conjunction
  #with the spillover adjustments).
  
  #A predictive model is fit to account for this 
  a[['data']]@data["est_max_spillovers"] <- 0
  a[['data']]@data[a[['data']]@data[t.name]==0,]["est_max_spillovers"] <- spillovers_c
  x.var.str <- strsplit(as.character(a[[1]]), "~")[3]
  pred.str <- paste("est_max_spillovers", "~", x.var.str)

  prediction.results <- lm(pred.str, data=a[['data']]@data)
  
  #Remaining unexplained spillover potential.
  #Defined as the total difference between predicted values and estimated spillover values.
  #More complex metric of variance remaining should be considered in the future.
  #Capped at 100% and 0%.
  
  #absolute deviation vector
  abs.dev.vec <- abs(a[['data']]@data[a[['data']]@data[t.name]==0,]["est_max_spillovers"] - predict(prediction.results, newdata=a[['data']]))
  
  #Deviance as a percentage of total spillover
  #A "1" (or greater) would indicate no spillover was explained by the covariates.
  #A "0" would indicate all spillover was explained by the covariates.
  prop.exp <- abs.dev.vec / a[['data']]@data[a[['data']]@data[t.name]==0,]["est_max_spillovers"]
  prop.exp[prop.exp > 1] <- 1
  
  #Covariate-weighted adjustement
  Yc.star <- Yc - (spillovers_c * prop.exp)
  

  
  #Update the original outcome values to remove spillover.
  outcome.variable.adjusted <- paste(outcome.variable,outcome.suffix,sep="")
  a[['data']]@data[outcome.variable.adjusted] <- a[['data']]@data[outcome.variable]
  a[['data']]@data[a[['data']]@data[t.name]==0,][outcome.variable.adjusted] <- Yc.star
  
  #Object with pre- and post statistics.
  pre.post.spatial <- 
    c(summary(a[['data']]@data[a[['data']]@data[t.name]==0,][outcome.variable]),
      summary(a[['data']]@data[a[['data']]@data[t.name]==0,][outcome.variable.adjusted]))
  
  #Run matchIt as usual for the user.
  spdfA <- a[['data']]
  dfA <- as.data.frame(a[['data']]@data)
  a[['data']] <- dfA
  m.res <- do.call("matchit", a)

  spdfA@data$weights <- m.res$weights
  spdfA@data$matched <- m.res$weights>0
  spdfA@data$distance <- m.res$distance
  
  m.res$spdf <- spdfA
  m.res$optim <- (Ut.optim$convergence == 0)
  
  return(m.res)
  
  #These adjusted values can then be used in a traditional matchIt.
  
  

}