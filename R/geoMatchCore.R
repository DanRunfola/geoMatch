#Core functionality for the geoMatch MatchIt wrapper.
geoMatch.Core <- function (..., outcome.variable,outcome.suffix="_adjusted"){
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
    S <- rep(Ut[length(Ut)],length(Ut)-1)
    D <- Ut[1:(length(Ut)-1)]
    Yc.spill.est.genA = S * ((3/2) * (Dct / D) - (1/2) * (Dct/D)^3)
    Yc.spill.est.genA[Yc.spill.est.genA < 0.0] <- 0
    Yc.spill.est.genB <- sweep(Yc.spill.est.genA,MARGIN=2,Yt[[1]],'*')
    Yc.spill.est <- rowSums(Yc.spill.est.genB)
    Yc.err = sum(abs(Yc - Yc.spill.est))
    return(Yc.err)
  }
  
  sf <- function(...)
  {
    S <- rep(Ut[length(Ut)],length(Ut)-1)
    D <- Ut[1:(length(Ut)-1)]
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
  Ut <- c(runif(nrow(Yt),(min(Dct)+.00001),max(Dct)),runif(1))
  #print(Ut)
  m_init <- round(max(Dct)*4,0)
  Ut.optim <- 
    optimx(par = Ut, 
        fn=sf.opt, 
        gr=NULL,
        hess=NULL,
        method = "spg",
        lower = 0.01,
        upper=m_init,
        itnmax=100000,
        #hessian=FALSE,
        control=list(trace=0),
        Dct)

  if(Ut.optim$convcode != 0)
  {
    warning("No optimal spatial decay function was found, which can indicate a lack of spatial autocorrelation or a highly complex system. Consider using the unadjusted estimate.")
    warning(Ut.optim$message)
  }
  

    Ut <- Ut.optim[1:length(Ut)]


  #Calculate adjusted Yc*, which - for each C - removes spatial spillover.
  #Yc* = Yc - (sf[Dct, Ut]*Yt) [Note: Yt multiplier is applied in the function
  #to make this code easier to read].
    
  spillovers_c <- sf(Dct, Ut)
  Yc.star <- Yc - spillovers_c
  

  
  #Update the original outcome values to remove spillover.
  outcome.variable.adjusted <- paste(outcome.variable,outcome.suffix,sep="")
  a[['data']]@data[outcome.variable.adjusted] <- a[['data']]@data[outcome.variable]
  a[['data']]@data[a[['data']]@data[t.name]==0,][outcome.variable.adjusted] <- Yc.star
  a[['data']]@data["est_spillovers"] <- 0
  a[['data']]@data[a[['data']]@data[t.name]==0,]["est_spillovers"] <- spillovers_c
  
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
  m.res$optim <- (Ut.optim$convcode == 0)
  
  return(m.res)
  
  #These adjusted values can then be used in a traditional matchIt.
  
  

}