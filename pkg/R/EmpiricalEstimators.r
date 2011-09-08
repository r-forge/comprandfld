####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Email: simone.padoan@unibg.it.
### Institute: Department of Information Technology
### and Mathematical Methods, University of Bergamo
### File name: EmpiricalEstimators.r
### Description:
### This file contains a set of procedures in order
### to estimate the empircal covariance or the extreme
### dependence structures for a given dataset.
### Last change: 06/09/2011.
####################################################

### Procedures are in alphabetical order.

EVariogram <- function(data, coordx, coordy=NULL, coordt=NULL, cloud=FALSE, grid=FALSE,
                       gev=c(0,1,0), lonlat=FALSE, maxdist=NULL, maxtime=NULL,
                       numbins=NULL, replicates=1, type='variogram')
  {
    call <- match.call()
    corrmodel <- 'gauss'
    ### Check the parameters given in input:
    if(is.null(type))
      type <- 'variogram'
    # Checks if its a variogram or a madogram
    if(!is.null(type) & all(type!='variogram', type!='madogram', type!='Fmadogram'))
      stop('the admitted types are: variogram or madogram\n')
    # Set the type of model:
    if(type=='variogram') model <- 'Gaussian'
    else model <- 'ExtGauss'
    # Checks if its a spatial or spatial-temporal random field:
    if(!is.null(coordt))
      if(is.numeric(coordt))
        if(length(coordt)>1) corrmodel <- 'gneiting'
    # Checks the input:
    checkinput <- CheckInput(coordx, coordy, coordt, corrmodel, data, NULL, grid, 'None',
                             lonlat, "Frechet", maxdist, maxtime, model, 'Nelder-Mead',
                             replicates, NULL, NULL, NULL, 'WLeastSquare', FALSE, 'SubSamp',
                             FALSE, NULL, NULL)
    # Checks if there are errors in the input:
    if(!is.null(checkinput$error))
      stop(checkinput$error)
    ### START -- Specific checks of the Empirical Variogram:
    if(!is.null(cloud) & !is.logical(cloud))
      stop('insert a logical value (TRUE/FALSE) for the cloud parameter\n')
    if(type=='madogram' & (!is.numeric(gev) || is.null(gev) || !length(gev)==3))
      stop('insert a numeric vector with the three GEV parameters\n')
    if(!is.null(numbins) & !is.integer(numbins))
      if(numbins < 0)
        stop('insert a positive integer value for the number of bins\n')
    if(is.null(numbins))
      numbins <- 13
    ### END -- Specific checks of the Empirical Variogram

    ### Initialization parameters:
    initparam <- InitParam(coordx, coordy, coordt, corrmodel, data, NULL, grid, 'None',
                           lonlat, "Frechet", maxdist, maxtime, model, NULL, FALSE,
                           replicates, NULL, NULL, 'WLeastSquare', FALSE, 'SubSamp', FALSE, 1)
    # Checks if there are inconsistences:
    if(!is.null(initparam$error))
      stop(initparam$error)

    numvario <- numbins-1
    fname <- 'Binned_Variogram'
    if(cloud){
        numbins <- numvario <- initparam$numpairs
        fname <- 'Cloud_Variogram'}
    ### Estimation of the empirical spatial or spatial-temporal variogram:
    bins <- double(numbins) # spatial bins
    moments <- double(numvario) # vector of spatial moments
    lenbins <- integer(numvario) # vector of spatial bin sizes
    bint <- NULL
    lenbinst <- NULL
    lenbint <- NULL
    variogramst <- NULL
    variogramt <- NULL
    if(type=='madogram'){ # Trasform to GEV margins:
      if(cloud) fname <- 'Cloud_Madogram' else fname <- 'Binned_Madogram'
        if(gev[3]>=1)
          stop('the shape parameter can not be greater or equal to 1')
        initparam$data <- Dist2Dist(initparam$data, to='Gev',
                                    loc=rep(gev[1], initparam$numcoord),
                                    scale=rep(gev[2], initparam$numcoord),
                                    shape=rep(gev[3], initparam$numcoord))}
    if(type=='Fmadogram'){ # Transform to Uniform margins:
      if(cloud) fname <- 'Cloud_Madogram' else fname <- 'Binned_Madogram'
      initparam$data <- Dist2Dist(initparam$data, to='Uniform')}
    if(initparam$spacetime){
      numbint <- initparam$numtime-1 # number of temporal bins
      bint <- double(numbint)        # temporal bins
      momentt <- double(numbint)     # vector of temporal moments
      lenbint <- integer(numbint)    # vector of temporal bin sizes
      numbinst <- numvario*numbint   # number of spatial-temporal bins
      binst <- double(numbinst)      # spatial-temporal bins
      momentst <- double(numbinst)   # vector of spatial-temporal moments
      lenbinst <- integer(numbinst)  # vector of spatial-temporal bin sizes
      if(cloud) fname <- 'Cloud_Variogram_st' else fname <- 'Binned_Variogram_st'
      # Compute the spatial-temporal moments:
      .C(fname, bins, bint, as.double(initparam$data), lenbins, lenbinst, lenbint,
         moments, momentst, momentt, as.integer(numbins), as.integer(numbint),
         PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)
      indbin <- lenbins>0
      indbint <- lenbint>0
      indbinst <- lenbinst>0
      centers <- bins[1:numvario]+diff(bins)/2
      bins <- bins[indbin]
      bint <- bint[indbint]
      centers <- centers[indbin]
      moments <- moments[indbin]
      lenbins <- lenbins[indbin]
      momentt <- momentt[indbint]
      lenbint <- lenbint[indbint]
      momentst <- momentst[indbinst]
      lenbinst <- lenbinst[indbinst]
      # Computes the spatial marginal variogram:
      variograms <- moments/lenbins
      # Computes the temporal marginal variogram:
      variogramt <- momentt/lenbint
      # Computes the spatial-temporal variogram:
      variogramst <- momentst/lenbinst}
    else{# Computes the spatial moments
      .C(fname, bins, as.double(initparam$data), lenbins, moments, as.integer(numbins),
         PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)
      # Computes the spatial variogram:
      indbin <- lenbins>0
      bins <- bins[indbin]
      numbins <- length(bins)
      # check if cloud or binned variogram:
      if(cloud) centers <- bins else centers <- bins[1:(numbins-1)]+diff(bins)/2
      centers <- centers[indbin]
      moments <- moments[indbin]
      lenbins <- lenbins[indbin]
      variograms <- moments/lenbins}
    # Start --- compute the extremal coefficient
    extcoeff <- NULL
    if(type == 'madogram'){ # Check if its a madogram:
      if(gev[3] == 0)
        extcoeff <- exp(variograms / gev[2])
      else
        extcoeff <- Dist2Dist(gev[1] + variograms / gamma(1 - gev[3]))
      extcoeff[extcoeff>2] <- NA}
    if(type == 'Fmadogram'){ # Check if its a Fmadogram:
      extcoeff <- (1 + 2 * variograms) / (1 - 2 * variograms)
      extcoeff[extcoeff>2] <- NA}

    .C('DeleteGlobalVar', PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)

    EVariogram <- list(bins=bins,
                       bint=bint,
                       cloud=cloud,
                       centers=centers,
                       extcoeff=extcoeff,
                       lenbins=lenbins,
                       lenbinst=lenbinst,
                       lenbint=lenbint,
                       srange=initparam$srange,
                       variograms=variograms,
                       variogramst=variogramst,
                       variogramt=variogramt,
                       trange=initparam$trange,
                       type=type)

    structure(c(EVariogram, call = call), class = c("Variogram"))

  }

