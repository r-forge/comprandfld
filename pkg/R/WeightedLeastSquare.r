####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Email: simone.padoan@epfl.ch.
### Institute: EPFL.
### File name: WeightedLeastSquare.r
### Description:
### This file contains a set of procedures in order
### to estimate the parameters of some covariance 
### function models for a given dataset.
### Last change: 28/11/2010.
####################################################


### Procedures are in alphabetical order.

print.WLS <- function(x, digits = max(3, getOption("digits") - 3), ...)
  {
    dimdata <- dim(x$data)
    
    if(x$grid)
      numdata <- dimdata[3]
    else
      numdata <- dimdata[1]
    
    numcoord <- nrow(x$coord)
    numparam <- length(x$param)

    if(x$weighted)
      method <- 'Weighted Least Square'
    else
      method <- method <- 'Least Square'

    cat('\n##############################################################')
    cat('\nResults:', method,'Fitting of Random Fields.\n')
    cat('\nCovariance model:', x$corrmodel, '\n')
    cat('Number of coordinates:', numcoord, '\n')
    cat('Number of observations per location:', numdata, '\n')
    cat('Number of bins', length(x$bins),'\n')
    cat('\nEmpirical variogram:\n')
    print.default(x$variogram, digits = digits, print.gap = 1, quote = FALSE)

    cat('\nEstimated parameters:\n')
    print.default(x$param, digits = digits, print.gap = 2,
                  quote = FALSE)

    cat('\n##############################################################\n')    
    invisible(x)
  }

WlsInit <- function(coordx, coordy, corrmodel, data, fixed, grid, likelihood,
                    lonlat, model, parscale, paramrange, replicates, start, type,
                    vartype, weighted)
  {
    
    ### Initialization parameters:
    
    initparam <- InitParam(coordx, coordy, corrmodel, data, fixed, grid,
                           likelihood, lonlat, model, parscale, paramrange,
                           replicates, start, 'WLeastSquare', vartype, weighted)

    if(!is.null(initparam$error))
      stop(initparam$error)

    initparam$type <- CheckType(type)

    if(initparam$numstart == initparam$numparam)
      {
        if(type == 'Standard' || type == 'Pairwise')
          {
            if(!any(names(fixed)=='mean'))
              {
                initparam$param <- c(initparam$fixed['mean'], initparam$param)
                initparam$namesparam <- sort(names(initparam$param))
                initparam$param <- initparam$param[initparam$namesparam]
                initparam$numparam <- length(initparam$param)
                initparam$flagnuis['mean'] <- 1
                initparam$numfixed <- initparam$numfixed - 1
                if(is.null(fixed))
                  initparam$fixed <- NULL
                else
                  initparam$fixed <- initparam$fixed[!names(initparam$fixed)=='mean']
              }
          }

        return(initparam)
      }
    if(initparam$numstart > 0)
      {
        for(i in 1 : initparam$numstart)
          {
            initparam$param <- initparam$param[!names(initparam$param)==initparam$namesstart[i]]
            initparam$fixed <- c(initparam$fixed, initparam$start[initparam$namesstart[i]])
          }
      }

    ### Define the object function:
    WLsquare <- function(bins, corrmodel, fixed, fun, lenbins, moments, numbins, param)
      {
        result <- -1.0e15
        if(!CheckParamRange(param))
          return(result)

        result <- double(1)
        param <- c(param, fixed)
        namesparam <- sort(names(param))
        param <- param[namesparam]
        nuisance <- c(param['nugget'], param['sill'])
        corrparam <- param[-match(c('mean', 'nugget', 'sill'), namesparam)]

        result <- .C(fun, as.double(bins), as.integer(corrmodel),
                     as.double(corrparam), as.integer(numbins), as.double(moments),
                     as.double(lenbins), as.double(nuisance), res=double(1),
                     PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)$res
        return(result)
      }

    fname <- NULL

    ### Estimation of the empirical variogram:
    numbins <- as.integer(13)
    maxdist <- double(1)

    bins <- double(numbins)
    moments <- double(numbins - 1)
    lenbins <- integer(numbins - 1)

    .C('Empiric_Variogram', bins, as.integer(FALSE), as.double(initparam$data), lenbins,
       maxdist, moments, as.integer(initparam$numdata), as.integer(initparam$numcoord),
       numbins, as.integer(1), PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)

    ###### ---------- START model fitting ----------######
    
    if(initparam$model == 1) # Gaussian case:
      fname <- 'LeastSquare_G'

    fitted <- optim(initparam$param, WLsquare, bins=bins, corrmodel=initparam$codecorrmodel,
                    fixed=initparam$fixed, fun=fname, lenbins=lenbins, moments=moments,
                    numbins=numbins, control=list(fnscale=-1, reltol=1e-14, maxit=1e8),
                    hessian=FALSE)

    ###### ---------- END model fitting ----------######

    if(initparam$numstart > 0)
      {
        for(i in 1 : initparam$numstart)
          {
            initparam$param <- c(initparam$param, initparam$start[initparam$namesstart[i]])
            initparam$fixed <- initparam$fixed[!initparam$namesfixed==initparam$namesstart[i]]
          }
        initparam$param <- initparam$param[initparam$namesparam]
      }

    if(fitted$convergence == 0)
      {
        if(type == 'Standard' || type == 'Pairwise')
          {
            if(!any(names(fixed)=='mean'))
              {
                initparam$param <- c(initparam$fixed['mean'], initparam$param)
                initparam$namesparam <- sort(names(initparam$param))
                initparam$param <- initparam$param[initparam$namesparam]
                initparam$numparam <- length(initparam$param)
                initparam$flagnuis['mean'] <- 1
                initparam$numfixed <- initparam$numfixed - 1
                if(is.null(fixed))
                  initparam$fixed <- NULL
                else
                  initparam$fixed <- initparam$fixed[!names(initparam$fixed)=='mean']
              }
          }
        initparam$param[names(fitted$par)] <- fitted$par
      }
  
   return(initparam)
  }

  
WLeastSquare <- function(coordx, coordy, corrmodel, data, fixed=NULL, grid=FALSE,
                         lonlat=FALSE, maxdist=NULL, model='Gaussian',
                         optimizer='Nelder-Mead', numbins=NULL, replicates=FALSE,
                         start=NULL, weighted=FALSE)
  {
    
    call <- match.call()
    
    ### Check the parameters given in input:

    checkinput <- CheckInput(coordx, coordy, corrmodel, data, fixed, grid, 'None',
                             lonlat, model, optimizer, replicates, start, 'WLeastSquare',
                             FALSE, 'SubSamp', weighted, NULL, NULL)

    if(!is.null(checkinput$error))
      stop(checkinput$error)

    if(!is.null(numbins) & !is.integer(numbins))
      if(numbins < 0)
        stop('insert a positive integer value for the number of bins\n')

    if(is.null(numbins))
      numbins <- 13

    if(!is.null(maxdist) & !is.numeric(maxdist))
      if(maxdist < 0)
        stop('insert a positive numeric value for maximum distance\n')

    if(is.null(maxdist))
      maxdist <- 0

    ### Define the object function:

    WLsquare <- function(bins, corrmodel, fixed, fun, lenbins, moments, numbins, param)
      {
        result <- -1.0e15
        if(!CheckParamRange(param))
          return(result)

        result <- double(1)
        param <- c(param, fixed)
        namesparam <- sort(names(param))
        param <- param[namesparam]
        nuisance <- c(param['nugget'], param['sill'])
        corrparam <- param[-match(c('mean', 'nugget', 'sill'), namesparam)]

        result <- .C(fun, as.double(bins), as.integer(corrmodel),
                     as.double(corrparam), as.integer(numbins), as.double(moments),
                     as.double(lenbins), as.double(nuisance), res=double(1),
                     PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)$res
        return(result)
      }
 
    ### Initialization global variables:
     
    WLeastSquare <- NULL
    fname <- NULL

    ### Initialization parameters:
    parscale <- NULL
    initparam <- InitParam(coordx, coordy, corrmodel, data, fixed, grid, 'None',
                           lonlat, model, parscale, optimizer=='L-BFGS-B',
                           replicates, start, 'WLeastSquare', 'SubSamp', FALSE)
  
    if(!is.null(initparam$error))
      stop(initparam$error)

    ### Estimation of the empirical variogram:

    bins <- double(numbins)
    moments <- double(numbins - 1)
    lenbins <- integer(numbins - 1)

    if(initparam$model > 1)
      initparam$data <- Dist2Dist(initparam$data, to='Uniform')

    .C('Empiric_Variogram', bins, as.integer(FALSE), as.double(initparam$data),
       lenbins, as.double(maxdist), moments, as.integer(initparam$numdata),
       as.integer(initparam$numcoord), as.integer(numbins), as.integer(initparam$model),
       PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)

    ###### ---------- START model fitting ----------######
    
    if(initparam$model == 1) # Gaussian case:
      {
        if(weighted)
          fname <- 'WLeastSquare_G'
        else
          fname <- 'LeastSquare_G'
      }

    if(initparam$model > 1) # Max-stable case:
      {
        variogram <- moments / lenbins
        extcoeff <- (1 + 2 * variogram) / (1 - 2 * variogram)
        moments <- extcoeff
      }

    if(optimizer=='L-BFGS-B')
      fitted <- optim(initparam$param, WLsquare, bins=bins, corrmodel=initparam$codecorrmodel,
                      fixed=initparam$fixed, fun=fname, lenbins=lenbins, method=optimizer,
                      moments=moments, numbins=numbins, control=list(fnscale=-1, factr=1,
                      pgtol=1e-14, maxit = 1e8), lower=initparam$lower, upper=initparam$upper,
                      hessian=FALSE)
    else
      fitted <- optim(initparam$param, WLsquare, bins=bins, corrmodel=initparam$codecorrmodel,
                      fixed=initparam$fixed, fun=fname, lenbins=lenbins, method=optimizer,
                      moments=moments, numbins=numbins, control=list(fnscale=-1, reltol=1e-14,
                      maxit=1e8), hessian=FALSE)

    ###### ---------- END model fitting ----------######

    .C('DelDistances', PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)

    WLeastSquare <- list(bins=bins,
                         coord = initparam$coord,
                         convergence = fitted$convergence,
                         corrmodel = corrmodel,
                         data = initparam$data,
                         fixed = initparam$fixed,
                         grid = grid,
                         iterations = fitted$counts,
                         message = fitted$message,
                         param = fitted$par,
                         variogram = moments / lenbins,
                         weighted = weighted)

    structure(c(WLeastSquare, call = call), class = c("WLS"))
  }
