####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Email: simone.padoan@eofl.ch.
### Institute: EPFL.
### File name: Likelihood.r
### Description:
### This file contains a set of procedures
### for maximum likelihood fitting of
### random fields.
### Last change: 18/06/2010.
####################################################


### Procedures are in alphabetical order.

print.WLS <- function(x, digits = max(3, getOption("digits") - 4), ...)
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

    cat('\n##############################################################')    
    invisible(x)
  }


Wls <- function(bins, corrmodel, fixed, lenbins, moments, numbins, param, weighted)
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

    result <- .C('Wls', as.double(bins), as.integer(corrmodel), as.double(corrparam),
                 as.integer(numbins), as.double(moments), as.double(lenbins),
                 as.double(nuisance), as.integer(weighted), res=double(1),
                 PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)$res

    return(result)
  }

WlsInit <- function(coordx, coordy, corrmodel, data, fixed, grid, likelihood,
                    lonlat, model, parscale, paramrange, start, time, type)
  {
    
    ### Initialization parameters:
    
    initparam <- InitParam(coordx, coordy, corrmodel, data, fixed, grid,
                           likelihood, lonlat, model, parscale, paramrange,
                           start, time, 'WLeastSquare')

    if(!is.null(initparam$error))
      stop(initparam$error)

    initparam$type <- CheckType(type)

    ### Estimation of empirical variogram:

    numbins <- 13
    maxdist <- double(1)

    bins <- double(numbins)
    moments <- double(numbins - 1)
    lenbins <- double(numbins - 1)

    .C('Empiric_Variogram', as.double(bins), as.double(initparam$coord[,1]), as.double(initparam$coord[,2]),
       as.double(initparam$data), as.double(initparam$lags), as.double(lenbins), as.double(maxdist),
       as.double(moments), as.integer(initparam$numpairs), as.integer(initparam$numcoord), as.integer(numbins),
       PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)


    ### Model fitting:
    fitted <- optim(initparam$param, Wls, bins=bins, corrmodel=initparam$corrmodel,
                    fixed=initparam$fixed, lenbins=lenbins, moments=moments, numbins=numbins,
                    weighted=FALSE, control=list(fnscale=-1, reltol=1e-14, maxit=1e8),
                    hessian=FALSE)

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
  
    if(is.list(start))
      {
        namesstart <- names(start)
        for(i in 1 : length(start))
          initparam$param[namesstart[i]] <- as.numeric(start[namesstart[i]])
      }

   return(initparam)
  }

  
WLeastSquare <- function(coordx, coordy, corrmodel, data, fixed=NULL, grid=FALSE,
                         lonlat=FALSE, maxdist=NULL, optimizer='Nelder-Mead',
                         numbins=NULL, start=NULL, time=FALSE, weighted=FALSE)
  {
    
    call <- match.call()
    
    ### Check the parameters given in input:

    checkinput <- CheckInput(coordx, coordy, corrmodel, data, fixed, grid,
                             'None', lonlat, 'None', optimizer, start, FALSE,
                             time, 'WLeastSquare', weighted, NULL)

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
 
    ### Initialization global variables:
     
    WLeastSquare <- NULL

    ### Initialization parameters:
    parscale <- NULL
    initparam <- WlsInit(coordx, coordy, corrmodel, data, fixed, grid,
                         'None', lonlat, 'None', parscale, optimizer=='L-BFGS-B',
                         start, time, 'WLeastSquare')
  
    if(!is.null(initparam$error))
      stop(initparam$error)

    ### Estimation of empirical variogram:

    bins <- double(numbins)
    moments <- double(numbins - 1)
    lenbins <- double(numbins - 1)

    .C('Empiric_Variogram', bins, as.double(initparam$coord[,1]), as.double(initparam$coord[,2]),
       as.double(initparam$data), as.double(initparam$lags), lenbins, as.double(maxdist), moments,
       as.integer(initparam$numpairs), as.integer(initparam$numcoord), as.integer(numbins),
       PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)

    ### Model fitting:(bins, corrmodel, fixed, lenbins, moments, numbins, param, weighted)
    if(optimizer=='L-BFGS-B')
      fitted <- optim(initparam$param, Wls, bins=bins, corrmodel=initparam$codecorrmodel,
                      fixed=initparam$fixed, lenbins=lenbins, method=optimizer, moments=moments,
                      numbins=numbins, weighted=weighted, control=list(fnscale=-1, factr=1, pgtol=1e-14,
                      maxit = 1e8), lower=initparam$lower, upper=initparam$upper, hessian=FALSE)
    else
      fitted <- optim(initparam$param, Wls, bins=bins, corrmodel=initparam$codecorrmodel,
                          fixed=initparam$fixed, lenbins=lenbins, method=optimizer,
                          moments=moments, numbins=numbins, weighted=weighted,
                          control=list(fnscale=-1, reltol=1e-14, maxit=1e8), hessian=FALSE)


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
