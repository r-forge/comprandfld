####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Email: simone.padoan@eofl.ch.
### Institute: EPFL.
### File name: Fitting.r
### Description:
### This file contains a set of procedures
### for maximum composite-likelihood fitting of
### random fields.
### Last change: 04/06/2010.
####################################################


### Procedures are in alphabetical order.


### Fitting procedure: 

FitComposite <- function(coordx, coordy=NULL, corrmodel, data, fixed=NULL, grid=FALSE, hessian=FALSE,
                         likelihood='Marginal', lonlat=FALSE, model='Gaussian', optimizer='Nelder-Mead',
                         start=NULL, varest=FALSE, time=FALSE, type='Pairwise', weighted=FALSE, weights=NULL)
{
    call <- match.call()
    
    ### Check the parameters given in input:
    
    checkinput <- CheckInput(coordx, coordy, corrmodel, data, fixed, grid,
                             likelihood, lonlat, model, optimizer, start,
                             varest, time, type, weighted, weights)

    if(!is.null(checkinput$error))
      stop(checkinput$error)

    ### Initialization global variables:
     
    FitComposite <- NULL
    clic <- parscale <- varcov <- stderr <- NULL
  
    ### Initialization parameters:

    initparam <- WlsInit(coordx, coordy, corrmodel, data, fixed, grid, likelihood,
                         lonlat, model, parscale, optimizer=='L-BFGS-B', start, time, type)
    
    if(!is.null(initparam$error))
      stop(initparam$error)

    ### Model fitting section

    # Full likelihood:
    
    if(initparam$likelihood == 2)
      {
        hessian <- varest
        # Fitting by log-likelihood maximization:
        fitted <- OptimLik(initparam$corrmodel, initparam$data, initparam$fixed, grid, hessian,
                           initparam$lags, initparam$lower, initparam$model, initparam$namescorr,
                           initparam$namesnuis, initparam$namesparam, initparam$numcoord,
                           initparam$numdata, initparam$numpairs, optimizer, initparam$param,
                           varest, initparam$type, initparam$upper) 
      }
    
    # Composite likelihood:

    if(initparam$likelihood == 3 || initparam$likelihood == 1)
      {
        fitted <- OptimCompLik(initparam$coord[,1], initparam$coord[,2], initparam$corrmodel, initparam$data,
                               initparam$flagcorr, initparam$flagnuis, initparam$fixed, grid, FALSE, initparam$lags,
                               initparam$likelihood, initparam$lower, initparam$model, initparam$namescorr,
                               initparam$namesnuis, initparam$namesparam, initparam$numcoord, initparam$numdata,
                               initparam$numparam, initparam$numparamcorr, optimizer, initparam$param, varest,
                               initparam$type, initparam$upper, weighted)
      }

    ### Set the output object:
     
    FitComposite <- list(clic = fitted$clic,
                         coord = initparam$coord,
                         convergence = fitted$convergence,
                         corrmodel = corrmodel,
                         data = initparam$data,
                         fixed = initparam$fixed,
                         grid = grid,
                         iterations = fitted$counts,
                         likelihood = likelihood,
                         logCompLik = fitted$value,
                         lonlat = lonlat,
                         message = fitted$message,
                         model = model,
                         param = fitted$par,
                         stderr = fitted$stderr,
                         sensmat = fitted$sensmat,
                         varcov = fitted$varcov,
                         varimat = fitted$varimat,
                         type = type)

    structure(c(FitComposite, call = call), class = c("FitComposite"))
  }

print.FitComposite <- function(x, digits = max(3, getOption("digits") - 3), ...)
  {
    dimdata <- dim(x$data)
    
    if(x$grid)
      numdata <- dimdata[3]
    else
      numdata <- dimdata[1]
    
    numcoord <- nrow(x$coord)
    numparam <- length(x$param)

    if(x$likelihood == 'Full')
      {
        method <- 'Likelihood'
        clic <- 'AIC'
      }
    else
      {
        method <- 'Composite-Likelihood'
        clic <- 'CLIC'
      }

    cat('\n##############################################################')
    cat('\nResults: Maximum', method,'Fitting of Random Fields.\n')
    cat('\nSettings:', x$likelihood, method, '\n')
    cat('\nThe density associated to the likelihood objects:', x$model, '\n')
    cat('\nType of the likelihood objects:', x$type, x$method,'\n')
    cat('\nCovariance model:', x$corrmodel, '\n')
    cat('Number of coordinates:', numcoord, '\n')
    cat('Number of observations per location:', numdata, '\n')
    cat('\nMaximum log-', method, ' value: ',
        format(x$logCompLik, digits = digits, nsmall = 2), '\n', sep='')

    if(!is.null(x$clic))
      cat(clic,':', format(x$clic, digits = digits),'\n')

    cat('\nEstimated parameters:\n')
    print.default(x$param, digits = digits, print.gap = 2,
                  quote = FALSE)

    if(!is.null(x$stderr))
      {
        cat('\nStandard errors:\n')
        print.default(x$stderr, digits = digits, print.gap = 2,
                      quote = FALSE)
      }
    
    if(!is.null(x$varcov))
      {
        cat('\nVariance-covariance matrix of the estimates:\n')
        print.default(x$varcov, digits = digits, print.gap = 3,
                      quote = FALSE)
      }

    cat('\n##############################################################')    
    invisible(x)
  }

