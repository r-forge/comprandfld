####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Email: simone.padoan@unibg.it.
### Institute: Department of Information Technology
### and Mathematical Methods, University of Bergamo
### File name: Fitting.r
### Description:
### This file contains a set of procedures
### for maximum composite-likelihood fitting of
### random fields.
### Last change: 07/09/2011.
####################################################


### Procedures are in alphabetical order.

### Fitting procedure:

FitComposite <- function(data, coordx, coordy=NULL, coordt=NULL, corrmodel, fixed=NULL,
                         grid=FALSE, likelihood='Marginal', lonlat=FALSE, margins='Gev',
                         maxdist=NULL, maxtime=NULL, model='Gaussian', optimizer='Nelder-Mead',
                         replicates=1, start=NULL, taper=NULL, threshold=NULL, type='Pairwise',
                         varest=FALSE, vartype='SubSamp', weighted=FALSE, weights=NULL, winconst=NULL)
{
    call <- match.call()
    ### Check the parameters given in input:
    checkinput <- CheckInput(coordx, coordy, coordt, corrmodel, data, fixed, grid,
                             likelihood, lonlat, margins, maxdist, maxtime, model,
                             optimizer, replicates, start, taper, threshold, type,
                             varest, vartype, weighted, weights, winconst)

    if(!is.null(checkinput$error))
      stop(checkinput$error)

    ### If the case set the sub-sampling parameter to the default value
    #if(varest & (vartype=='SubSamp') & (missing(winconst) || !is.numeric(winconst)))
    #  winconst <- 1

    ### Initialization global variables:
    FitComposite <- NULL
    clic <- parscale <- varcov <- stderr <- NULL

    ### Initialization parameters:
    initparam <- WlsInit(coordx, coordy, coordt, corrmodel, data, fixed, grid, likelihood,
                         lonlat, margins, maxdist, maxtime, model, parscale, optimizer=='L-BFGS-B',
                         replicates, start, threshold, type, varest, vartype, weighted, winconst)

    if(!is.null(initparam$error))
      stop(initparam$error)

    ### Model fitting section
    # Full likelihood:
    if(likelihood=='Full')
      {
        # Fitting by log-likelihood maximization:
        fitted <- Likelihood(initparam$corrmodel, initparam$data, initparam$fixed, grid, initparam$lower,
                             initparam$model, initparam$namescorr, initparam$namesnuis, initparam$namesparam,
                             initparam$numcoord, initparam$numrep, initparam$numtime, optimizer, initparam$param,
                             initparam$spacetime, varest, taper, initparam$type, initparam$upper)
      }

    # Composite likelihood:
    if(likelihood=='Marginal' || likelihood=='Conditional')
      {
        vartype <- CheckVarType(vartype)
        fitted <- CompLikelihood(initparam$coordx, initparam$coordy, initparam$corrmodel, initparam$data,
                                 initparam$flagcorr, initparam$flagnuis, initparam$fixed, grid, initparam$likelihood,
                                 lonlat, initparam$lower, initparam$model, initparam$namescorr, initparam$namesnuis,
                                 initparam$namesparam, initparam$numparam, initparam$numparamcorr, optimizer,
                                 initparam$param, initparam$spacetime, initparam$threshold, initparam$type,
                                 initparam$upper, varest, vartype, initparam$winconst)
      }

    # Delete the global variables:
    .C('DeleteGlobalVar', PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)

    ### Set the output object:
    FitComposite <- list(clic = fitted$clic,
                         coordx = initparam$coordx,
                         coordy = initparam$coordy,
                         coordt = initparam$coordt,
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
                         numcoord=initparam$numcoord,
                         numrep=initparam$numrep,
                         numtime=initparam$numtime,
                         param = fitted$par,
                         srange = initparam$srange,
                         stderr = fitted$stderr,
                         sensmat = fitted$sensmat,
                         varcov = fitted$varcov,
                         varimat = fitted$varimat,
                         trange = initparam$trange,
                         type = type)

    structure(c(FitComposite, call = call), class = c("FitComposite"))
  }

print.FitComposite <- function(x, digits = max(3, getOption("digits") - 3), ...)
  {
    if(x$likelihood=='Full')
      {
        method <- 'Likelihood'
        clic <- 'AIC'
      }
    else
      {
        method <- 'Composite-Likelihood'
        clic <- 'CLIC'
      }
    if(x$model=='Gaussian'){ process <- x$model
                             model <- x$model}
     if(x$model=='BinaryGaussian'){ process <- 'Binary'
                             model <- 'Binary Gaussian'}
    if(x$model=='ExtGauss'){ process <- 'Max-Stable'
                             model <- 'Extremal Gaussian'}
    if(x$model=='BrowResn'){ process <- 'Max-Stable'
                             model <- 'Brown-Resnick'}
    if(x$model=='ExtT'){ process <- 'Max-Stable'
                             model <- 'Extremal T'}

    cat('\n##################################################################')
    cat('\nMaximum', method, 'Fitting of', process, 'Random Fields\n')
    cat('\nSetting:', x$likelihood, method, '\n')
    cat('\nModel associated to the likelihood objects:', model, '\n')
    cat('\nType of the likelihood objects:', x$type, x$method,'\n')
    cat('\nCovariance model:', x$corrmodel, '\n')
    cat('Number of spatial coordinates:', x$numcoord, '\n')
    cat('Number of dependent temporal realisations:', x$numtime, '\n')
    cat('Number of replicates of the random field:', x$numrep, '\n')
    cat('Number of estimated parameters:', length(x$param), '\n')
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

    cat('\n##################################################################\n')
    invisible(x)
  }

