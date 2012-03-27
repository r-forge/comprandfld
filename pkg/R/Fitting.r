####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Emails: simone.padoan@stat.unipd.it,
### moreno.bevilacqua@unibg.it
### Institutions: Department of Statistical Science,
### University of Padua and Department of Information
### Technology and Mathematical Methods, University
### of Bergamo.
### File name: Fitting.r
### Description:
### This file contains a set of procedures
### for maximum composite-likelihood fitting of
### random fields.
### Last change: 02/03/2012.
####################################################


### Procedures are in alphabetical order.

### Fitting procedure:

FitComposite <- function(data, coordx, coordy=NULL, coordt=NULL, corrmodel, fixed=NULL,
                         grid=FALSE, likelihood='Marginal', lonlat=FALSE, margins='Gev',
                         maxdist=NULL, maxtime=NULL, model='Gaussian', optimizer='Nelder-Mead',
                         replicates=1, start=NULL, taper=NULL, threshold=NULL, type='Pairwise',
                         varest=FALSE, vartype='SubSamp', weighted=FALSE, winconst, winstp)
{
    call <- match.call()
    ### Check the parameters given in input:
    checkinput <- CheckInput(coordx, coordy, coordt, corrmodel, data, "Fitting",
                             fixed, grid, likelihood, lonlat, margins, maxdist,
                             maxtime, model, NULL, optimizer, NULL, replicates,
                             start, taper, threshold, type, varest, vartype, weighted)

    if(!is.null(checkinput$error))
      stop(checkinput$error)

    ### Initialization global variables:
    FitComposite <- NULL
    score <- sensmat <- varcov <- varimat <- parscale <- NULL

    ### Initialization parameters:
    initparam <- WlsInit(coordx, coordy, coordt, corrmodel, data, "Fitting", fixed, grid,
                         likelihood, lonlat, margins, maxdist, maxtime, model, NULL, NULL,
                         parscale, optimizer=='L-BFGS-B', replicates, start, threshold,
                         type, varest, vartype, weighted, winconst, winstp)

    if(!is.null(initparam$error))
      stop(initparam$error)
    ### Model fitting section
    # Full likelihood:
    if(likelihood=='Full')
      {
          # Fitting by log-likelihood maximization:
          fitted <- Likelihood(initparam$corrmodel,initparam$data,initparam$fixed,initparam$flagcorr,
                               initparam$flagnuis,grid,initparam$lower,initparam$model,initparam$namescorr,
                               initparam$namesnuis,initparam$namesparam,initparam$numcoord,initparam$numpairs,
                               initparam$numparamcorr,initparam$numrep,initparam$numtime,optimizer,
                               initparam$param,initparam$setup,initparam$spacetime,varest,taper,initparam$type,
                               initparam$upper)
      }
    # Composite likelihood:
    if(likelihood=='Marginal' || likelihood=='Conditional')
      {
          fitted <- CompLikelihood(initparam$coordx, initparam$coordy, initparam$corrmodel, initparam$data,
                                   initparam$flagcorr, initparam$flagnuis, initparam$fixed, grid, initparam$likelihood,
                                   lonlat, initparam$lower, initparam$model, initparam$namescorr, initparam$namesnuis,
                                   initparam$namesparam, initparam$numparam, initparam$numparamcorr, optimizer,
                                   initparam$param, initparam$spacetime, initparam$threshold, initparam$type,
                                   initparam$upper, varest, initparam$vartype, initparam$winconst, initparam$winstp)
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
                         score = fitted$score,
                         srange = initparam$srange,
                         stderr = fitted$stderr,
                         sensmat = fitted$sensmat,
                         varcov = fitted$varcov,
                         varimat = fitted$varimat,
                         vartype = vartype,
                         trange = initparam$trange,
                         threshold = initparam$threshold,
                         type = type,
                         winconst = initparam$winconst,
                         winstp = initparam$winstp)

    structure(c(FitComposite, call = call), class = c("FitComposite"))
  }

print.FitComposite <- function(x, digits = max(3, getOption("digits") - 3), ...)
  {
    if(x$likelihood=='Full'){
        method <- 'Likelihood'
        if(x$type=="Tapering") clic <- "CLIC"
        else clic <- "AIC"}
    else{
        method <- 'Composite-Likelihood'
        clic <- 'CLIC'}
    if(x$model=='Gaussian'){ process <- x$model
                             model <- x$model}
     if(x$model=='BinaryGauss'){ process <- 'Binary'
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

