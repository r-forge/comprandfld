####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Email: simone.padoan@unibg.it.
### Department of Information Technology
### and Mathematical Methods, University of Bergamo.
### File name: Covariogram.r
### Description:
### This file contains a set of procedures
### to compute and plot the estimated covariance
### function and the variogram after fitting a
### random field by composite-likelihood.
### Last change: 07/09/2011.
####################################################


### Procedures are in alphabetical order.

### Compute and plot the (estimated) covariance function and the variogram
### from a fitted model obtain from the FitComposite or the WLeastSquare procedure
Covariogram <- function(fitted, lags=NULL, answer.cov=FALSE, answer.vario=FALSE,
                        answer.extc=FALSE, answer.range=FALSE, show.cov=FALSE,
                        show.vario=FALSE, show.extc=FALSE, show.range=FALSE,
                        add.cov=FALSE, add.vario=FALSE, add.extc=FALSE, pract.range=95,
                        vario=NULL, ...)
  {
    result <- NULL
    # START ---- check the user input --- #
    if(!class(fitted)=='FitComposite' & !class(fitted)=='WLS')
        stop("Enter an object obtained from fitting a random field with the composite-likelihood or the weigthed least square method")

    if(!is.null(vario) & !class(vario)=='Variogram')
        stop("Enter an object obtained from the function EmpVariogram")

    if(!is.numeric(pract.range) & answer.range)
        stop("Enter a number for the parameter % of sill")
    else{
        if(pract.range < 0 || pract.range > 100)
            stop("Entered an incorrect value for the % of sill")
        else
          pract.range <- pract.range / 100}
    # determine the distances:
    if(is.null(lags)) lags <- seq(0,fitted$srange[2],length=100)
    # END ---- check the user input --- #
    ExtremalCoeff <- function(corrmodel, lags, model, nuisance, param)
    {
        numlags <- length(lags)
        extc <- double(numlags)

        .C('ExtCoeff', as.integer(corrmodel), extc, as.double(lags),
           as.integer(model), as.integer(numlags), as.double(nuisance),
           as.double(param), PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)
        return(extc)
     }
    CorrelationFct <- function(corrmodel, lags, param)
    {
        numlags <- length(lags)
        corr <- double(numlags)
        .C('VectCorrelation', corr, as.integer(corrmodel), as.double(lags),
           as.integer(numlags), as.double(param), PACKAGE='CompRandFld',
           DUP = FALSE, NAOK=TRUE)
        return(corr)
    }
    # Check the type of process:
    # START code for the Gaussian random field
    if(fitted$model=="Gaussian"){
        param <- c(fitted$fixed,fitted$param)[CorrelationParam(fitted$corrmodel)]
        nuisance <- c(fitted$fixed,fitted$param)[NuisanceParam(fitted$model)]
        corrmodel <- CheckCorrModel(fitted$corrmodel)
        correlation <- CorrelationFct(corrmodel, lags, param)
        covariance <- nuisance["nugget"]+nuisance["sill"]*correlation
        # set range
        if(answer.range || show.range){
            if(fitted$corrmodel=='whittlematern'){
                lower <- 1e-10
                upper <- 1e20}
            else{ lower <- 0
                  upper <- 1e100}
            PracRangeNorm <- function(corrmodel, lags, nuisance, param, pract.range)
                return(nuisance["nugget"]+nuisance["sill"]*CorrelationFct(corrmodel, lags, param) -
                       (nuisance["nugget"]+nuisance["sill"]*(1-pract.range)))
            Range <- uniroot(PracRangeNorm, c(lower, upper), corrmodel=corrmodel,
                             nuisance=nuisance, param=param, pract.range=pract.range)$root}
        # compute the variogram function
        if(answer.vario || show.vario)
            variogram <- nuisance["nugget"]+nuisance["sill"]*(1-correlation)
        # display the covariance function
        if(show.cov){
            if(add.cov & dev.cur()!=1){
                lines(lags, covariance)
                if(show.range) abline(v=Range)}
            else{
                plot(lags, covariance, type='l', ylim=c(min(covariance), max(covariance)))
                if(show.range) abline(v=Range)}}
        # display the variogram function
        if(show.vario){
            if(add.vario & dev.cur()!=1){
                if(class(vario)=='Variogram') points(vario$centers, vario$variograms)
                lines(lags, variogram)
                if(show.range) abline(v=Range)}
            else{
                bnds <- range(variogram)
                if(class(vario)=='Variogram'){
                    bnds[1] <- min(bnds[1], vario$variograms)
                    bnds[2] <- max(bnds[2], vario$variograms)}
                plot(lags, variogram, type='l', ylim=c(bnds[1], bnds[2]))
                if(class(vario)=='Variogram')
                    points(vario$centers, vario$variograms)

            if(show.range)
              abline(v=Range)}}
        # return the estimated covariance function
        if(answer.cov)
            result <- list(lags=lags, covariance=covariance)
        # return the estimated variogram function
        if(answer.vario)
            if(!is.list(result)) result <- list(lags=lags, variogram=variogram)
            else result$variogram <- variogram
        # return the range
        if(answer.range)
            if(!is.list(result)) result <- list(range=Range)
            else result$range <- Range}
    # END code for the Gaussian random field
    # START code for the max-stable random field
    if(fitted$model %in% c("BrownResn", "ExtGauss", "ExtT")){# Extremal processes
        param <- c(fitted$fixed,fitted$param)[CorrelationParam(fitted$corrmodel)]
        nuisance <- c(fitted$fixed,fitted$param)[NuisanceParam(fitted$model)]
        corrmodel <- CheckCorrModel(fitted$corrmodel)
        model <- CheckModel(fitted$model)
        # compute the extremal coefficient function
        extrcoeff <- ExtremalCoeff(corrmodel, lags, model, nuisance, param)
        # set range
        if(answer.range || show.range){
            if(fitted$corrmodel=='whittlematern'){
                lower <- 1e-10
                upper <- 1e20}
            else{ lower <- 0
                  upper <- 1e100}
            PracRangeExt <- function(corrmodel, lags, model, nuisance, param, pract.range)
                return(ExtremalCoeff(corrmodel, lags, model, nuisance, param)-(2*pract.range))
            Range <- uniroot(PracRangeExt, c(lower, upper), corrmodel=corrmodel, model=model,
                             nuisance=nuisance, param=param, pract.range=pract.range)$root}

        if(answer.cov || show.cov){
            correlation <- CorrelationFct(corrmodel, lags, param)
            covariance <- 1-nuisance["sill"]+nuisance["sill"]*correlation}

        # return the estimated covariance function
        if(answer.cov)
            result <- list(lags=lags, covariance=covariance)
        # return the estimated extremal coefficient function
        if(answer.extc)
            if(!is.list(result)) result <- list(lags=lags, extrcoeff=extrcoeff)
            else result$extrcoeff <- extrcoeff
        # display the covariance function
        if(show.cov){
            if(add.cov & dev.cur()!=1){
                lines(lags, covariance)
                if(show.range) abline(v=Range)}
            else{
                plot(lags, covariance, type='l', ylim=c(min(covariance), max(covariance)))
                if(show.range) abline(v=Range)}}
        # display the extremal function
        if(show.extc){
            if(add.extc & dev.cur()!=1){
                if(class(vario)=='Variogram') points(vario$centers, vario$extcoeff)
                lines(lags, extrcoeff)
                if(show.range) abline(v=Range)}
            else{
                bnds <- range(extrcoeff)
                if(class(vario)=='Variogram'){
                    bnds[1] <- min(bnds[1], vario$extcoeff)
                    bnds[2] <- max(bnds[2], vario$extcoeff)}
                plot(lags, extrcoeff, type='l', ylim=c(bnds[1], bnds[2]))
                if(class(vario)=='Variogram')
                    points(vario$centers, vario$extcoeff)

            if(show.range)
              abline(v=Range)}}
        }
    # END code for the max-stable random field
    if(!is.null(result))
      return(result)
  }
