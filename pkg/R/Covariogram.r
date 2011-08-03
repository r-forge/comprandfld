####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Email: simone.padoan@unibg.it.
### Institute: University of Bergamo.
### File name: Covariogram.r
### Description:
### This file contains a set of procedures
### to compute and plot the estimated covariance
### function and the variogram after fitting a
### random field by composite-likelihood.
### Last change: 2011/08/03.
####################################################


### Procedures are in alphabetical order.

### Compute and plot the (estimated) covariance function and the variogram
### from a fitted model obtain from the FitComposite or the WLeastSquare procedure
Covariogram <- function(fitted, lags=NULL, answer.cov=FALSE, answer.vario=FALSE,
                        answer.range=FALSE, show.cov=FALSE, show.vario=FALSE,
                        show.range=FALSE, add.cov=FALSE, add.vario=FALSE,
                        pract.range=95, vario=NULL, ...)
  {
    result <- NULL

    if(!class(fitted)=='FitComposite' & !class(fitted)=='WLS')
      {
        cat('Enter an object obtained from fitting a random field with the composite-likelihood or the weigthed least square method\n')
        return(result)
      }

    if(!is.null(vario) & !class(vario)=='Variogram')
      {
        cat('Enter an object obtained from the function EmpVariogram\n')
        return(result)
      }

    if(!is.numeric(pract.range) & answer.range)
      {
        cat('Enter a number for the parameter % of sill\n')
        return(result)
      }
    else
      {
        if(pract.range < 0 || pract.range > 100)
          {
            cat('Entered an incorrect value for the % of sill\n')
            return(result)
          }
        else
          pract.range <- pract.range / 100
      }

    if(is.null(lags))
      {
        numcoord <- nrow(fitted$coord)
        maxlags <- double(1)
        minlags <- double(1)
        .C('RangeDist',  maxlags, minlags,  PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)
        lags <- seq(minlags, maxlags, length=100)
      }

    detect <- DetectParam(fitted$corrmodel, fitted$fixed, fitted$param)
    correlation <- CorrelationFct(detect$corrmodel, lags, detect$param)
    param <- c(fitted$fixed, fitted$param)
    sill <- param['sill']
    nugget <- param['nugget']
    covariance <- nugget + sill * correlation

    if(answer.range || show.range)
      {
        if(fitted$corrmodel=='whittlematern')
          {
            lower <- 1e-10
            upper <- 1e20
          }
        else
          {
            lower <- 0
            upper <- 1e100
          }
        PracticalRange <- function(corrmodel, lags, param, pract.range)
          return(nugget + sill * CorrelationFct(corrmodel, lags, param) -
                 (nugget + sill * (1 - pract.range)))

        Range <- uniroot(PracticalRange, c(lower, upper), corrmodel=detect$corrmodel,
                         param=detect$param, pract.range=pract.range)$root
      }

    if(answer.vario || show.vario)
      variogram <- nugget + sill * (1 - correlation)

    if(show.cov)
      {
        if(add.cov & dev.cur()!=1)
          {
            lines(lags, covariance)
            if(show.range)
              abline(v=Range)
          }
        else
          {
            plot(lags, covariance, type='l', ylim=c(min(covariance), max(covariance)))
            if(show.range)
              abline(v=Range)
          }
      }

    if(show.vario)
      {
        if(add.vario & dev.cur()!=1)
          {
            if(class(vario)=='Variogram')
              points(vario$centers, vario$variogram)

            lines(lags, variogram)
            if(show.range)
              abline(v=Range)
          }
        else
          {
            bnds <- range(variogram)
            if(class(vario)=='Variogram')
              {
                bnds[1] <- min(bnds[1], vario$variogram)
                bnds[2] <- max(bnds[2], vario$variogram)
              }
            plot(lags, variogram, type='l', ylim=c(bnds[1], bnds[2]))
            if(class(vario)=='Variogram')
              points(vario$centers, vario$variogram)

            if(show.range)
              abline(v=Range)
          }
      }

    if(answer.cov)
      result <- list(lags=lags, covariance=covariance)

    if(answer.vario)
      if(!is.list(result))
        result <- list(lags=lags, variogram=variogram)
      else
        result$variogram <- variogram

    if(answer.range)
      if(!is.list(result))
        result <- list(range=Range)
      else
        result$range <- Range

    if(!is.null(result))
      return(result)
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
