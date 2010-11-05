####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Email: simone.padoan@epfl.ch.
### Institute: EPFL.
### File name: Likelihood.r
### Description:
### This file contains a set of procedures
### for maximum likelihood fitting of
### random fields.
### Last change: 06/04/2010.
####################################################


### Procedures are in alphabetical order.



### Log-likelihood procedure for Random Fields

Likelihood <- function(corrmodel, data, fixed, grid, lags, model, namescorr, namesnuis,
                       numcoord, numdata, numpairs, param, type)
  {
    result <- -1.0e8

    if(!CheckParamRange(param))
      return(result)

    param <- c(param, fixed)
    paramcorr <- param[namescorr]
    nuisance <- param[namesnuis]
    
    ### Gaussian model:
    
    if(model == 1)
      {
        stdata <- data - nuisance['mean']
        corr <- double(numpairs)

        .C('CorrelationMat', corr, as.integer(corrmodel), as.double(lags),
          as.integer(numpairs), as.double(paramcorr), PACKAGE='CompRandFld',
          DUP = FALSE, NAOK=TRUE)

        corr <- corr * nuisance['sill']

        varcov <- (nuisance['nugget'] + nuisance['sill']) * diag(numcoord)
        varcov[lower.tri(varcov)] <- corr
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- corr

        cholvarcov <- try(chol(varcov), silent = TRUE)
        if(!is.matrix(cholvarcov))
          return(result)
    
        ivarcov <- chol2inv(cholvarcov)
        if(!is.matrix(ivarcov))
          return(result)

        detvarcov <- sum(log(diag(cholvarcov)))

        if(!grid)
          result <- sum(apply(stdata, 1, LogNormDen, detvarcov=detvarcov,
                              ivarcov=ivarcov, numcoord=numcoord, type=type))

        if(grid)
          {
            result <- NULL
            for(i in 1:numdata)
              result <- result + sum(LogNormDen(c(stdata[,,i]), detvarcov, ivarcov, numcoord, type))
          }
    
        return(result)
      } 
  }

### log of the normal density:

LogNormDen <- function(stdata, detvarcov, ivarcov, numcoord, type)
  {
    if(type == 3)# Restricted log Gaussian density
      {
        sumvarcov <- sum(ivarcov)
        p <- ivarcov - array(rowSums(ivarcov), c(numcoord, 1)) %*% colSums(ivarcov) / sumvarcov
        LogNormDen <- -detvarcov - .5 * (log(sumvarcov) + crossprod(t(crossprod(stdata, p)), stdata) +
                                         (numcoord - 1) * log(2 * pi))
      }

    if(type == 4) # Standard log Gaussian density
      LogNormDen <- -detvarcov - .5 * (crossprod(t(crossprod(stdata, ivarcov)),stdata) + numcoord * log(2 * pi))
    
    return(LogNormDen)
  }

### Optim call for log-likelihood maximization

OptimLik <- function(corrmodel, data, fixed, grid, lags, lower, model,
                     namescorr, namesnuis, namesparam, numcoord, numdata,
                     numpairs, optimizer, param, varest, type, upper)
  {
    if(optimizer=='L-BFGS-B')
      OptimLik <- optim(param, Likelihood, corrmodel=corrmodel, control=list(fnscale=-1,
                        factr=1, pgtol=1e-14, maxit = 1e8), data=data, fixed=fixed,
                        grid=grid, hessian=varest, lags=lags, lower=lower,
                        method=optimizer, model=model, namescorr=namescorr,
                        namesnuis=namesnuis, numcoord=numcoord, numdata=numdata,
                        numpairs=numpairs, upper=upper, type=type)
    else
      OptimLik <- optim(param, Likelihood, corrmodel=corrmodel, control=list(fnscale=-1,
                        reltol=1e-14, maxit=1e8), data=data, fixed=fixed, grid=grid,
                        hessian=varest, lags=lags, method=optimizer, model=model,
                        namescorr=namescorr, namesnuis=namesnuis, numcoord=numcoord,
                        numdata=numdata, numpairs=numpairs, type=type)

    ### Some checks of the output from the optimization procedure:
    if(OptimLik$convergence == 0)
      OptimLik$convergence <- 'Successful'
    else
      if(OptimLik$convergence == 1)
        OptimLik$convergence <- 'Iteration limit reached'
      else
        OptimLik$convergence <- "Optimization may have failed"

    penalty <- length(OptimLik$par)
    OptimLik$clic <- -2 * (OptimLik$value - penalty)

    ### Computation of the variance covariance matrix:
    if(varest)
      {
        # Compute the variance-covariance matrix:
        OptimLik$varcov <- try(solve(-OptimLik$hessian), silent = TRUE)
        OptimLik$sensmat <- NULL
        OptimLik$varimat <- NULL
    
        if(!is.matrix(OptimLik$varcov))
          {
            warning("observed information matrix is singular")
            OptimLik$varcov <- 'none'
            OptimLik$stderr <- 'none'
          }
        else
          {
            dimnames(OptimLik$varcov) <- list(namesparam, namesparam)
            OptimLik$stderr <- diag(OptimLik$varcov)
            if(any(OptimLik$stderr < 0))
              OptimLik$stderr <- 'none'
            else
              {
                OptimLik$stderr <- sqrt(OptimLik$stderr)
                names(OptimLik$stderr) <- namesparam
              }
          }
      }
    return(OptimLik)
  }


