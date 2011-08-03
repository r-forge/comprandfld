####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Email: simone.padoan@unibg.it.
### Institute: University of Bergamo.
### File name: Likelihood.r
### Description:
### This file contains a set of procedures
### for maximum likelihood fitting of
### random fields.
### Last change: 2011/08/03.
####################################################


### Procedures are in alphabetical order.

### Optim call for log-likelihood maximization
Likelihood <- function(corrmodel, data, fixed, grid, lower, model, namescorr,
                     namesnuis, namesparam, numcoord, numdata, numpairs,
                     optimizer, param, varest, type, upper)
  {
    ### Define the object function:
    loglik <- function(corrmodel, data, fixed, grid, model, namescorr,
                       namesnuis, numcoord, numdata, numpairs, param, type)
      {
        loglik <- -1.0e8

        if(!CheckParamRange(param))
          return(loglik)

        param <- c(param, fixed)
        paramcorr <- param[namescorr]
        nuisance <- param[namesnuis]

        stdata <- data - nuisance['mean']
        corr <- double(numpairs)

        .C('CorrelationMat', corr, as.integer(corrmodel), as.integer(numpairs),
           as.double(paramcorr), PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)

        corr <- corr * nuisance['sill']

        varcov <- (nuisance['nugget'] + nuisance['sill']) * diag(numcoord)
        varcov[lower.tri(varcov)] <- corr
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- corr

        cholvarcov <- try(chol(varcov), silent = TRUE)
        if(!is.matrix(cholvarcov))
          return(loglik)

        ivarcov <- chol2inv(cholvarcov)
        if(!is.matrix(ivarcov))
          return(loglik)

        detvarcov <- sum(log(diag(cholvarcov)))

        if(!grid)
          loglik <- sum(apply(stdata, 1, LogNormDen, detvarcov=detvarcov,
                              ivarcov=ivarcov, numcoord=numcoord, type=type))

        if(grid)
          {
            loglik <- NULL
            for(i in 1:numdata)
              loglik <- loglik + sum(LogNormDen(c(stdata[,,i]), detvarcov, ivarcov, numcoord, type))
          }

        return(loglik)
      }

    if(optimizer=='L-BFGS-B')
      Likelihood <- optim(param, loglik, corrmodel=corrmodel, control=list(fnscale=-1,
                          factr=1, pgtol=1e-14, maxit = 1e8), data=data, fixed=fixed,
                          grid=grid, hessian=varest, lower=lower, method=optimizer,
                          model=model, namescorr=namescorr, namesnuis=namesnuis,
                          numcoord=numcoord, numdata=numdata, numpairs=numpairs,
                          upper=upper, type=type)
    else
      Likelihood <- optim(param, loglik, corrmodel=corrmodel, control=list(fnscale=-1,
                          reltol=1e-14, maxit=1e8), data=data, fixed=fixed, grid=grid,
                          hessian=varest, method=optimizer, model=model, namescorr=namescorr,
                          namesnuis=namesnuis, numcoord=numcoord, numdata=numdata,
                          numpairs=numpairs, type=type)

    ### Some checks of the output from the optimization procedure:
    if(Likelihood$convergence == 0)
      Likelihood$convergence <- 'Successful'
    else
      if(Likelihood$convergence == 1)
        Likelihood$convergence <- 'Iteration limit reached'
      else
        Likelihood$convergence <- "Optimization may have failed"

    penalty <- length(Likelihood$par)
    Likelihood$clic <- -2 * (Likelihood$value - penalty)

    ### Computation of the variance covariance matrix:
    if(varest)
      {
        # Compute the variance-covariance matrix:
        Likelihood$varcov <- try(solve(-Likelihood$hessian), silent = TRUE)
        Likelihood$sensmat <- NULL
        Likelihood$varimat <- NULL

        if(!is.matrix(Likelihood$varcov))
          {
            warning("observed information matrix is singular")
            Likelihood$varcov <- 'none'
            Likelihood$stderr <- 'none'
          }
        else
          {
            dimnames(Likelihood$varcov) <- list(namesparam, namesparam)
            Likelihood$stderr <- diag(Likelihood$varcov)
            if(any(Likelihood$stderr < 0))
              Likelihood$stderr <- 'none'
            else
              {
                Likelihood$stderr <- sqrt(Likelihood$stderr)
                names(Likelihood$stderr) <- namesparam
              }
          }
      }
    return(Likelihood)
  }

### Standard and restricted log-likelihood for multivariate normal density:
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



