####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Email: simone.padoan@unibg.it.
### Institute: Department of Information Technology
### and Mathematical Methods, University of Bergamo
### File name: Likelihood.r
### Description:
### This file contains a set of procedures
### for maximum likelihood fitting of
### random fields.
### Last change: 08/09/2011.
####################################################

### Procedures are in alphabetical order.

### Optim call for log-likelihood maximization
Likelihood <- function(corrmodel, data, fixed, grid, lower, model, namescorr,
                       namesnuis, namesparam, numcoord, numrep, numtime, optimizer,
                       param, spacetime, varest, taper, type, upper)
  {
    ### START Define the objective functions
    # Restricted log-likelihood for multivariate normal density:
    LogNormDenRestr <- function(data, varcov, dimat, tapmat)
    {
        llik <- -1.0e8
        # Cholesky decomposition of the covariance matrix:
        cholvarcov <- try(chol(varcov),silent = TRUE)
        if(!is.matrix(cholvarcov)) return(llik)
        # Computes the determinat of the covariance matrix:
        detvarcov <- sum(log(diag(cholvarcov)))
        ivarcov <- chol2inv(cholvarcov)
        sumvarcov <- sum(ivarcov)
        p <- ivarcov - array(rowSums(ivarcov),c(dimat, 1))%*%colSums(ivarcov)/sumvarcov
        llik <- -detvarcov-.5*(log(sumvarcov)+crossprod(t(crossprod(data, p)), data))
        return(llik)
    }
    # Standard log-likelihood function for full normal density:
    LogNormDenStand <- function(data, varcov, dimat, tapmat)
    {
        llik <- -1.0e8
        # Cholesky decomposition of the covariance matrix:
        cholvarcov <- try(chol(varcov),silent = TRUE)
        if(!is.matrix(cholvarcov)) return(llik)
        # Computes the determinat of the covariance matrix:
        detvarcov <- sum(log(diag(cholvarcov)))
        llik <- -detvarcov -.5*(sum(data*backsolve(cholvarcov,
                 forwardsolve(cholvarcov,data,transpose=TRUE,upper.tri=TRUE))))
        return(llik)
    }
    # Tapering log-likelihood for multivariate normal density:
    LogNormDenTap <- function(data, varcov, dimat, tapmat)
    {
        llik <- -1.0e8
        varcov <- varcov*tapmat #tapering matrix
        #using spam package
        varcov <- try(as.spam(varcov),silent=T)
        if(!is.spam(varcov)) return(llik)
        cholvarcov <- try(chol.spam(varcov),silent=T)
        if(!class(cholvarcov)=="spam.chol.NgPeyton") return(llik)
        detvarcov <- c(determinant.spam.chol.NgPeyton(cholvarcov)$modulus)
        ivarcov <- solve.spam(cholvarcov)
        llik <- -detvarcov -.5*(crossprod(t(crossprod(data, ivarcov*tapmat)), data))
        return(llik)
    }
    ### END Define the objective functions
    ### Prepare the data for the objective functions:
    loglik <- function(corr, corrmat, corrmodel, dimat, data, fixed, fname,
                       grid, ident, model, namescorr, namesnuis, param, tapmat)
      {
        loglik <- -1.0e8
        # Set the parameter vector:
        param <- c(param, fixed)
        paramcorr <- param[namescorr]
        nuisance <- param[namesnuis]
        # Standardizes the data:
        stdata <- data-nuisance['mean']
        # Computes the verctor of the correlations:
        .C(corrmat, corr, as.integer(corrmodel), as.double(nuisance),
           as.double(paramcorr), PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)
        if(corr[1]==-2) return(loglik)
        # Set the the format of the correlation matrix:
        corr <- corr*nuisance['sill']
        # Computes the covariance matrix:
        varcov <- (nuisance['nugget']+nuisance['sill'])*ident
        varcov[lower.tri(varcov)] <- corr
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- corr
        # Computes the log-likelihood:
        loglik <- sum(apply(stdata, 1, fname, varcov=varcov, dimat=dimat, tapmat=tapmat))
        return(loglik)
      }
    ### START the main code of the function:
    # Set the data format:
    dimat <- numcoord*numtime
    # Set the constants:
    #nl2pi <- dimat*log(2*pi)
    corr <- double(.5*dimat*(dimat-1))
    corrmat <- 'CorrelationMat'
    if(spacetime){
      dim(data) <- c(numrep, dimat)
      corrmat <- 'CorrelationMat_st'}
    # detects the type of likelihood:
    tapmat <- 1
    if(type==3) fname <- 'LogNormDenRestr'
    if(type==4) fname <- 'LogNormDenStand'
    if(type==5){
        fname <- 'LogNormDenTap'
    ### computing taper matrix
        tapmat <- diag(dimat)
        tapcorr <- double(.5*dimat*(dimat-1))
        corrtap <- switch(taper,
                          Wendland1=15,
                          Wendland2=16,
                          Wendland3=17)
        .C(corrmat, tapcorr, as.integer(corrtap), as.double(c(0,0,1)),
           as.double(1), PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)
        tapmat[lower.tri(tapmat)] <- tapcorr
        tapmat <- t(tapmat)
        tapmat[lower.tri(tapmat)] <- tapcorr}

    if(optimizer=='L-BFGS-B')
      Likelihood <- optim(param, loglik, corr=corr, corrmat=corrmat, corrmodel=corrmodel,
                          control=list(fnscale=-1, factr=1, pgtol=1e-14, maxit=1e8),
                          dimat=dimat, data=data, fname=fname, fixed=fixed, grid=grid,
                          hessian=varest, ident=diag(dimat), lower=lower, method=optimizer,
                          model=model, namescorr=namescorr, namesnuis=namesnuis, upper=upper,
                          tapmat=tapmat)
    else
      Likelihood <- optim(param, loglik, corr=corr, corrmat=corrmat, corrmodel=corrmodel,
                          control=list(fnscale=-1, reltol=1e-14, maxit=1e8), data=data,
                          dimat=dimat, fixed=fixed, fname=fname, grid=grid, hessian=varest,
                          ident=diag(dimat), method=optimizer, model=model, namescorr=namescorr,
                          namesnuis=namesnuis, tapmat=tapmat)

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
    ### END the main code of the function:
    return(Likelihood)
  }


