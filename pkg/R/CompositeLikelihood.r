####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Email: simone.padoan@epfl.ch.
### Institute: EPFL.
### File name: CompositeLikelihood.r
### Description:
### This file contains a set of procedures
### for maximum composite-likelihood fitting of
### random fields.
### Last change: 12/11/2010.
####################################################


### Procedures are in alphabetical order.


CompLikelihood <- function(coordx, coordy, corrmodel, data, fixed, likelihood,
                           model, namescorr, namesnuis, numcoord, numdata, param, type)
  {
    result <- -1.0e15

    if(!CheckParamRange(param))
      return(result)

    result <- double(1)
    param <- c(param, fixed)
    paramcorr <- param[namescorr]
    nuisance <- param[namesnuis]

    .C('CompLikelihood', as.double(coordx), as.double(coordy), as.integer(corrmodel),
       as.double(data), as.integer(likelihood), as.integer(model), as.double(nuisance),
       as.integer(numdata), as.integer(numcoord), as.double(paramcorr), result,
       as.integer(type), PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)
    
    return(result)
  }

### Optim call for Composite log-likelihood maximization

OptimCompLik <- function(coordx, coordy, corrmodel, data, flagcorr, flagnuis, fixed, grid,
                         likelihood, lonlat, lower, model, namescorr, namesnuis, namesparam,
                         numcoord, numdata, numparam, numparamcorr, optimizer, param, type,
                         upper, varest, vartype, winconst)
  {
   
    if(optimizer=='L-BFGS-B')
      OptimCompLik <- optim(param, CompLikelihood, coordx=coordx, coordy=coordy, corrmodel=corrmodel,
                            control=list(fnscale=-1, factr=1, pgtol=1e-14, maxit = 1e8), data=data,
                            fixed=fixed, hessian=FALSE, likelihood=likelihood, lower=lower,
                            method=optimizer, model=model, namescorr=namescorr, namesnuis=namesnuis,
                            numcoord=numcoord, numdata=numdata, type=type, upper=upper)
    else
      OptimCompLik <- optim(param, CompLikelihood, coordx=coordx, coordy=coordy, corrmodel=corrmodel,
                            control=list(fnscale=-1, reltol=1e-14, maxit=1e8), data=data, fixed=fixed,
                            hessian=FALSE, likelihood=likelihood, method=optimizer, model=model,
                            namescorr=namescorr, namesnuis=namesnuis, numcoord=numcoord, numdata=numdata,
                            type=type)
    
    if(OptimCompLik$convergence == 0)
      OptimCompLik$convergence <- 'Successful'
    else
      if(OptimCompLik$convergence == 1)
        OptimCompLik$convergence <- 'Iteration limit reached'
      else
        OptimCompLik$convergence <- "Optimization may have failed"

        ### Computation of the variance-covariance matrix:

        if(varest)
          {
            # The sensitivity (H) and the variability (J) matrices
            # are estimated by the sample estimators contro-parts:

            dimmat <- numparam^2
            dmat <- numparam * (numparam + 1) / 2
            eps <- (.Machine$double.eps)^(1/3)
            param <- c(OptimCompLik$par, fixed)
                
            paramcorr <- param[namescorr]
            nuisance <- param[namesnuis]
            
            sensmat <- double(dmat)
            varimat <- double(dmat)
 
            # Set the window parameter:

            .C('GodambeMat', as.double(coordx), as.double(coordy), as.integer(corrmodel), as.double(data),
               as.double(eps), as.integer(flagcorr), as.integer(flagnuis), as.integer(likelihood), as.integer(lonlat),
               as.integer(model), as.integer(numdata), as.integer(numparam), as.integer(numparamcorr),
               as.integer(numcoord), as.double(paramcorr), as.double(nuisance), sensmat, as.integer(type),
               varimat, as.integer(vartype), as.double(winconst), PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)
          
            # Set sensitivity matrix:
            OptimCompLik$sensmat <- matrix(double(dimmat), ncol=numparam)
            OptimCompLik$sensmat[lower.tri(OptimCompLik$sensmat, diag=TRUE)] <- sensmat
            OptimCompLik$sensmat <- t(OptimCompLik$sensmat)
            OptimCompLik$sensmat[lower.tri(OptimCompLik$sensmat, diag=TRUE)] <- sensmat
            
            # Set variability matrix:
            OptimCompLik$varimat <- matrix(double(dimmat), ncol=numparam)
            OptimCompLik$varimat[lower.tri(OptimCompLik$varimat, diag=TRUE)] <- varimat
            OptimCompLik$varimat <- t(OptimCompLik$varimat)
            OptimCompLik$varimat[lower.tri(OptimCompLik$varimat, diag=TRUE)] <- varimat
                 
            namesgod <- c(namesnuis[as.logical(flagnuis)], namescorr[as.logical(flagcorr)])
            dimnames(OptimCompLik$sensmat) <- list(namesgod, namesgod)
            dimnames(OptimCompLik$varimat) <- list(namesgod, namesgod)
            OptimCompLik$sensmat <- OptimCompLik$sensmat[namesparam, namesparam]
            OptimCompLik$varimat <- OptimCompLik$varimat[namesparam, namesparam]
            isensmat <- try(solve(OptimCompLik$sensmat), silent = TRUE)
    
            if(!is.matrix(isensmat) || !is.matrix(OptimCompLik$varimat))
              {
                warning("observed information matrix is singular")
                OptimCompLik$varcov <- 'none'
                OptimCompLik$stderr <- 'none'
              }
            else
              {
                penalty <- OptimCompLik$varimat %*% isensmat
                OptimCompLik$clic <- -2 * (OptimCompLik$value - sum(diag(penalty)))
        
                OptimCompLik$varcov <- isensmat %*% penalty
                OptimCompLik$stderr <- diag(OptimCompLik$varcov)
                if(any(OptimCompLik$stderr < 0))
                  OptimCompLik$stderr <- 'none'
                else
                  OptimCompLik$stderr <- sqrt(OptimCompLik$stderr)
              }
          }
    return(OptimCompLik)
  }

