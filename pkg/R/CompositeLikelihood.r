####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Email: simone.padoan@eofl.ch.
### Institute: EPFL.
### File name: CompositeLikelihood.r
### Description:
### This file contains a set of procedures
### for maximum composite-likelihood fitting of
### random fields.
### Last change: 04/06/2010.
####################################################


### Procedures are in alphabetical order.


CompLikelihood <- function(coordx, coordy, corrmodel, data, fixed, model,
                           namescorr, namesnuis, numcoord, numdata, param, type)
  {
    result <- -1.0e15

    if(!CheckParamRange(param))
      return(result)

    result <- double(1)
    param <- c(param, fixed)
    paramcorr <- param[namescorr]
    nuisance <- param[namesnuis]

    .C('CompLikelihood', as.double(coordx), as.double(coordy),
       as.integer(corrmodel), as.double(data), as.integer(model),
       as.double(nuisance), as.integer(numdata), as.integer(numcoord),
       as.double(paramcorr), result, as.integer(type), PACKAGE='CompRandFld',
       DUP = FALSE, NAOK=TRUE)
    
    return(result)
  }

CompScore <- function(coordx, coordy, corrmodel, data, fixed, flag, flagcorr,
                      numcoord, numdata, numparam, numparamcorr, param, type, weight)
  {
    result <- rep(1.0e15, numparam)

    if(!CheckParamRange(param))
      return(result)

    eps <- (.Machine$double.eps)^(1/3)
    result <- double(numparam)
    param <- c(param, fixed)
    namesparam <- sort(names(param))
    param <- param[namesparam]
    par <- param[c('mean','nugget','sill')]
    param <- param[-match(c('mean','nugget','sill'), namesparam)]
 
    .C('CompScore', as.double(coordx), as.double(coordy), as.integer(corrmodel),
       as.double(data), as.double(eps), as.integer(flag), as.integer(flagcorr),
       as.integer(numdata), as.integer(numparamcorr), as.integer(numparam),
       as.integer(numcoord), as.double(par), as.double(param), result, as.integer(type),
       as.integer(weight), PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)
    
    return(result)
  }

### Optim call for Composite log-likelihood maximization

OptimCompLik <- function(coordx, coordy, corrmodel, data, flagcorr, flagnuis, fixed, grid,
                         hessian, lower, model, namemodel, namescorr, namesnuis, namesparam,
                         namessim, numcoord, numdata, numparam, numparamcorr, optimizer, param,
                         varest, type, upper)
  {
    if(optimizer=='L-BFGS-B')
      OptimCompLik <- optim(param, CompLikelihood, coordx=coordx, coordy=coordy, corrmodel=corrmodel,
                            control=list(fnscale=-1, factr=1, pgtol=1e-14, maxit = 1e8), data=data,
                            fixed=fixed, hessian=hessian, lower=lower, method=optimizer, model=model,
                            namescorr=namescorr, namesnuis=namesnuis, numcoord=numcoord, numdata=numdata,
                            type=type, upper=upper)
    else
      OptimCompLik <- optim(param, CompLikelihood, coordx=coordx, coordy=coordy, corrmodel=corrmodel,
                            control=list(fnscale=-1, reltol=1e-14, maxit=1e8), data=data, flagcorr,
                            fixed=fixed, hessian=hessian, method=optimizer, model=model,
                            namescorr=namescorr, namesnuis=namesnuis, numcoord=numcoord,
                            numdata=numdata, type=type)

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
            eps <- (.Machine$double.eps)^(1/3)
            godambe <- double(2 * dimmat)
            param <- c(OptimCompLik$par, fixed)
                
            paramcorr <- param[namescorr]
            nuisance <- param[namesnuis]
                   
            if((numdata / numcoord) > 1)
              {                
                # Compute the sensitivity and variability matrices:
                .C('GodambeMat', as.double(coordx), as.double(coordy), as.integer(corrmodel), as.double(data),
                   as.double(eps), as.integer(flagcorr), as.integer(flagnuis), as.integer(model), as.integer(numdata),
                   as.integer(numparam), as.integer(numparamcorr),as.integer(numcoord), as.double(paramcorr),
                   as.double(nuisance), godambe, as.integer(type), PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)
              }
            else
              {
                # The sensitivity (H) and the variability matrix (J)
                # are estimated by Monte Carlo methods (bootstrap and Jackknife):
                
                numsim <- 500
                param <- param[namessim]
                mcgodambe <- double(2 * dimmat)
                #mcsample <- matrix(numeric(numsim * numparam), ncol=numparam, nrow=numsim)
                for(i in 1 : numsim)
                  {
                    # Simulate the random field:
                    sim <- GaussRF(x=coordx, y=coordy, model=namemodel, grid=grid,
                                   param=param, n=numdata, pch='')
                     # Compute the sensitivity and variability matrices:
                    .C('GodambeMat', as.double(coordx), as.double(coordy), as.integer(corrmodel), as.double(sim),
                       as.double(eps), as.integer(flagcorr), as.integer(flagnuis), as.integer(model), as.integer(numdata),
                       as.integer(numparam), as.integer(numparamcorr),as.integer(numcoord), as.double(paramcorr),
                       as.double(nuisance), godambe, as.integer(type), PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)

                    mcgodambe <- mcgodambe + godambe
                  }
                godambe <- mcgodambe / numsim
              }
            # Set sensitivity matrix:
            OptimCompLik$sensmat <- matrix(godambe[1 : dimmat], ncol=numparam)
            # Set variability matrix:
            OptimCompLik$varmat <- matrix(godambe[(dimmat + 1) : (2 * dimmat)], ncol=numparam)
            namesgod <- c(namesnuis[as.logical(flagnuis)], namescorr[as.logical(flagcorr)])
            dimnames(OptimCompLik$sensmat) <- list(namesgod, namesgod)
            dimnames(OptimCompLik$varmat) <- list(namesgod, namesgod)
            OptimCompLik$sensmat <- OptimCompLik$sensmat[namesparam, namesparam]
            OptimCompLik$varmat <- OptimCompLik$varmat[namesparam, namesparam]
            isensmat <- try(solve(OptimCompLik$sensmat), silent = TRUE)
    
            if(!is.matrix(isensmat) || !is.matrix(OptimCompLik$varmat))
              {
                warning("observed information matrix is singular")
                OptimCompLik$varcov <- 'none'
                OptimCompLik$stderr <- 'none'
              }
            else
              {
                penalty <- OptimCompLik$varmat %*% isensmat
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


