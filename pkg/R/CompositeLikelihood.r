####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Email: simone.padoan@unibg.it.
### Institute: University of Bergamo.
### File name: CompositeLikelihood.r
### Description:
### This file contains a set of procedures
### for maximum composite-likelihood fitting of
### random fields.
### Last change: 2011/08/03.
####################################################


### Procedures are in alphabetical order.

### Optim call for Composite log-likelihood maximization

CompLikelihood <- function(coordx, coordy, corrmodel, data, flagcorr, flagnuis, fixed, grid,
                         likelihood, lonlat, lower, model, namescorr, namesnuis, namesparam,
                         numcoord, numdata, numparam, numparamcorr, optimizer, param, type,
                         upper, varest, vartype, winconst)
  {
    ### Define the object function:
    comploglik <- function(corrmodel, data, fixed, fun, namescorr,
                           namesnuis, numcoord, numdata, param)
      {
        result <- -1.0e15

        if(!CheckParamRange(param))
          return(result)

        result <- double(1)
        param <- c(param, fixed)
        paramcorr <- param[namescorr]
        nuisance <- param[namesnuis]

        .C(fun, as.integer(corrmodel), as.double(data), as.double(nuisance),
           as.integer(numdata), as.integer(numcoord), as.double(paramcorr),
           result, PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)

        return(result)
      }
    fname <- NULL

    if(likelihood == 1 & type == 2)
      fname <- 'Comp_Cond_Gauss'
    if(likelihood == 3 & type == 1)
      fname <- 'Comp_Diff_Gauss'
    if(likelihood == 3 & type == 2)
      fname <- 'Comp_Pair_Gauss'

    if(optimizer=='L-BFGS-B')
      CompLikelihood <- optim(param, comploglik, corrmodel=corrmodel, control=list(fnscale=-1,
                            factr=1, pgtol=1e-14, maxit=1e8), data=data, fixed=fixed,
                            fun=fname, hessian=FALSE, lower=lower, method=optimizer,
                            namescorr=namescorr, namesnuis=namesnuis, numcoord=numcoord,
                            numdata=numdata, upper=upper)
    else
      CompLikelihood <- optim(param, comploglik, corrmodel=corrmodel, control=list(fnscale=-1,
                            reltol=1e-14, maxit=1e8), data=data, fixed=fixed, fun=fname,
                            hessian=FALSE, method=optimizer, namescorr=namescorr,
                            namesnuis=namesnuis, numcoord=numcoord, numdata=numdata)

    if(CompLikelihood$convergence == 0)
      CompLikelihood$convergence <- 'Successful'
    else
      if(CompLikelihood$convergence == 1)
        CompLikelihood$convergence <- 'Iteration limit reached'
      else
        CompLikelihood$convergence <- "Optimization may have failed"

        ### Computation of the variance-covariance matrix:

        if(varest)
          {
            # The sensitivity (H) and the variability (J) matrices
            # are estimated by the sample estimators contro-parts:

            dimmat <- numparam^2
            dmat <- numparam * (numparam + 1) / 2
            eps <- (.Machine$double.eps)^(1/3)
            param <- c(CompLikelihood$par, fixed)

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
            CompLikelihood$sensmat <- matrix(double(dimmat), ncol=numparam)
            CompLikelihood$sensmat[lower.tri(CompLikelihood$sensmat, diag=TRUE)] <- sensmat
            CompLikelihood$sensmat <- t(CompLikelihood$sensmat)
            CompLikelihood$sensmat[lower.tri(CompLikelihood$sensmat, diag=TRUE)] <- sensmat

            # Set variability matrix:
            CompLikelihood$varimat <- matrix(double(dimmat), ncol=numparam)
            CompLikelihood$varimat[lower.tri(CompLikelihood$varimat, diag=TRUE)] <- varimat
            CompLikelihood$varimat <- t(CompLikelihood$varimat)
            CompLikelihood$varimat[lower.tri(CompLikelihood$varimat, diag=TRUE)] <- varimat

            namesgod <- c(namesnuis[as.logical(flagnuis)], namescorr[as.logical(flagcorr)])
            dimnames(CompLikelihood$sensmat) <- list(namesgod, namesgod)
            dimnames(CompLikelihood$varimat) <- list(namesgod, namesgod)
            CompLikelihood$sensmat <- CompLikelihood$sensmat[namesparam, namesparam]
            CompLikelihood$varimat <- CompLikelihood$varimat[namesparam, namesparam]
            isensmat <- try(solve(CompLikelihood$sensmat), silent = TRUE)

            if(!is.matrix(isensmat) || !is.matrix(CompLikelihood$varimat))
              {
                warning("observed information matrix is singular")
                CompLikelihood$varcov <- 'none'
                CompLikelihood$stderr <- 'none'
              }
            else
              {
                penalty <- CompLikelihood$varimat %*% isensmat
                CompLikelihood$clic <- -2 * (CompLikelihood$value - sum(diag(penalty)))

                CompLikelihood$varcov <- isensmat %*% penalty
                CompLikelihood$stderr <- diag(CompLikelihood$varcov)
                if(any(CompLikelihood$stderr < 0))
                  CompLikelihood$stderr <- 'none'
                else
                  CompLikelihood$stderr <- sqrt(CompLikelihood$stderr)
              }
          }
    return(CompLikelihood)
  }

