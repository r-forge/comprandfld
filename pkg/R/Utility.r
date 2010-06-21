####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Email: simone.padoan@eofl.ch.
### Institute: EPFL.
### File name: Utility.r
### Description:
### This file contains a set of procedures
### for supporting all the other functions.
### Last change: 01/04/2010.
####################################################

### Procedures are in alphabetical order.

CheckCorrModel <- function(corrmodel)
  {
    CheckCorrModel <- NULL
    # Correlation function are in alphabetical order
    CheckCorrModel <- switch(corrmodel,
                      cauchy=1,
                      exponential=2,
                      gauss=3,
                      gencauchy=4,
                      stable=5,
                      whittlematern=6)

    
    return(CheckCorrModel)
  }

CheckInput <- function(coordx, coordy, corrmodel, data, fixed, grid, likelihood,
                       lonlat, model, optimizer, start, varest, time, type, weighted,
                       weights)
  {
    error <- NULL

    # Check if the input is inserted correctly
    
    if(missing(coordx) || !is.numeric(coordx))
      {
        error <- 'insert a suitable set of numeric coordinates\n'
        return(list(error=error))
      }

    if(!is.null(coordy) & !is.numeric(coordy))
      {
        error <- 'insert a suitable set of numeric coordinates\n'
        return(list(error=error))
      }
    
    if(missing(corrmodel) || !is.character(corrmodel))
      {
        error <- 'insert the correlation model\n'
        return(list(error=error))
      }
      
    if(missing(data) || !is.numeric(data))
      {
        error <- 'insert a numeric vector or matrix of data\n'
        return(list(error=error))
      }

    if(!is.null(fixed) & !is.list(fixed))
      {
        error <- 'insert fixed values as a list of parameters\n'
        return(list(error=error))
      }

    if(!is.null(grid) & !is.logical(grid))
      {
        error <- 'the parameter grid need to be a logic value\n'
        return(list(error=error))
      }
        
    if(!is.null(likelihood) & !is.character(likelihood))
      {
        error <- 'insert the type of likelihood objects\n'
        return(list(error=error))
      }

    if(!is.null(model) & !is.character(model))
      {
        error <- 'insert the model of the marginals or conditionals likelihoods\n'
        return(list(error=error))
      }

    if(!is.null(optimizer) & !is.character(optimizer))
      {
        error <- 'insert the type of maximising algorithm\n'
        return(list(error=error))
      }

    if(!is.null(start))
      {
        if(!is.list(start))
          {
            error <- 'insert starting values as a list of parameters\n'
            return(list(error=error))
          }

        if(any(names(start)=='mean') & (type=='Difference' || type=='Restricted'))
          {
            error <- 'the mean parameter is not allow with the difference composite likelihood\n'
            return(list(error=error))
          }
      }

    if(!is.null(varest) & !is.logical(varest))
      {
        error <- 'the parameter std.err need to be a logical value\n'
        return(list(error=error))
      }

    if(!is.null(time) & !is.logical(time))
      {
        error <- 'the parameter time need to be a logical value\n'
        return(list(error=error))
      }

    if(!is.null(type) & !is.character(type))
      {
        error <- 'insert the configuration of the likelihood objects\n'
        return(list(error=error))
      }

    if(!is.null(lonlat) & !is.logical(lonlat))
      {
        error <- 'insert the type of coordinates'
        return(list(error=error))
      }
    
    if(!is.null(weighted) & !is.logical(weighted))
      {
        error <- 'insert if the composite likelihood need to be weighted'
        return(list(error=error))
      }

    if(!is.null(weights) & !is.numeric(weights))
      {
        error <- 'insert a numeric vector or matrix of weights'
        return(list(error=error))
      }

    # Check the correctness of the inserted input
 
    if(is.null(CheckCorrModel(corrmodel)))
      {
        error <- 'the name of the correlation model is not correct\n'
        return(list(error=error))
      }

    if(!is.null(fixed))
      {
        if(!CheckParam(corrmodel, names(fixed), length(fixed)))
          {
            error <- 'some names of the fixed parameters is/are not correct\n'
            return(list(error=error))
          }
        
        if(!CheckParamRange(unlist(fixed)))
          {
            error <- 'some fixed values are out of the range\n'
            return(list(error=error))
          }
      }
    
    checklik <- CheckLikelihood(likelihood)
    
    if(is.null(checklik))
      {
        error <- 'the setting name of the likelihood objects is not correct\n'
        return(list(error=error))
      }

    if(is.null(CheckModel(model)))
      {
        error <- 'the model name of the likelihood objects is not correct\n'
        return(list(error=error))
      }

    if(!is.null(start))
      {
        if(!CheckParam(corrmodel, names(start), length(start)))
          {
            error <- 'some names of the starting parameters is/are not correct\n'
            return(list(error=error))
          }
        
        if(!CheckParamRange(unlist(start)))
          {
            error <- 'some starting values are out of the range\n'
            return(list(error=error))
          }
      }

    checktype <- CheckType(type)
    
    if(is.null(checktype))
      {
        error <- 'the type name of the likelihood objects is not correct\n'
        return(list(error=error))
      }

    if(checklik == 2)
      {
        if(!any(checktype == c(3, 4)))
          {
            error <- 'insert a type name of the likelihood objects compatible with the full likelihood'
            return(list(error=error))
          }
      }

    if(checklik == 3)
      {
        if(!any(checktype == c(1, 2)))
          {
            error <- 'insert a type name of the likelihood objects compatible with the composite-likelihood'
            return(list(error=error))
          }
      }

    dimdata <- dim(data)
    if(is.null(dimdata) & time)
      {
        error <- c('insert a numeric matrix of observations')
        return(list(error=error))
      }

    if(is.null(coordy))
      {
        if(ncol(coordx) != 2)
          {
            error <- ('insert a d x 2 matrix of coordinates\n')
            return(list(error=error))
          }
                
        if(time)
          {
            numcoord <- nrow(coordx)      
            if(grid)
              {
                numcoordx <- sqrt(numcoord)
                if(!is.integer(numcoordx))
                  {
                    error <- c('the format of the coordinates is not correct\n')
                    return(list(error=error))
                  }
                if(numcoordx != dimdata[1] || numcoordx != dimdata[2])
                  {
                    error <- c('the format of the data is not correct\n')
                    return(list(error=error))
                  }
              }
            else
              if(numcoord != dimdata[2])
                {
                  error <- c('the format of the data is not correct\n')
                  return(list(error=error))
                 }
          }
      }
    else
      {
        if(time)
          {
            numcoordx <- length(coordx)
            numcoordy <- length(coordy)
            
            if(grid)
              {
                if(numcoordx != dimdata[1] || numcoordy != dimdata[2])
                  {
                    error <- c('the format of the data is not correct\n')
                    return(list(error=error))
                  }
              }
            else
              {
                if(numcoordx != numcoordy)
                  {
                    error <- c('the number of coordinates is not the same\n')
                    return(list(error=error))
                  }
                if(numcoordx != dimdata[2])
                  {
                    error <- c('the format of the data is not correct\n')
                    return(list(error=error))
                  }
              }
          }
      }
  }

CheckLikelihood <- function(likelihood)
  {
    CheckLikelihood <- switch(likelihood,
                              None=0,
                              Conditional=1,
                              Full=2,
                              Marginal=3)

    return(CheckLikelihood)
  }

CheckModel <- function(model)
  {
    CheckModel <- switch(model,
                         None=0,
                         Gaussian=1)

    return(CheckModel)
  }

CheckParam <- function(corrmodel, namesparam, numparam)
  {
    for(i in 1 : numparam)
      {
        if(corrmodel=='exponential' || corrmodel=='gauss')
          if(is.null(switch(namesparam[i],
                            mean=2,
                            nugget=3,
                            scale=4,
                            sill=5)))
            return(FALSE)

        if(corrmodel=='stable')
          if(is.null(switch(namesparam[i],
                            mean=2,
                            nugget=3,
                            power=4,
                            scale=5,
                            sill=6)))
            return(FALSE)

        if(corrmodel=='cauchy')
          if(is.null(switch(namesparam[i],
                            mean=2,
                            nugget=3,
                            power2=4,
                            scale=5,
                            sill=6)))
            return(FALSE)
        
        if(corrmodel=='gencauchy')
          if(is.null(switch(namesparam[i],
                            mean=2,
                            nugget=3,
                            power1=4,
                            power2=5,
                            scale=6,
                            sill=7)))
            return(FALSE)

        if(corrmodel=='whittlematern')
          if(is.null(switch(namesparam[i],
                            mean=2,
                            nugget=3,
                            scale=4,
                            sill=5,
                            smooth=6)))
            return(FALSE)
      }

    return(TRUE)
  }

CheckParamRange <- function(param)
  {
    if(!is.na(param['nugget'])) if(param['nugget'] < 0) return(FALSE)
    if(!is.na(param['power'])) if(param['power'] <=0 || param['power'] > 2) return(FALSE)
    if(!is.na(param['power1'])) if(param['power1'] <=0 || param['power1'] > 2) return(FALSE)
    if(!is.na(param['power2'])) if(param['power2'] <= 0) return(FALSE)
    if(!is.na(param['scale'])) if(param['scale'] <= 0) return(FALSE)
    if(!is.na(param['sill'])) if(param['sill'] <= 0) return(FALSE)
    if(!is.na(param['smooth'])) if(param['smooth'] <= 0) return(FALSE)

    return(TRUE)
  }

CheckType <- function(type)
  {
    CheckType <- switch(type,
                        Difference=1,
                        Pairwise=2,
                        Restricted=3,
                        Standard=4,
                        WLeastSquare=5)
    return(CheckType)
  }

CorrelationParam <- function(corrmodel)
  {
    namesparam <- NULL
    
    if(corrmodel=='cauchy')
      namesparam <- c('power2', 'scale')
    
    if(corrmodel=='exponential' || corrmodel=='gauss')
      namesparam <- c('scale')

    if(corrmodel=='gencauchy')
      namesparam <- c('power1', 'power2','scale')

    if(corrmodel=='stable')
      namesparam <- c('power', 'scale')

    if(corrmodel=='whittlematern')
      namesparam <- c('scale', 'smooth')

    return(namesparam)
  }

DetectParam <- function(corrmodel, fixed, param)
  {
    param <- c(fixed, param)
    param <- param[CorrelationParam(corrmodel)]
    corrmodel <- CheckCorrModel(corrmodel)

    return(list(corrmodel=corrmodel, param=param))
  }

InitParam <- function(coordx, coordy, corrmodel, data, fixed, grid, likelihood,
                      lonlat, model, parscale, paramrange, start, time, type)
  {    
    ### Initialize the model parameters:
    error <- NULL
    mean <- mean(data)
    nugget <- 0
    scale <- 10
    smooth <- 1
    sill <- var(data)
    numfixed <- numstart <- 0

    ### Set returning variables:
    
    codecorrmodel <- CheckCorrModel(corrmodel)
    likelihood <- CheckLikelihood(likelihood)
    model <- CheckModel(model)
    type <- CheckType(type)
    
    ### Set the names of the parameters
    
    param <- c(mean, nugget, scale, sill, rep(smooth, 4))
    namesparam <- c("mean", "nugget", "scale", "sill", "power",
                    "power1", "power2", "smooth")
    names(param) <- namesparam

    namescorr <- CorrelationParam(corrmodel)
    namesnuis <- c('mean', 'nugget', 'sill')
    namesparam <- sort(c(namescorr, namesnuis))
    namessim <- c('mean', 'sill', 'nugget', 'scale', namescorr[!namescorr == 'scale'])
    
    param <- param[namesparam]

    numparam <- length(param)
    numparamcorr <- length(namescorr)
    flag <- rep(1, numparam)
    namesflag <- namesparam
    names(flag) <- namesflag

    if(any(type == c(1, 3, 5)))
      {
        if(is.list(fixed))
          fixed$mean <- mean
        else
          fixed <- list(mean=mean)
      }

    ### Update the  parameters with fixed values:

    if(!is.null(fixed))
      {
        fixed <- unlist(fixed)
        namesfixed <- names(fixed)
        numfixed <- length(namesfixed)
    
        if(numfixed==numparam)
          {
            error <- 'the are not parameters left to estimate\n'
            return(list(error=error))
          }
        
        for(i in 1 : numfixed)
          {
            flag[namesflag==namesfixed[i]] <- 0
            param <- param[!namesparam==namesfixed[i]]
            if(any(namescorr==namesfixed[i]))
              numparamcorr <- numparamcorr - 1
            namesparam <- names(param)
          }
        numparam <- length(param)
      }
    flagcorr <- flag[namescorr]
    flagnuis <- flag[namesnuis]
    
    ### Update the parameters with starting values:
    
    if(!is.null(start))
      {
        start <- unlist(start)
        namesstart <- names(start)

        if(any(type == c(1, 3, 5)))
          if(any(namesstart == 'mean'))
            start <- start[!namesstart == 'mean']

        namesstart <- names(start)
        numstart <- length(start)
        
        for(i in 1 : numstart)
          param[namesstart[i]] <- start[namesstart[i]]
      }

    ### Check the consistency between fixed and starting values
    
    if(numstart > 0 && numfixed > 0)
      for(i in 1 : numstart)
        for(j in 1 : numfixed)
          if(namesstart[i]==namesfixed[j])
            {
              error <- ('some fixed parameter name/s is/are matching with starting parameter name/s\n')
              return(list(error=error))
            }

    ### set the scale of the parameters:

    # Insert here!
    
    ### set the range of the parameters if its the case

    if(paramrange)
      paramrange <- SetRangeParam(namesparam, numparam)
    else
      paramrange <- list(lower=NULL, upper=NULL)

    ### Set the data format:

    dimdata <- dim(data)
    numdata <- 1
    if(time)
      {
        if(grid)
          numdata <- dim(data)[3] # number of observations
        else
          numdata <- nrow(data)    # number of observations
      }

    if(!is.numeric(coordy))
      {
        coord <- coordx
        numcoord <- nrow(coord)   # number of coordinates

        if(grid)
          {
            numcoordx <- sqrt(numcoord)
            dim(data) <- c(numcoordx, numcoordx, numdata)
          }
        else
          dim(data) <- c(numdata, numcoord)              
      }
    else
      {        
        numcoordx <- length(coordx)
        numcoordy <- length(coordy)

        if(grid)
          {
            coord <- expand.grid(coordx, coordy)
            numcoord <- numcoordx * numcoordy  # number of coordinates
            dim(data) <- c(numcoordx, numcoordy, numdata)
          }
        else
          {
            coord <- cbind(coordx, coordy)
            numcoord <- numcoordx   # number of coordinates
            dim(data) <- c(numdata, numcoord)
          }
      }

    ### Compute distances:
    numpairs <- numcoord * (numcoord - 1) / 2
    lags <- double(numpairs)
    .C('Distances', as.double(coord[,1]), as.double(coord[,2]), lags, as.integer(numcoord),
       as.integer(lonlat), PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)


    return(list(corrmodel=codecorrmodel, coord=coord, data=data, error=error, flagcorr=flagcorr, flagnuis=flagnuis,
                fixed=fixed, lags=lags, likelihood=likelihood, lower=paramrange$lower, model=model, namescorr=namescorr,
                namesnuis=namesnuis, namesparam=namesparam, namessim=namessim, numcoord=numcoord, numdata=numdata,
                numpairs=numpairs, numparam=numparam, numparamcorr=numparamcorr, numfixed=numfixed, param=param,
                upper=paramrange$upper, type=type))
  }

SetRangeParam <- function(namesparam, numparam)
  {
    low <- 1e-12
    lower <- NULL
    upper <- NULL
    
    for(i in 1 : numparam)
      {
        if(namesparam[i]=='mean')
          {
            lower <- c(lower, -Inf)
            upper <- c(upper, Inf)
          }
        
        if(namesparam[i]=='nugget')
          {
            lower <- c(lower, 0)
            upper <- c(upper, Inf)
          }

        if(namesparam[i]=='power')
          {
            lower <- c(lower, low)
            upper <- c(upper, 2)
          }

        if(namesparam[i]=='power1')
          {
            lower <- c(lower, low)
            upper <- c(upper, 2)
          }

        if(namesparam[i]=='power2')
          {
            lower <- c(lower, low)
            upper <- c(upper, Inf)
          }

        if(namesparam[i]=='scale')
          {
            lower <- c(lower, low)
            upper <- c(upper, Inf)
          }

        if(namesparam[i]=='sill')
          {
            lower <- c(lower, low)
            upper <- c(upper, Inf)
          }
        
        if(namesparam[i]=='smooth')
          {
            lower <- c(lower, low)
            upper <- c(upper, Inf)
          }

      }
    return(list(lower=lower, upper=upper))
  }

SetComp <- function(type)
  {
    SetComp <- switch(type, pair=1, diff=2, ML=3, REML=4, wls=5)

    return(SetComp)
  }
