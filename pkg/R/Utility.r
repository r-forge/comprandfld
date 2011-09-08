####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Email: simone.padoan@unibg.it.
### Institute: Department of Information Technology
### and Mathematical Methods, University of Bergamo
### File name: Utility.r
### Description:
### This file contains a set of procedures
### for supporting all the other functions.
### Last change: 08/09/2011.
####################################################

### Procedures are in alphabetical order.

CheckCorrModel <- function(corrmodel)
  {
    CheckCorrModel <- NULL
    # Correlation function are in alphabetical order
    CheckCorrModel <- switch(# spatial or temporal correlations
                             corrmodel,
                             cauchy=1,
                             exponential=2,
                             gauss=3,
                             gencauchy=4,
                             spherical=5,
                             stable=6,
                             whittlematern=7,
                             # spatial-temporal non-separable correlations
                             gneiting=21,
                             iacocesare=22,
                             porcu=23,
                             stein=24,
                             # spatial-temporal non-separable correlations
                             exp_cauchy=41,
                             exp_exp=42,
                             exp_gauss=43,
                             gauss_gauss=44,##not implemented
                             matern_cauchy=45,
                             matern_exp=46)


    return(CheckCorrModel)
  }

CheckInput <- function(coordx, coordy, coordt, corrmodel, data, fixed, grid, likelihood,
                       lonlat, margins, maxdist, maxtime, model, optimizer, replicates,
                       start, taper, threshold, type, varest, vartype, weighted, weights,
                       winconst)
  {
    error <- NULL
    # START Include internal functions:
    CheckParamRange <- function(param)
    {
        if(!is.na(param['df'])) if(param['df'] <= 0) return(FALSE)
        if(!is.na(param['nugget'])) if(param['nugget'] < 0) return(FALSE)
        if(!is.na(param['power'])) if(param['power'] <=0 || param['power'] > 2) return(FALSE)
        if(!is.na(param['power_s'])) if(param['power_s'] <=0 || param['power_s'] > 2) return(FALSE)
        if(!is.na(param['power_t'])) if(param['power_t'] <=0 || param['power_t'] > 2) return(FALSE)
        if(!is.na(param['power1'])) if(param['power1'] <=0 || param['power1'] > 2) return(FALSE)
        if(!is.na(param['power2'])) if(param['power2'] <= 0) return(FALSE)
        if(!is.na(param['sep'])) if(param['sep'] < 0 || param['sep'] > 1) return(FALSE)
        if(!is.na(param['scale'])) if(param['scale'] <= 0) return(FALSE)
        if(!is.na(param['scale_s'])) if(param['scale_s'] <= 0) return(FALSE)
        if(!is.na(param['scale_t'])) if(param['scale_t'] <= 0) return(FALSE)
        if(!is.na(param['sill'])) if(param['sill'] <= 0) return(FALSE)
        if(!is.na(param['smooth'])) if(param['smooth'] <= 0) return(FALSE)
        if(!is.na(param['smooth_s'])) if(param['smooth_s'] <= 0) return(FALSE)
        if(!is.na(param['smooth_t'])) if(param['smooth_t'] <= 0) return(FALSE)
        return(TRUE)
    }
    # Check if the correlation is spatial or spatial-temporal
    CheckSpaceTime <- function(corrmodel)
    {
        CheckSpaceTime <- NULL
        if(corrmodel >= 1 & corrmodel <= 20) CheckSpaceTime <- FALSE
        else CheckSpaceTime <- TRUE
        return(CheckSpaceTime)
    }
    # END Include internal functions
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

    if(!is.null(maxdist))
      {
        error <- "insert a positive numeric value for the maximum spatial distance\n"
        if(!is.numeric(maxdist))
          return(list(error=error))
        else
          if(maxdist<0)
            return(list(error=error))
      }

    if(!is.null(maxtime))
      {
        error <- "insert a positive numeric value for the maximum time interval\n"
        if(!is.numeric(maxtime))
          return(list(error=error))
        else
          if(maxtime<0)
            return(list(error=error))
      }

    if(!is.null(model) & !is.character(model))
      {
        error <- 'insert the name of the random field\n'
        return(list(error=error))
      }

     if(is.null(CheckModel(model)))
      {
        error <- 'the model name of the random field is not correct\n'
        return(list(error=error))
      }

    if(model=="BinaryGauss" || model=="BinaryExt"){
      # check the threshold
      error <- 'insert a numeric value for the threshold'
      if(is.null(threshold)) return(list(error=error))
      else if(!is.numeric(threshold)) return(list(error=error))
      # check the data
      if(length(unique(c(data)))!=2){
        error <- 'the data are not binary'
        return(list(error=error))}}

    if(model %in% c("BrownResn","ExtGauss","ExtT"))
        if(!margins %in% c("Frechet","Gev")){
            error <- 'insert the correct type of marginal distributions\n'
            return(list(error=error))}

    if(!is.null(optimizer) & !is.character(optimizer))
      {
        error <- 'insert the type of maximising algorithm\n'
        return(list(error=error))
      }

    if(!is.null(varest) & !is.logical(varest))
      {
        error <- 'the parameter std.err need to be a logical value\n'
        return(list(error=error))
      }

    if(is.null(replicates) || (abs(replicates-round(replicates))>0) || replicates<1)
      {
        error <- 'the parameter replicates need to be a positive integer\n'
        return(list(error=error))
      }

    if(type=="Tapering"){
        if(is.null(taper) || is.null(maxdist)){
          error <- 'tapering need a taper correlation model and/or a compact support\n'
          return(list(error=error))}
        if(!taper %in% c("Wendland1","Wendland2","Wendland3")){
           error <- 'insert a correct name for the taper correlation model\n'
           return(list(error=error))}
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
    codcorr <- CheckCorrModel(corrmodel)
    if(is.null(corrmodel))
      {
        error <- 'the name of the correlation model is not correct\n'
        return(list(error=error))
      }

    if(!is.null(fixed))
      {
        if(!all(names(fixed) %in% c(NuisanceParam(model), CorrelationParam(corrmodel))))
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

    if(!is.null(start))
      {
        if(!is.list(start))
          {
            error <- 'insert starting values as a list of parameters\n'
            return(list(error=error))
          }
        namesstart <- names(start)
        if(!all(namesstart %in% c(NuisanceParam(model), CorrelationParam(corrmodel))))
          {
            error <- 'some names of the starting parameters is/are not correct\n'
            return(list(error=error))
          }
        if(any(namesstart=='mean') & (type=='Difference' || type=='Restricted'))
          {
            error <- 'the mean parameter is not allow with the difference composite likelihood\n'
            return(list(error=error))
          }
        if(any(namesstart=='sill') & (model=='BrowResn'))
          {
            error <- 'the sill parameter is not allow with Brown-Renick model\n'
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
        if(!any(checktype == c(3, 4, 5)))
          {
            error <- 'insert a type name of the likelihood objects compatible with the full likelihood\n'
            return(list(error=error))
          }
      }

    if(checklik == 3)
      {
        if(!any(checktype == c(1, 2)))
          {
            error <- 'insert a type name of the likelihood objects compatible with the composite-likelihood\n'
            return(list(error=error))
          }
      }

    if(varest & (checklik == 1 || checklik == 3)  & (!is.null(vartype) & !is.character(vartype)))
      {
        error <- 'insert the type of estimation method for the variances\n'
        return(list(error=error))
      }

    vartp <- CheckVarType(vartype)

    if(varest & is.null(vartp) & (checklik == 1 || checklik == 3))
      {
        error <- 'the name of the estimation type for the variances is not correct\n'
        return(list(error=error))
      }

    # START - check the format of the inserted coordinates and dataset
    dimdata <- dim(data) # set the data dimension
    if(is.null(coordt)) # START 1) spatial random field
      {
        if(CheckSpaceTime(codcorr))
          {
            error <- 'temporal coordinates are missing\n'
            return(list(error=error))
          }
        if(replicates>1) # a) n iid replicates of a spatial random field
          {
            if(grid) # START regular grid
              {
                if(is.null(dimdata))
                  {
                    error <- c('insert an array d x d of n iid spatial observations\n')
                    return(list(error=error))
                  }
                if(length(dimdata)!=3)
                  {
                    error <- c('the dimension of the data matrix is not correct\n')
                    return(list(error=error))
                  }
                if(length(coordx)!=dimdata[1] || length(coordy)!=dimdata[2])
                  {
                    error <- c('the number of coordinates does not match with the number of spatial observations\n')
                    return(list(error=error))
                  }
                if(dimdata[3]!=replicates)
                  {
                    error <- c('the number of replications does not match with the data observations\n')
                    return(list(error=error))
                  }
              } # END regular grid
            else # START irregular grid
              {
                if(is.null(dimdata))
                  {
                    error <- c('insert a matrix n x d of spatial observations\n')
                    return(list(error=error))
                  }
                if(length(dimdata)!=2)
                  {
                    error <- c('the dimension of the data matrix is not correct\n')
                    return(list(error=error))
                  }
                if(is.null(coordy))
                  {
                    if(is.null(dim(coordx)))
                      {
                        error <- c('insert a matrix d x 2 of spatial coordinates\n')
                        return(list(error=error))
                      }
                    if(dimdata[2]!=nrow(coordx) || ncol(coordx)!=2)
                      {
                        error <- c('the number of coordinates does not match with the number of spatial observations\n')
                        return(list(error=error))
                      }
                  }
                else
                  if(length(coordx)!=length(coordy))
                    {
                      error <- c('the number of the two coordinates does not match\n')
                      return(list(error=error))
                    }
                if(dimdata[1]!=replicates)
                  {
                    error <- c('the number of replications does not match with the data observations\n')
                    return(list(error=error))
                  }
              } # END irregular grid
          }
        else # b) START one realisation of a spatial random field
          {
            if(grid) # START regular grid
              {
                if(is.null(dimdata))
                  {
                    error <- c('insert a matrix d x d of spatial observations\n')
                    return(list(error=error))
                  }
                if(length(dimdata)!=2)
                  {
                    error <- c('the dimension of the data matrix is not correct\n')
                    return(list(error=error))
                  }
                if(length(coordx)!=dimdata[1] || length(coordy)!=dimdata[2])
                  {
                    error <- c('the number of coordinates does not match with the number of spatial observations\n')
                    return(list(error=error))
                  }
              } # END regular grid
            else # START irregular grid
              {
                numsite <- length(data)
                if(is.null(numsite))
                  {
                    error <- c('insert a vector of spatial observations\n')
                    return(list(error=error))
                  }
                if(is.null(coordy))
                  {
                    dimcoord <- dim(coordx)
                    if(is.null(dimcoord))
                      {
                        error <- c('insert a suitable set of coordinates\n')
                        return(list(error=error))
                      }
                    else
                      {
                        if(dimcoord[1]!=numsite || dimcoord[2]!=2)
                          {
                            error <- c('the number of coordinates does not match with the number of spatial observations\n')
                            return(list(error=error))
                          }
                      }
                  }
                else
                  {
                    if(length(coordx)!=length(coordy))
                      {
                        error <- c('the number of the two coordinates does not match\n')
                        return(list(error=error))
                      }
                    if(length(coordx)!=numsite)
                      {
                        error <- c('the number of coordinates does not match with the number of spatial observations\n')
                        return(list(error=error))
                      }
                  }
              }
          } # b) END one realisation of a spatial random field
      } # END 1) spatial random field
    else # 2) case: spatial-temporal random field
      {
        if(!is.numeric(coordt))
          {
            error <- 'insert a numerical vector of temporal coordinates\n'
            return(list(error=error))
          }
        if(length(coordt)<=1)
          {
            error <- 'insert a numerical vector of temporal coordinates\n'
            return(list(error=error))
          }
        if(replicates>1) # START a) n iid replicates of a spatial-temporal random field
          {
            if(grid) # START regular grid
              {
                if(is.null(dimdata))
                  {
                    error <- c('insert an array d x d x t of n iid spatial-temporal observations\n')
                    return(list(error=error))
                  }
                if(length(dimdata)!=4)
                  {
                    error <- c('the dimension of the data matrix is not correct\n')
                    return(list(error=error))
                  }
                if(length(coordx)!=dimdata[1] || length(coordy)!=dimdata[2] )
                  {
                    error <- c('the number of coordinates does not match with the number of spatial observations\n')
                    return(list(error=error))
                  }
                if(length(coordt)!=dimdata[3])
                  {
                    error <- c('the number of the temporal coordinate does not match with the third dimensiona of the data matrix\n')
                    return(list(error=error))
                  }
                if(dimdata[4]!=replicates)
                  {
                    error <- c('the number of replications doen not match with the data observations\n')
                    return(list(error=error))
                  }
              } # END regular grid
            else # START irregular grid
              {
                if(is.null(dimdata))
                  {
                    error <- c('insert an array t x d of n iid spatial-temporal observations\n')
                    return(list(error=error))
                  }
                if(length(dimdata)!=3)
                  {
                    error <- c('the dimension of the data matrix is not correct\n')
                    return(list(error=error))
                  }
                if(is.null(coordy))
                  {
                    if(dimdata[2]!=nrow(coordx) || ncol(coordx)!=2)
                      {
                        error <- c('the number of coordinates does not match with the number of spatial observations\n')
                        return(list(error=error))
                      }
                  }
                else
                  {
                    if(length(coordx)!=length(coordy))
                      {
                        error <- c('the number of the two spatial coordinates does not match\n')
                        return(list(error=error))
                      }
                    if(length(coordx)!=dimdata[2])
                      {
                        error <- c('the number of coordinates does not match with the number of spatial observations\n')
                        return(list(error=error))
                      }
                  }
                if(dimdata[1]!=length(coordt))
                  {
                    error <- c('the time coordinate does not match with the number of rows of the data array\n')
                    return(list(error=error))
                  }
                if(dimdata[3]!=replicates)
                  {
                    error <- c('the number of replications does not match with the data observations\n')
                    return(list(error=error))
                  }
              } # END irregular grid
          }# END a) n iid replicates of a spatial-temporal random field
        else # START b) one realisation of a spatial-temporal random field
          {
            if(grid) # START regular grid
              {
                if(is.null(dimdata))
                  {
                    error <- c('insert an array of d x d x t spatial-temporal observations\n')
                    return(list(error=error))
                  }
                if(length(dimdata)!=3)
                  {
                    error <- c('the dimension of the data matrix is not correct\n')
                    return(list(error=error))
                  }
                if(length(coordx)!=dimdata[1] || length(coordy)!=dimdata[2])
                  {
                    error <- c('the number of coordinates does not match with the number of spatial observations\n')
                    return(list(error=error))
                  }
                if(dimdata[3]!=length(coordt))
                  {
                    error <- c('the time coordinate does not match with the third dimension of the data array\n')
                    return(list(error=error))
                  }
              } # END regular grid
            else # START irregular grid
              {
                if(is.null(dimdata))
                  {
                    error <- c('insert a matrix of t x d spatial-temporal observations\n')
                    return(list(error=error))
                  }
                if(length(dimdata)!=2)
                  {
                    error <- c('the dimension of the data matrix is not correct\n')
                    return(list(error=error))
                  }
                if(is.null(coordy))
                  {
                    if(dimdata[2]!=nrow(coordx) || ncol(coordx)!=2)
                      {
                        error <- c('the number of coordinates does not match with the number of spatial observations\n')
                        return(list(error=error))
                      }
                  }
                else
                  {
                    if(length(coordx)!=length(coordy))
                      {
                        error <- c('the number of the two spatial coordinates does not match\n')
                        return(list(error=error))
                      }
                    if(length(coordx)!=dimdata[2])
                      {
                        error <- c('the number of coordinates does not match with the number of spatial observations\n')
                        return(list(error=error))
                      }
                  }
                if(dimdata[1]!=length(coordt))
                  {
                    error <- c('the time coordinate does not match with the number of the matrix rows\n')
                    return(list(error=error))
                  }
              } # END irregular grid
          } # END b) one realisation of a spatial-temporal random field
      }
    # END - check the format of the inserted coordinates and dataset
    # Check the range of validity for the sub-sampling parameter:
    if(!CheckSpaceTime(codcorr))
    if(varest & (vartp==2) & is.numeric(winconst))
      {
        if(!grid)
          {
            rcoordx <- range(coordx[, 1])
            rcoordy <- range(coordx[, 2])
          }
        else
          {
            rcoordx <- range(coordx)
            rcoordy <- range(coordy)
          }
         delta <- min(rcoordx[2] - rcoordx[1], rcoordy[2] - rcoordy[1])
         wincup <- delta / sqrt(delta)

         if(winconst < 0 || winconst > wincup)
           {
             error <- paste('for the sub-sampling constant insert a positive real value less or equal than:', wincup,'\n')
             return(list(error=error))
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
                         Gaussian=1,
                         BinaryGauss=2,
                         BrowResn=3,
                         ExtGauss=4,
                         ExtT=5)
    return(CheckModel)
  }

CheckVarType <- function(type)
  {
    CheckVarType <- switch(type,
                           Sampling=1,
                           SubSamp=2,
                           Theoretical=3)
    return(CheckVarType)
  }

CheckType <- function(type)
  {
    CheckType <- switch(type,
                        Difference=1,
                        Pairwise=2,
                        Restricted=3,
                        Standard=4,
                        Tapering=5,
                        WLeastSquare=6)
    return(CheckType)
  }

CorrelationParam <- function(corrmodel)
  {
    param <- NULL
    # Cauchy correlation model:
    if(corrmodel=='cauchy'){
      param <- c('power2', 'scale')
      return(param)}
    # Exponential and Gaussian correlation models:
    if(corrmodel=='exponential' || corrmodel=='gauss' || corrmodel=='spherical'){
      param <- c('scale')
      return(param)}
    # Generalised Cauchy correlation model:
    if(corrmodel=='gencauchy'){
      param <- c('power1', 'power2','scale')
      return(param)}
    # Stable correlation model:
    if(corrmodel=='stable'){
      param <- c('power', 'scale')
      return(param)}
    # Whittle-Matern correlation model:
    if(corrmodel=='whittlematern'){
      param <- c('scale', 'smooth')
      return(param)}
    # Non-separable spatial-temporal correlations:
    # Gneiting model:
    if(corrmodel=='gneiting'){
      param <- c('power_s', 'power_t','scale_s','scale_t','sep')
      return(param)}
    # Iaco-Cesare model:
    if(corrmodel=='iacocesare'){
      param <- c('power2','power_s', 'power_t','scale_s','scale_t')
      return(param)}
    # Porcu model:
    if(corrmodel=='porcu'){
      param <- c('power_s', 'power_t','scale_s','scale_t','sep')
      return(param)}
    # Stein model:
    if(corrmodel=='stein'){
      namesparam <- c('power_t','scale_s','scale_t','smooth_s')
      return(param)}
    # Separable spatial-temporal correlations:
    # Exponential-exponential and exponential-Gaussian models:
    if(corrmodel=='exp_exp'||corrmodel=='exp_gauss'){
      param <- c('scale_s','scale_t')
      return(param)}
    # Exponential-Cauchy model:
     if(corrmodel=='exp_cauchy'){
       param <- c('power2','scale_s','scale_t')
       return(param)}
    # Whittle-Matern-exponential model:
     if(corrmodel=='matern_exp'){
       param <- c('scale_s','scale_t','smooth_s')
       return(param)}
    # Whittle-Matern-Cauchy model:
     if(corrmodel=='matern_cauchy'){
       param <- c('power2','scale_s','scale_t','smooth_s')
       return(param)}

    return(param)
  }

NuisanceParam <- function(model)
{
  param <- NULL
  # Gaussian random field:
  if(model=='Gaussian' || model=='BinaryGauss'){
    param <- c('mean', 'nugget', 'sill')
    return(param)}
  # Max-stable random field (Extremal Gaussian):
  if(model=='ExtGauss'){
    param <- c('sill')
    return(param)}
  # Max-stable random field (Brown Resnick):
  if(model=='BrowResn'){
    param <- c('sill')
    return(param)}
  # Max-stable random field (Extremal T):
  if(model=='ExtT'){
    param <- c('df', 'sill')
    return(param)}
  return(param)
}

InitParam <- function(coordx, coordy, coordt, corrmodel, data, fixed, grid, likelihood,
                      lonlat, margins, maxdist, maxtime, model, parscale, paramrange,
                      replicates, start, threshold, type, varest, vartype, weighted,
                      winconst)
{
    ### START Includes internal functions:
    # Check if the correlation is spatial or spatial-temporal
    CheckSpaceTime <- function(corrmodel)
    {
        CheckSpaceTime <- NULL
        if(corrmodel >= 1 & corrmodel <= 20) CheckSpaceTime <- FALSE
        else CheckSpaceTime <- TRUE
        return(CheckSpaceTime)
    }
    # Determines the range of the parameters for a given correlation
    SetRangeParam <- function(namesparam, numparam)
    {
        low <- 1e-12
        lower <- NULL
        upper <- NULL
        # Check for the param set:
        for(i in 1:numparam){
            if(namesparam[i]=='mean'){
                lower <- c(lower, -Inf)
                upper <- c(upper, Inf)}
            if(namesparam[i]=='nugget'){
                lower <- c(lower, 0)
                upper <- c(upper, Inf)}
            if(namesparam[i]=='power'){
                lower <- c(lower, low)
                upper <- c(upper, 2)}
            if(namesparam[i]=='power1'){
                lower <- c(lower, low)
                upper <- c(upper, 2)}
            if(namesparam[i]=='power2'){
                lower <- c(lower, low)
                upper <- c(upper, Inf)}
            if(namesparam[i]=='scale'){
                lower <- c(lower, low)
                upper <- c(upper, Inf)}
            if(namesparam[i]=='sill'){
                lower <- c(lower, low)
                upper <- c(upper, Inf)}
            if(namesparam[i]=='smooth'){
                lower <- c(lower, low)
                upper <- c(upper, Inf)}}
        return(list(lower=lower, upper=upper))
    }
    ### END Includes internal functions
    ### Set returning variables and initialize the model parameters:
    namesnuis <- NuisanceParam(model)
    likelihood <- CheckLikelihood(likelihood)
    model <- CheckModel(model)
    vartype <- CheckVarType(vartype)
    type <- CheckType(type)
    # Initialises the starting and fixed parameters' names
    error <- NULL
    namesfixed <- NULL
    namesstart <- NULL
    numfixed <- numstart <- 0
    # Set the correlation and the nuisance parameters:
    namescorr <- CorrelationParam(corrmodel)
    numparamcorr <- length(namescorr)
    paramcorr <- rep(1, numparamcorr)
    names(paramcorr) <- namescorr
    nuisance <- NULL
    # Set the if the correlation is space-time:
    corrmodel <- CheckCorrModel(corrmodel)
    spacetime <- CheckSpaceTime(corrmodel)
    ### Parameters' settings:
    if(model==1){ # Gaussian or Binary Gaussian random field:
        mu <- mean(data)
        if(any(type==c(1, 3, 6)))# Checks the type of likelihood
          if(is.list(fixed)) fixed$mean <- mu# Fixs the mean
          else fixed <- list(mean=mu)
        nuisance <- c(mu, 0, var(c(data)))}
    if(model==2){ # Binary Gaussian random field:
        p <- mean(data)
        mu <- threshold+qnorm(p)
        nuisance <- c(mu, 0, 1)
        if(!is.null(start$nugget))
            if(length(start)>1) start<-start[!names(start)%in%"nugget"]
            else start<-NULL
        if(is.list(fixed)) fixed$nugget<-0# Fixs the nugget
          else fixed<-list(nugget=0)}
    if(model>2){ # Max-stable random field:
        if(model==3)# Checks if its the Brown-Resnick model
          if(is.list(fixed)) fixed$sill <- 1 # Fixs the sill
          else fixed <- list(sill=1)
        # Cheks if its the Extremal-t model
        if(model==5) nuisance <- c(nuisance,1)
        nuisance <- c(nuisance,0.5)
        if(margins=="Gev") data <- Dist2Dist(data)}
    names(nuisance) <- namesnuis
    namesparam <- sort(c(namescorr, namesnuis))
    param <- c(nuisance, paramcorr)
    param <- param[namesparam]
    namessim <- c('mean', 'sill', 'nugget', 'scale', namescorr[!namescorr == 'scale'])
    numparam <- length(param)
    flag <- rep(1, numparam)
    namesflag <- namesparam
    names(flag) <- namesflag
    ### Update the parameters with fixed values:
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

        for(i in 1:numfixed)
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

    ### Check the consistency between fixed and starting values:
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
    ### If the case set the sub-sampling parameter to the default value
    if(varest & (vartype==2) & (missing(winconst) || !is.numeric(winconst)))
      if(spacetime) winconst <- 0
      else winconst <- 1
    ### START settings the data structure:
    # set the coordinates sizes:
    if(is.null(coordy)){coordy <- coordx[,2]
                        coordx <- coordx[,1]}
    # checks is the data are on regular grid:
    if(grid) {numcoordx <- length(coordx)
              numcoordy <- length(coordy)
              numcoord <- numcoordx*numcoordy}
    else numcoord <- numcoordx <- numcoordy <- length(coordx)

    if(spacetime)
      { # set the number of temporal realisations:
        numtime <- length(coordt)
        if(grid) # if the data are in regular grid:
          if(replicates>1){ # checks if there are iid replicates:
              dim(data) <- c(numcoord, numtime, replicates)
              data <- aperm(data, c(2,1,3))}
          else # if there are not iid replicates:
            data <- matrix(data, ncol=numcoord, nrow=numtime, byrow=T)}
    else{numtime <- 1
         coordt <- 0
         if(grid)
           data <- matrix(data, ncol=numcoord, nrow=replicates, byrow=T)
         else
           data <- matrix(data, ncol=numcoord, nrow=replicates)}
    ### END settings the data structure

    ### Compute distances:
    numpairs=.5*numcoord*(numcoord-1)
    srange <- double(1)
    trange <- double(1)
    if(is.null(maxdist)) srange <- c(srange, double(1)) else srange <- c(srange, as.double(maxdist))
    if(is.null(maxtime)) trange <- c(trange, double(1)) else trange <- c(trange, as.double(maxtime))
    isinit <- as.integer(1)
    .C('SetGlobalVar', as.double(coordx), as.double(coordy), as.double(coordt),
       as.integer(grid), isinit, as.integer(numcoord), as.integer(numcoordx),
       as.integer(numcoordy), as.integer(replicates), as.integer(spacetime),
       srange, as.integer(numtime), trange, as.integer(lonlat), as.integer(weighted),
       PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)
    ### Returned list of objects:
    return(list(coordx=coordx, coordy=coordy, coordt=coordt, corrmodel=corrmodel, data=data,
                error=error, flagcorr=flagcorr, flagnuis=flagnuis, fixed=fixed, likelihood=likelihood,
                lower=paramrange$lower, model=model, namescorr=namescorr, namesfixed=namesfixed,
                namesnuis=namesnuis, namesparam=namesparam, namessim=namessim, namesstart=namesstart,
                numcoord=numcoord, numfixed=numfixed, numpairs=numpairs, numparam=numparam,
                numparamcorr=numparamcorr, numrep=replicates, numstart=numstart, numtime=numtime,
                param=param, spacetime=spacetime, srange=srange, start=start, upper=paramrange$upper,
                type=type, threshold=threshold, trange=trange, vartype=vartype, winconst=winconst))
  }
