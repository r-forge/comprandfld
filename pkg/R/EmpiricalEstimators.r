####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Email: simone.padoan@unibg.it.
### Institute: University of Bergamo.
### File name: EmpiricalEstimators.r
### Description:
### This file contains a set of procedures in order
### to estimate the empircal covariance or the extreme
### dependence structures for a given dataset.
### Last change: 2011/08/03.
####################################################

### Procedures are in alphabetical order.


EVariogram <- function(coordx, coordy, data, cloud=FALSE, extcoeff=FALSE, grid=FALSE, gev=c(0, 1, 0),
                       lonlat=FALSE, maxdist=NULL, numbins=NULL, replicates=FALSE, type='variogram')
  {

    call <- match.call()

    ### Check the parameters given in input:

    checkinput <- CheckInput(coordx, coordy, 'gauss', data, NULL, grid, 'None',
                             lonlat, 'None', 'Nelder-Mead', replicates, NULL,
                             'WLeastSquare', FALSE, 'SubSamp', FALSE, NULL, NULL)

    if(!is.null(checkinput$error))
      stop(checkinput$error)

    if(!is.null(cloud) & !is.logical(cloud))
      stop('insert a logical value (TRUE/FALSE) for the cloud parameter\n')

    if(!is.null(extcoeff) & !is.logical(extcoeff))
      stop('insert a logical value (TRUE/FALSE) for the extcoeff parameter\n')

    if(extcoeff & (!is.numeric(gev) || is.null(gev) || !length(gev)==3))
      stop('insert a numeric vector with the three GEV parameters\n')

    if(!is.null(numbins) & !is.integer(numbins))
      if(numbins < 0)
        stop('insert a positive integer value for the number of bins\n')

    if(is.null(numbins))
      numbins <- 13

    if(!is.null(maxdist) & !is.numeric(maxdist))
      if(maxdist < 0)
        stop('insert a positive numeric value for maximum distance\n')

    if(is.null(maxdist))
      maxdist <- double(1)

    if(is.null(type))
      type <- 'variogram'

    if(!is.null(type) & all(type!='variogram', type!='madogram', type!='Fmadogram'))
      stop('the admitted types are: variogram or madogram\n')

    ### Initialization parameters:
    initparam <- InitParam(coordx, coordy, 'gauss', data, NULL, grid, 'None',
                           lonlat, 'None', NULL, FALSE, replicates, NULL,
                           'WLeastSquare', 'SubSamp', FALSE)

    if(!is.null(initparam$error))
      stop(initparam$error)

    numvario <- numbins - 1
    if(cloud)
      {
        numbins <- numvario <- initparam$numpairs
        maxdist <- double(1)
      }

    ### Estimation of the empirical variogram:

    bins <- double(numbins)
    moments <- double(numvario)
    lenbins <- integer(numvario)
    vtype <- as.integer(1)
    if(type=='madogram') # Trasform to GEV margins:
      {
        if(gev[3]>=1)
          stop('the shape parameter can not be greater or equal to 1')
        vtype <- as.integer(2)
        initparam$data <- Dist2Dist(initparam$data, to='Gev',
                                    loc=rep(gev[1], initparam$numcoord),
                                    scale=rep(gev[2], initparam$numcoord),
                                    shape=rep(gev[3], initparam$numcoord))
      }
    if(type=='Fmadogram') # Transform to Uniform margins:
      {
        vtype <- as.integer(2)
        initparam$data <- Dist2Dist(initparam$data, to='Uniform')
      }
    # Compute the variogram:
    .C('Empiric_Variogram', bins, as.integer(cloud), as.double(initparam$data),
       lenbins, as.double(maxdist), moments, as.integer(initparam$numdata),
       as.integer(initparam$numcoord), as.integer(numbins), vtype,
       PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)

    variogram <- moments / lenbins

    if(cloud)
      centers <- bins
    else
      centers <- bins[1:numvario] + diff(bins) / 2

    extremalcoeff <- NULL
    if(extcoeff)
      {
        if(type == 'madogram') # Check if its a madogram:
          if(gev[3] == 0)
            extremalcoeff <- exp(variogram / gev[2])
          else
            extremalcoeff <- Dist2Dist(gev[1] + variogram / gamma(1 - gev[3]))
        if(type == 'Fmadogram') # Check if its a Fmadogram:
          extremalcoeff <- (1 + 2 * variogram) / (1 - 2 * variogram)
      }

    .C('DelDistances', PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)

    EVariogram <- list(bins=bins,
                       cloud=cloud,
                       centers=centers,
                       extremalcoeff=extremalcoeff,
                       lenbins=lenbins,
                       maxdist=maxdist,
                       variogram=variogram,
                       type=type)

    structure(c(EVariogram, call = call), class = c("Variogram"))

  }

