####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Emails: simone.padoan@stat.unipd.it,
### moreno.bevilacqua@unibg.it
### Institutions: Department of Statistical Science,
### University of Padua and Department of Information
### Technology and Mathematical Methods, University
### of Bergamo.
### File name: Simulation.r
### Description:
### This file contains a set of procedures
### for the simulation of Gaussian random fields and
### related functions.
### Last change: 26/03/2012.
####################################################


### Procedures are in alphabetical order.

# Simulate spatial and spatio-temporal random felds:
RFsim <- function(coordx, coordy=NULL, coordt=NULL, corrmodel, grid=FALSE,
                  model='Gaussian', numblock=NULL, param, replicates=1, threshold=NULL)
{
    ### START -- exact direct simulation based on the Cholesky decomposition
    CholRFsim<- function(corrmodel, grid,  model, nuisance, numcoord, numcoordx,
                         numcoordy, numtime, paramcorr, replicates, spacetime,
                         threshold)
    {
        # Compute the correlation function:
        corr <- double(.5*(numcoord*numtime)*(numcoord*numtime-1))
        if(!spacetime) fname <- "CorrelationMat" else fname <- "CorrelationMat_st"
        .C(fname, corr, as.integer(corrmodel), as.double(nuisance), as.double(paramcorr),
           PACKAGE='CompRandFld', DUP=FALSE, NAOK=TRUE)
        # Builds the covariance matrix:
        varcov <- (nuisance['nugget'] + nuisance['sill']) * diag(numcoord*numtime)
        corr <- corr * nuisance['sill']
        varcov[lower.tri(varcov)] <- corr
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- corr
        # Simulation based on the Cholesky decomposition:
        if(grid){
            sim <- array(double(numcoord*numtime*replicates), c(numcoordx, numcoordy, numtime, replicates))
            for(i in 1:replicates)
                sim[,,,i] <- array(nuisance['mean']+t(chol(varcov))%*%rnorm(numcoord*numtime),
                                   c(numcoordx,numcoordy, numtime))
            if(replicates==1)
                if(!spacetime) sim <- array(sim, c(numcoordx, numcoordy))
                else sim <- array(sim, c(numcoordx, numcoordy, numtime))
            else if(!spacetime) sim <- array(sim, c(numcoordx, numcoordy, replicates))}
        else{
            sim <- array(double(numcoord*numtime*replicates), c(numtime, numcoord, replicates))
            for(i in 1:replicates)
                sim[,,i] <- nuisance['mean']+t(chol(varcov))%*%rnorm(numcoord*numtime)
            if(replicates==1)
                if(!spacetime) sim <- c(sim)
                else sim <- matrix(sim, nrow=numtime, ncol=numcoord)
            else if(!spacetime) sim <- matrix(sim, nrow=replicates, ncol=numcoord, byrow=TRUE)}
        # If a Binary-Gaussian field is required, the Gaussian field is transformed into a Binary version:
        if(model==2){
            simdim <- dim(sim)
            sim <- as.numeric(sim>threshold)
            dim(sim) <- simdim}
      return(sim)
    }
    ### END -- exact direct simulation based on the Cholesky decomposition
    # Check the user input
    checkinput <- CheckInput(coordx, coordy, coordt, corrmodel, NULL, "Simulation",
                             NULL, grid, NULL, NULL, NULL, NULL, NULL, model,
                             numblock, NULL, param, replicates, NULL, NULL,
                             threshold, NULL, NULL, NULL, NULL)

    if(!is.null(checkinput$error)) stop(checkinput$error)
    # Initialising the simulation parameters:
    initparam <- InitParam(coordx, coordy, coordt, corrmodel, NULL, "Simulation",
                           NULL, grid, NULL, FALSE, NULL, NULL, NULL, model, numblock,
                           param, NULL, NULL, replicates, NULL, threshold, NULL, NULL,
                           NULL, NULL, FALSE, NULL, NULL)
    if(!is.null(initparam$error)) stop(initparam$error)
    ### Simulation of Gaussian or Binary-Gaussian random fields:
    if(model %in% c("Gaussian","BinaryGauss")){
        # Direct method based on Cholesky decomposition:
        sim <- CholRFsim(initparam$corrmodel,grid,initparam$model,initparam$param[initparam$namesnuis],
                         initparam$numcoord,initparam$numcoordx,initparam$numcoordy,initparam$numtime,
                         initparam$param[initparam$namescorr],initparam$numrep,initparam$spacetime,
                         initparam$threshold)}
    if(model=="ExtGauss"){
        # Simulate directly Extremal Gaussian random fields from Random Felds:
        sim <- RandomFields::MaxStableRF(x=initparam$coordx, y=initparam$coordy, model=corrmodel, grid=grid,
                                         maxstable="extr",param=initparam$param[initparam$namessim],
                                         n=replicates,pch='')
        # Formatting of output:
        if(grid){
            if(replicates==1) sim <- array(sim, c(initparam$numcoordx, initparam$numcoordy))
            else sim <- array(sim, c(initparam$numcoordx, initparam$numcoordy, replicates))}
        else{
            if(replicates==1) sim <- c(sim)
            else sim <- matrix(sim, nrow=replicates,ncol=initparam$numcoord, byrow=TRUE)}}
    if(model=="ExtT" || model=="BrowResn"){
        # Set simulation variables:
        sim <- NULL
        RandomFields::DeleteAllRegisters()
        RandomFields::RFparameters(Storing=TRUE, PrintLevel=1)
        for(i in 1:replicates){
            maxima <- double(initparam$numcoord)
            # Simulate underlying Gaussian random fields (to get then Student-t fields):
            onesim <- RandomFields::GaussRF(x=initparam$coordx, y=initparam$coordy,model=corrmodel,
                                            param=initparam$param[initparam$namessim],grid=grid,
                                            n=initparam$numblock,pch='')
            # Compute compotentwise maxima:
              .C("ComputeMaxima",as.double(initparam$param["df"]),maxima,
                 as.integer(initparam$model),as.integer(initparam$numblock),
                 as.integer(initparam$numcoord),onesim,
                 PACKAGE='CompRandFld',DUP=FALSE,NAOK=TRUE)
            # Update:
            sim <- c(sim,maxima)}
        # Formatting of output:
        if(grid){
            if(replicates==1) sim <- array(sim, c(initparam$numcoordx, initparam$numcoordy))
            else sim <- array(sim, c(initparam$numcoordx, initparam$numcoordy, replicates))}
        else{
            if(replicates==1) sim <- c(sim)
            else sim <- matrix(sim, nrow=replicates,ncol=initparam$numcoord, byrow=TRUE)}}
    # Delete the global variables:
    .C('DeleteGlobalVar', PACKAGE='CompRandFld', DUP = FALSE, NAOK=TRUE)
    # Return the objects list:
    RFsim <- list(coordx = initparam$coordx,
                  coordy = initparam$coordy,
                  coordt = initparam$coordt,
                  corrmodel = corrmodel,
                  data = sim,
                  grid = grid,
                  model = model,
                  numcoord = initparam$numcoord,
                  numtime = initparam$numtime,
                  param = initparam$param,
                  randseed=.Random.seed,
                  replicates = initparam$replicates,
                  spacetime = initparam$spacetime,
                  threshold = initparam$threshold)

    structure(c(RFsim, call = call), class = c("RFsim"))
}

