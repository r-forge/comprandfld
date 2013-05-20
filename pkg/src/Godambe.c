/*###################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Emails: simone.padoan@stat.unipd.it,
### moreno.bevilacqua@unibg.it
### Institutions: Department of Statistical Science,
### University of Padua and Department of Information
### Technology and Mathematical Methods, University
### of Bergamo.
### File name: Godambe.c
### Description:
### This file contains a set of procedures
### for the computation of the Godambe matrix of
### random fields.
### Last change: 12/06/2012.
##################################################*/

#include "header.h"
// Empirical Godambe matrix in the conditional spatial Gaussian case:
void God_Cond_Gauss(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
		    int *npar, int *nparc, double *parcor, double *nuis, double *score,
		    double *sensmat, double *varimat)
{
  int d=0,i=0,j=0,h=0,k=0,m=0,n=0;
  double *gradcor, *grad, *gradient,rho=0.0;
  // Set the gradient vectors:
  gradcor=(double *) R_alloc(*nparc,sizeof(double));
  grad=(double *) R_alloc(*npar,sizeof(double));
  gradient=(double *) R_alloc(*npar,sizeof(double));

  for(n=0;n<*nrep;n++){
    // Initialize the gradient vector
    for(i=0;i<*npar;i++) gradient[i]=0;
    // Compute the gradient of a log likelihood object
    // Main loop, respect the pairs:
   for(i=0; i<(*ncoord-1);i++){
    for(j=(i+1); j<*ncoord;j++){
         if(mlags[i][j]<*maxdist){
      //Compute the correlation function
      rho=CorFct(cormod,mlags[i][j],0,parcor);
      //Compute the gradient for a given correlation model
	  GradCorrFct(rho,cormod,eps,flagcor,gradcor,mlags[i][j],0,parcor);
	  // Gradient of the log conditional likelihood
	  Grad_Cond_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,
			  data[(n+i* *nrep)], data[(n+j* *nrep)]);
	  m=0;
	  // Set the sensitivity matrix:
	  for(d=0;d<*npar;d++){
	    gradient[d]=gradient[d]+grad[d];
	    score[d]=score[d]+grad[d];
	    for(k=d;k<*npar;k++){
	      sensmat[m]=sensmat[m]+grad[d]*grad[k];
	      m++;}}}}}
    m=0;
    // Set the variability matrix:
    for(h=0;h<*npar;h++)
      for(k=h;k<*npar;k++){
	varimat[m]=varimat[m]+gradient[h]*gradient[k];
	m++;}}
  return;
}
// Empirical Godambe matrix in the conditional spatio-temporal Gaussian case:
void God_Cond_Gauss_st(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
		       int *npar, int *nparc, double *parcor, double *nuis, double *score,
		       double *sensmat, double *varimat)
{
  int d=0,h=0,i=0,j=0,k=0,m=0,n=0,t=0,v=0;
  double *gradcor, *grad, *gradient,rho=0.0;

  gradcor=(double *) R_alloc(*nparc,sizeof(double));
  grad=(double *) R_alloc(*npar,sizeof(double));
  gradient=(double *) R_alloc(*npar,sizeof(double));

  for(n=0;n<*nrep;n++){
    // Initialize the gradient vector
    for(i=0;i<*npar;i++) gradient[i]=0;
    //Main loop, respect the pairs:
 for(i=0;i<*ncoord;i++){
   for(t=0;t<*ntime;t++){
     for(j=i;j<*ncoord;j++){
       if(i==j){// marginal temporal log-likelihood:
	 for(v=t+1;v<*ntime;v++){
	   if(mlagt[t][v]<=*maxtime){

      //Compute the temporal correlation function for a given pair:
      rho=CorFct(cormod,0,mlagt[t][v],parcor);
      //Compute the gradient of the temporal correlation function for a given pair:
      GradCorrFct(rho,cormod,eps,flagcor,gradcor,0,mlagt[t][v],parcor);
      //Compute the gradient of the conditional likelihood:
      Grad_Cond_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,
		      data[(t+*ntime*i)+n* *ntime* *ncoord],
		      data[(v+*ntime*i)+n* *ntime* *ncoord]);
      m=0;
      //Set the sensitivity matrix:
      for(d=0;d<*npar;d++){
	gradient[d]=gradient[d]+grad[d];
	score[d]=score[d]+grad[d];
	for(k=d;k<*npar;k++){
	  sensmat[m]=sensmat[m]+grad[d]*grad[k];
	  m++;}}}}}
	  else {
	   for(v=0;v<*ntime;v++){// spatial-temporal log-likelihood:
	   if(mlagt[t][v]<=*maxtime && mlags[i][j]<=*maxdist){
	        //Compute the temporal correlation function for a given pair:
      rho=CorFct(cormod,mlags[i][j],mlagt[t][v],parcor);
      //Compute the gradient of the temporal correlation function for a given pair:
      GradCorrFct(rho,cormod,eps,flagcor,gradcor,mlags[i][j],mlagt[t][v],parcor);
      //Compute the gradient of the conditional likelihood:
      Grad_Cond_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,
		      data[(t+*ntime*i)+n* *ntime* *ncoord],
		      data[(v+*ntime*j)+n* *ntime* *ncoord]);
      m=0;
      //Set the sensitivity matrix:
      for(d=0;d<*npar;d++){
	gradient[d]=gradient[d]+grad[d];
	score[d]=score[d]+grad[d];
	for(k=d;k<*npar;k++){
	  sensmat[m]=sensmat[m]+grad[d]*grad[k];
	  m++;}}}}}}}}



    m=0;
    //Set the variability matrix:
    for(h=0;h<*npar;h++)
      for(k=h;k<*npar;k++){
	varimat[m]=varimat[m]+gradient[h]*gradient[k];
	m++;}}
  return;
}
// Empirical Godambe matrix in the difference Gaussian case:
void God_Diff_Gauss(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
		    int *npar, int *nparc, double *parcor, double *nuis, double *score,
		    double *sensmat, double *varimat)
{
  int d=0,i=0,j=0,h=0,k=0,m=0,n=0;
  double *gradcor, *grad, *gradient,rho=0.0;
  // Set the gradient vectors:
  gradcor=(double *) R_alloc(*nparc,sizeof(double));
  grad=(double *) R_alloc(*npar,sizeof(double));
  gradient=(double *) R_alloc(*npar,sizeof(double));

  for(n=0;n<*nrep;n++){
    // Initialize the gradient vector
    for(i=0;i<*npar;i++) gradient[i]=0;
    // Compute the gradient of a log likelihood object
    // Main loop, respect the pairs:
         for(i=0; i<(*ncoord-1);i++){
    for(j=(i+1); j<*ncoord;j++){
         if(mlags[i][j]<*maxdist){
      //Compute the correlation function
      rho=CorFct(cormod,mlags[i][j],0,parcor);
      //Compute the gradient for a given correlation model
      GradCorrFct(rho,cormod,eps,flagcor,gradcor,mlags[i][j],0,parcor);
      //Gradient of the log difference likelihood
      Grad_Diff_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,
		      data[(n+i* *nrep)], data[(n+j* *nrep)]);
      m=0;
      //Set the sensitivity matrix:
      for(d=0;d<*npar;d++){
	gradient[d]=gradient[d]+grad[d];
	score[d]=score[d]+grad[d];
	for(k=d;k<*npar;k++){
	  sensmat[m]=sensmat[m]+grad[d]*grad[k];
	  m++;}}}}}
    m=0;
    //Set the variability matrix:

    for(h=0;h<*npar;h++)
      for(k=h;k<*npar;k++){
	varimat[m]=varimat[m]+gradient[h]*gradient[k];
	m++;}}
  return;
}
// Empirical Godambe matrix in the difference spatio-temporal Gaussian case:
void God_Diff_Gauss_st(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
		       int *npar, int *nparc, double *parcor, double *nuis, double *score,
		       double *sensmat, double *varimat)
{
  int d=0,h=0,i=0,j=0,k=0,m=0,n=0,t=0,v=0;
  double *gradcor, *grad, *gradient,rho=0.0;

  gradcor=(double *) R_alloc(*nparc,sizeof(double));
  grad=(double *) R_alloc(*npar,sizeof(double));
  gradient=(double *) R_alloc(*npar,sizeof(double));

  for(n=0;n<*nrep;n++){
    // Initialize the gradient vector
    for(i=0;i<*npar;i++) gradient[i]=0;
    //Main loop, respect the pairs:
    for(i=0;i<*ncoord;i++){
   for(t=0;t<*ntime;t++){
     for(j=i;j<*ncoord;j++){
       if(i==j){// marginal temporal log-likelihood:
	 for(v=t+1;v<*ntime;v++){
	   if(mlagt[t][v]<=*maxtime){

      //Compute the temporal correlation function for a given pair:
      rho=CorFct(cormod,0,mlagt[t][v],parcor);
      //Compute the gradient of the temporal correlation function for a given pair:
      GradCorrFct(rho,cormod,eps,flagcor,gradcor,0,mlagt[t][v],parcor);
      //Compute the gradient of the conditional likelihood:
      Grad_Diff_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,
		      data[(t+*ntime*i)+n* *ntime* *ncoord],
		      data[(v+*ntime*i)+n* *ntime* *ncoord]);
      m=0;
      //Set the sensitivity matrix:
      for(d=0;d<*npar;d++){
	gradient[d]=gradient[d]+grad[d];
	score[d]=score[d]+grad[d];
	for(k=d;k<*npar;k++){
	  sensmat[m]=sensmat[m]+grad[d]*grad[k];
	  m++;}}}}}
else {
     for(v=0;v<*ntime;v++){// spatial-temporal log-likelihood:
	   if(mlagt[t][v]<=*maxtime && mlags[i][j]<=*maxdist){
	        //Compute the temporal correlation function for a given pair:
      rho=CorFct(cormod,0,mlagt[t][v],parcor);
      //Compute the gradient of the temporal correlation function for a given pair:
      GradCorrFct(rho,cormod,eps,flagcor,gradcor,0,mlagt[t][v],parcor);
      //Compute the gradient of the conditional likelihood:
      Grad_Diff_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,
		      data[(t+*ntime*i)+n* *ntime* *ncoord],
		      data[(v+*ntime*i)+n* *ntime* *ncoord]);
      m=0;
      //Set the sensitivity matrix:
      for(d=0;d<*npar;d++){
	gradient[d]=gradient[d]+grad[d];
	score[d]=score[d]+grad[d];
	for(k=d;k<*npar;k++){
	  sensmat[m]=sensmat[m]+grad[d]*grad[k];
	  m++;}}}}}
	  }}}
    m=0;
    //Set the variability matrix:
    for(h=0;h<*npar;h++)
      for(k=h;k<*npar;k++){
	varimat[m]=varimat[m]+gradient[h]*gradient[k];
	m++;}}
  return;
}
// Empirical Godambe matrix in the pairwise Gaussian case:
void God_Pair_Gauss(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
		    int *npar, int *nparc, double *parcor, double *nuis, double *score,
		    double *sensmat, double *varimat)
{
  int d=0,i=0,j=0,h=0,k=0,m=0,n=0;
  double *gradcor, *grad, *gradient,rho=0.0;
  // Set the gradient vectors:
  gradcor=(double *) R_alloc(*nparc,sizeof(double));
  grad=(double *) R_alloc(*npar,sizeof(double));
  gradient=(double *) R_alloc(*npar,sizeof(double));

  for(n=0;n<*nrep;n++){
    // Initialize the gradient vector
    for(i=0;i<*npar;i++) gradient[i]=0;
    // Compute the gradient of a log likelihood object
    // Main loop, respect the pairs:
    for(i=0; i<(*ncoord-1);i++){
    for(j=(i+1); j<*ncoord;j++){
           if(mlags[i][j]<*maxdist){
      //Compute the correlation function
      rho=CorFct(cormod,mlags[i][j],0,parcor);
      //Compute the gradient for a given correlation model
      GradCorrFct(rho,cormod,eps,flagcor,gradcor,mlags[i][j],0,parcor);
      //Gradient of the log pairwise likelihood
      Grad_Pair_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,
		      data[(n+i* *nrep)], data[(n+j* *nrep)]);
      m=0;
      //Set the sensitivity matrix:
      for(d=0;d<*npar;d++){
	gradient[d]=gradient[d]+grad[d];
	score[d]=score[d]+grad[d];
	for(k=d;k<*npar;k++){
	  sensmat[m]=sensmat[m]+grad[d]*grad[k];
	  m++;}}}}}
    m=0;
    // Set the variability matrix:

    for(h=0;h<*npar;h++)
      for(k=h;k<*npar;k++){
	varimat[m]=varimat[m]+gradient[h]*gradient[k];
	m++;}}
  return;
}
// Empirical Godambe matrix in the pairwise spatio-temporal Gaussian case:
void God_Pair_Gauss_st(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
		       int *npar, int *nparc, double *parcor, double *nuis, double *score,
		       double *sensmat, double *varimat)
{
  int d=0,h=0,i=0,j=0,k=0,m=0,n=0,t=0,v=0;
  double *gradcor, *grad, *gradient,rho=0.0;

  gradcor=(double *) R_alloc(*nparc,sizeof(double));
  grad=(double *) R_alloc(*npar,sizeof(double));
  gradient=(double *) R_alloc(*npar,sizeof(double));

  for(n=0;n<*nrep;n++){
    // Initialize the gradient vector
    for(i=0;i<*npar;i++) gradient[i]=0;
    //Main loop, respect the pairs:
    for(i=0;i<*ncoord;i++){
   for(t=0;t<*ntime;t++){
     for(j=i;j<*ncoord;j++){
       if(i==j){// marginal temporal log-likelihood:
	 for(v=t+1;v<*ntime;v++){
	   if(mlagt[t][v]<=*maxtime){
      //Compute the temporal correlation function for a given pair:
      rho=CorFct(cormod,mlags[i][j],mlagt[t][v],parcor);
      //Compute the gradient of the temporal correlation function for a given pair:
      GradCorrFct(rho,cormod,eps,flagcor,gradcor,mlags[i][j],mlagt[t][v],parcor);
      //Compute the gradient of the conditional likelihood:
      Grad_Pair_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,
		      data[(t+*ntime*i)+n* *ntime* *ncoord],
		      data[(v+*ntime*j)+n* *ntime* *ncoord]);
      m=0;
      //Set the sensitivity matrix:
      for(d=0;d<*npar;d++){
	gradient[d]=gradient[d]+grad[d];
	score[d]=score[d]+grad[d];
	for(k=d;k<*npar;k++){
	  sensmat[m]=sensmat[m]+grad[d]*grad[k];
	  m++;}}}}}
else{
      for(v=0;v<*ntime;v++){// spatial-temporal log-likelihood:
	   if(mlagt[t][v]<=*maxtime && mlags[i][j]<=*maxdist){
      //Compute the temporal correlation function for a given pair:
      rho=CorFct(cormod,mlags[i][j],mlagt[t][v],parcor);
      //Compute the gradient of the temporal correlation function for a given pair:
      GradCorrFct(rho,cormod,eps,flagcor,gradcor,mlags[i][j],mlagt[t][v],parcor);
      //Compute the gradient of the conditional likelihood:
      Grad_Pair_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,
		      data[(t+*ntime*i)+n* *ntime* *ncoord],
		      data[(v+*ntime*j)+n* *ntime* *ncoord]);
      m=0;
      //Set the sensitivity matrix:
      for(d=0;d<*npar;d++){
	gradient[d]=gradient[d]+grad[d];
	score[d]=score[d]+grad[d];
	for(k=d;k<*npar;k++){
	  sensmat[m]=sensmat[m]+grad[d]*grad[k];
	  m++;}}}}}
}}}
    m=0;
    //Set the variability matrix:
    for(h=0;h<*npar;h++)
      for(k=h;k<*npar;k++){
	varimat[m]=varimat[m]+gradient[h]*gradient[k];
	m++;}}
  return;
}
// Empirical Godambe matrix in the pairwise Brown-Resnick case:
void God_BrowResn(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
		  int *npar, int *nparc, double *parcor, double *nuis, double *score,
		  double *sensmat, double *varimat)
{
  int d=0,i=0,j=0,h=0,k=0,m=0,n=0;
  double *gradcor, *grad, *gradient,vario=0.0;
  // Set the gradient vectors:
  gradcor=(double *) R_alloc(*nparc,sizeof(double));
  grad=(double *) R_alloc(*npar,sizeof(double));
  gradient=(double *) R_alloc(*npar,sizeof(double));

  for(n=0;n<*nrep;n++){
    // Initialize the gradient vector
    for(i=0;i<*npar;i++) gradient[i]=0;
    // Compute the gradient of a log likelihood object
    // Main loop, respect the pairs:
    for(i=0; i<(*ncoord-1);i++){
    for(j=(i+1); j<*ncoord;j++){
         if(mlags[i][j]<*maxdist){
      //Compute the variogram function
      vario=VarioFct(cormod,mlags[i][j],parcor,0);
      //Compute the gradient for a given variogram model
      GradVarioFct(vario,cormod,eps,flagcor,gradcor,mlags[i][j],parcor,0);
      //Gradient of the log pairwise likelihood
      Grad_Brow_Resn(vario,flagnuis,gradcor,grad,npar,nuis,
		     data[(n+i * *nrep)], data[(n+j * *nrep)]);
      m=0;
      //Set the sensitivity matrix:
      for(d=0;d<*npar;d++){
	gradient[d]=gradient[d]+grad[d];
	score[d]=score[d]+grad[d];
	for(k=d;k<*npar;k++){
	  sensmat[m]=sensmat[m]+grad[d]*grad[k];
	  m++;}}}}}
    m=0;
    // Set the variability matrix:
    for(h=0;h<*npar;h++)
      for(k=h;k<*npar;k++){
	varimat[m]=varimat[m]+gradient[h]*gradient[k];
	m++;}}
  return;
}
// Empirical Godambe matrix in the pairwise spatio-temporal Brown-Resnick case:
void God_BrowResn_st(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
		     int *npar, int *nparc, double *parcor, double *nuis, double *score,
		     double *sensmat, double *varimat)
{
  int d=0,h=0,i=0,j=0,k=0,m=0,n=0,t=0,v=0;
  double *gradcor,*grad,*gradient,vario=0.0;

  gradcor=(double *) R_alloc(*nparc,sizeof(double));
  grad=(double *) R_alloc(*npar,sizeof(double));
  gradient=(double *) R_alloc(*npar,sizeof(double));

  for(n=0;n<*nrep;n++){
    // Initialize the gradient vector
    for(i=0;i<*npar;i++) gradient[i]=0;
    //Main loop, respect the pairs:
    for(h=0;h<*npairs;h++){
      //Compute the temporal variogram function for a given pair:
      vario=VarioFct(cormod,mlags[i][j],parcor,mlagt[t][v]);
      //Compute the gradient of the temporal variogram function for a given pair:
      GradVarioFct(vario,cormod,eps,flagcor,gradcor,mlags[i][j],parcor,mlagt[t][v]);
      //Compute the gradient of the Brown-Resnick likelihood:
      Grad_Brow_Resn(vario,flagnuis,gradcor,grad,npar,nuis,
		     data[(t+*ntime*i)+n* *ntime* *ncoord],
		     data[(v+*ntime*j)+n* *ntime* *ncoord]);
      m=0;
      //Set the sensitivity matrix:
      for(d=0;d<*npar;d++){
	gradient[d]=gradient[d]+grad[d];
	score[d]=score[d]+grad[d];
	for(k=d;k<*npar;k++){
	  sensmat[m]=sensmat[m]+grad[d]*grad[k];
	  m++;}}}
    m=0;
    //Set the variability matrix:
        for(h=0;h<*npar;h++)
      for(k=h;k<*npar;k++){
	varimat[m]=varimat[m]+gradient[h]*gradient[k];
	m++;}}
  return;
}
// Empirical Godambe matrix in the pairwise Extremal Gaussian case:
void God_Ext_Gauss(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
		   int *npar, int *nparc, double *parcor, double *nuis, double *score,
		   double *sensmat, double *varimat)
{
  int d=0,i=0,j=0,h=0,k=0,m=0,n=0;
  double *gradcor, *grad, *gradient,rho=0.0;
  // Set the gradient vectors:
  gradcor=(double *) R_alloc(*nparc,sizeof(double));
  grad=(double *) R_alloc(*npar,sizeof(double));
  gradient=(double *) R_alloc(*npar,sizeof(double));

  for(n=0;n<*nrep;n++){
    // Initialize the gradient vector
    for(i=0;i<*npar;i++) gradient[i]=0;
    // Compute the gradient of a log likelihood object
    // Main loop, respect the pairs:
      for(i=0; i<(*ncoord-1);i++){
    for(j=(i+1); j<*ncoord;j++){
         if(mlags[i][j]<*maxdist){
      //Compute the correlation function
      rho=CorFct(cormod,mlags[i][j],0,parcor);
      //Compute the gradient for a given correlation model
      GradCorrFct(rho,cormod,eps,flagcor,gradcor,mlags[i][j],0,parcor);
      //Gradient of the log pairwise likelihood
      Grad_Ext_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,
		     data[(n+i* *nrep)],data[(n+j* *nrep)]);
      m=0;
      //Set the sensitivity matrix:
      for(d=0;d<*npar;d++){
	gradient[d]=gradient[d]+grad[d];
	score[d]=score[d]+grad[d];
	for(k=d;k<*npar;k++){
	  sensmat[m]=sensmat[m]+grad[d]*grad[k];
	  m++;}}}}}
    m=0;
    // Set the variability matrix:
    for(h=0;h<*npar;h++)
      for(k=h;k<*npar;k++){
	varimat[m]=varimat[m]+gradient[h]*gradient[k];
	m++;}}
  return;
}
// Empirical Godambe matrix in the pairwise Extremal t case:
void God_Ext_T(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
	       int *npar, int *nparc, double *parcor, double *nuis, double *score,
	       double *sensmat, double *varimat)
{
  int d=0,i=0,j=0,h=0,k=0,m=0,n=0;
  double *gradcor, *grad, *gradient,rho=0.0;
  // Set the gradient vectors:
  gradcor=(double *) R_alloc(*nparc,sizeof(double));
  grad=(double *) R_alloc(*npar,sizeof(double));
  gradient=(double *) R_alloc(*npar,sizeof(double));

  for(n=0;n<*nrep;n++){
    // Initialize the gradient vector
    for(i=0;i<*npar;i++) gradient[i]=0;
    // Compute the gradient of a log likelihood object
    // Main loop, respect the pairs:
     for(i=0; i<(*ncoord-1);i++){
    for(j=(i+1); j<*ncoord;j++){
             if(mlags[i][j]<*maxdist){
      //Compute the correlation function
      rho=CorFct(cormod,mlags[i][j],0,parcor);
      //Compute the gradient for a given correlation model
      GradCorrFct(rho,cormod,eps,flagcor,gradcor,mlags[i][j],0,parcor);
      //Gradient of the log pairwise likelihood
      Grad_Ext_T(rho,flagnuis,gradcor,grad,npar,nuis,
		 data[(n+i* *nrep)], data[(n+j* *nrep)]);
      m=0;
      //Set the sensitivity matrix:
      for(d=0;d<*npar;d++){
	gradient[d]=gradient[d]+grad[d];
	score[d]=score[d]+grad[d];
	for(k=d;k<*npar;k++){
	  sensmat[m]=sensmat[m]+grad[d]*grad[k];
	  m++;}}}}}
    m=0;
    // Set the variability matrix:
    for(h=0;h<*npar;h++)
      for(k=h;k<*npar;k++){
	varimat[m]=varimat[m]+gradient[h]*gradient[k];
	m++;}}
  return;
}
// Empirical estimation of the Senstive (H) and Variability (J) components of
// the Godambe matrix:
void GodambeMat_emp(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
		    int *like, int *model, int *npar, int *nparc, double *parcor, double *nuis,
		    double *score, double *sensmat, int *spt, double *varimat, int *type)
{
  //Select the model for computing the matrices H and J
  switch(*model){
  case 1:// Gaussian random field:
    switch(*like){
    case 1:// Conditional likelihood:
      if(*spt) God_Cond_Gauss_st(cormod,data,eps,flagcor,flagnuis,npar,
			         nparc,parcor,nuis,score,sensmat,varimat);
      else God_Cond_Gauss(cormod,data,eps,flagcor,flagnuis,npar,
			  nparc,parcor,nuis,score,sensmat,varimat);
      break;
    case 3:// Marginal likelihood:
      switch(*type){
      case 1:// Gradient of the log difference likelihood
	if(*spt) God_Diff_Gauss_st(cormod,data,eps,flagcor,flagnuis,npar,
				   nparc,parcor,nuis,score,sensmat,varimat);
	else God_Diff_Gauss(cormod,data,eps,flagcor,flagnuis,npar,
			    nparc,parcor,nuis,score,sensmat,varimat);
	break;
      case 2:// Gradient of the log pairwise likelihood
	if(*spt) God_Pair_Gauss_st(cormod,data,eps,flagcor,flagnuis,npar,
				   nparc,parcor,nuis,score,sensmat,varimat);
	else God_Pair_Gauss(cormod,data,eps,flagcor,flagnuis,npar,
			    nparc,parcor,nuis,score,sensmat,varimat);
	break;}
      break;}
    break;
  case 3:// Max-stable random field (Brown Resnick):
    if(*spt) God_BrowResn_st(cormod,data,eps,flagcor,flagnuis,npar,
			     nparc,parcor,nuis,score,sensmat,varimat);
    else God_BrowResn(cormod,data,eps,flagcor,flagnuis,npar,
		      nparc,parcor,nuis,score,sensmat,varimat);
    break;
  case 4:// Max-stable random field (Extremal Gaussian):
    God_Ext_Gauss(cormod,data,eps,flagcor,flagnuis,npar,
		   nparc,parcor,nuis,score,sensmat,varimat);
    break;
  case 5:// Max-stable random field (Extremal-t):
    God_Ext_T(cormod,data,eps,flagcor,flagnuis,npar,
	       nparc,parcor,nuis,score,sensmat,varimat);
    break;}
  return;
}

// The exact Senstive (H) and Variability (J) components of
// the Godambe matrix:
void GodambeMat(double *coordx, double *coordy, int *cormod, double *data, int *dist,
		      double *eps,int *flagcor, int *flagnuis, int *grid, int *like, int *model,
		      int *npar, int *nparc, double *parcor, double *nuis, double *score,
		      double *sensmat, int *spt, double *thr, int *type, double *varimat,
		      int *vartype, double *winc, double *winstp)
{
  //---------- COMPUTATION OF THE GODAMBE MATRIX ---------//
  int *np;np=(int *) R_alloc(1, sizeof(int));  //numbers of effective pairs
  switch(*vartype)
    {
    case 1://------------ START EMPIRICAL ESTIMATION ------------//
      GodambeMat_emp(cormod,data,eps,flagcor,flagnuis,like,
		     model,npar,nparc,parcor,nuis,score,
		     sensmat,spt,varimat,type);
      break;//------------ END EMPIRICAL ESTIMATION ------------//
    case 2://------------ START SUB-SAMPLE ESTIMATION ------------//
      Sensitivity(cormod,data,eps,flagcor,flagnuis,like,model,
		  npar,nparc,parcor,nuis,np,score,sensmat,spt,type);
      if(*spt) Vari_SubSamp_st(cormod,data,dist,eps,flagcor,flagnuis,
			       like,npar,nparc,nuis,np,parcor,type,varimat,winc,winstp);
      else Vari_SubSamp(coordx,coordy,cormod,data,dist,eps,flagcor,flagnuis,
			grid,like,model,npar,nparc,nuis,np,parcor,thr, type, varimat,winc,winstp);
      break;//------------ END SUB-SAMPLE ESTIMATION ------------//
    case 3://------------ START THEORETICAL COMPUTATION ------------//
      switch(*like)
	{
	case 1://----------- CONDITIONAL LIKELIHOOD ------------------------//
	  //------- GODAMBE FOR THE CONDITIONAL CASE -------//
	  break;
	case 3://----------- MARGINAL LIKELIHOOD ---------------------------//
	  switch(*type)
	    {
	    case 1://------- GODAMBE FOR THE DIFFERENCE CASE -------//
	      GodambeMat_Diff(coordx,coordy,cormod,eps,flagcor,flagnuis,
			      model,npar,nparc,parcor,nuis,sensmat,varimat);
	      break;
	    case 2://------- GODAMBE FOR THE PAIRWISE CASE -------//
	      break;
	    }
	  break;
	}
      break;//------------ END THEORETICAL COMPUTATION ------------//

    }
  return;
}
// Theoretical H and J matrices for the difference case:
void GodambeMat_Diff(double *coordx, double *coordy, int *cormod, double *eps, int *flagcor,
		     int *flagnuis, int *model, int *npar, int *nparc, double *parcor,
		     double *nuis, double *sensmat, double *varimat)
{

  double *gradcor_ij, *grad_ij, *gradcor_lk, *grad_lk;
  double rho_ij=0.0, rho_lk=0.0, crosscor=0.0, *vario;
  int s=0, l=0, k=0 , i=0, j=0, m=0, n=0;
  gradcor_ij=(double *) R_alloc(*nparc, sizeof(double));
  gradcor_lk=(double *) R_alloc(*nparc, sizeof(double));
  grad_ij=(double *) R_alloc(*npar, sizeof(double));
  grad_lk=(double *) R_alloc(*npar, sizeof(double));
  vario=(double *) R_alloc(6, sizeof(double));

  //--------------- GAUSSIAN MODEL --------------------------------//
  for(i=0;i<(*ncoord-1);i++){
  for(j=(i+1);j<*ncoord;j++){
	  if(mlags[i][j]<=*maxdist)
	    {// Compute the correlation function
	      rho_ij=CorFct(cormod,mlags[i][j],0,parcor);
	      // Compute the gradient for a given correlation value
	      GradCorrFct(rho_ij,cormod,eps,flagcor,gradcor_ij,mlags[i][j],0,parcor);
	      // Compute the gradient of the variogram for the difference Gaussian likelihood
	      Grad_Diff_Vario(rho_ij,flagnuis,gradcor_ij,grad_ij,npar,nuis);
	      // COMPUTATION OF THE SENSITIVITY MATRIX
	      Sens_Diff_Gauss_ij(grad_ij,npar,sensmat);

	      // START COMPUTATION OF THE VARIABILITY MATRIX
	      vario[0]=Variogram(cormod,mlags[i][j],0,nuis,parcor);

	      for(l=0;l<(*ncoord-1);l++)
		{// Insert the case of great circle distance
		  vario[2]=Variogram(cormod,hypot(coordx[i]-coordx[l],coordy[i]-coordy[l]),0,nuis,parcor);
		  vario[3]=Variogram(cormod,hypot(coordx[j]-coordx[l],coordy[j]-coordy[l]),0,nuis,parcor);

		  for(k=(l+1);k<*ncoord;k++)
		    {
		      if(mlags[l][k]<=*maxdist)
			{// Compute the correlation function
			  rho_lk=CorFct(cormod,mlags[l][k],0,parcor);
			  // Compute the gradient of a given correlation model

			  GradCorrFct(rho_lk,cormod,eps,flagcor,gradcor_lk,mlags[l][k],0,parcor);
			  // Compute the gradient of the variogram for the LK difference Gaussian likelihood
			  Grad_Diff_Vario(rho_lk,flagnuis,gradcor_lk,grad_lk,npar,nuis);

			  vario[1]=Variogram(cormod,mlags[l][k],0,nuis,parcor);
			  vario[4]=Variogram(cormod,hypot(coordx[i]-coordx[k],coordy[i]-coordy[k]),0,nuis,parcor);
			  vario[5]=Variogram(cormod,hypot(coordx[j]- coordx[k],coordy[j]-coordy[k]),0,nuis,parcor);

			  crosscor = (R_pow(vario[2]-vario[3]-vario[4]+vario[5],2))/(4*vario[0]*vario[1]);
			  s=0;

			  for(m=0;m<*npar;m++)
			    for(n=m;n<*npar;n++)
			      {
				varimat[s]=varimat[s]+0.5*grad_ij[m]*grad_lk[n]*crosscor;
				s++;
			      }
			  // END COMPUTATION OF THE VARIABILITY MATRIX
			}}}}}}

  return;
}
// Compute the Sensitivity matrix of a random field:
void Sensitivity(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis, int *like,
		 int *model, int *npar, int *nparc, double *parcor, double *nuis, int *np,double *score,
		 double *sensmat, int *spt, int *type)
{
  int b=0,d=0,i=0,j=0,n=0;
  double rho=0,*grad,*gradcor,vario=0;
  grad=(double *) R_alloc(*npar,sizeof(double));// gradient of the ijth composite log-likelihood
  gradcor=(double *) R_alloc(*nparc, sizeof(double));// gradient of the correlation
  switch(*model)// Compute the Sensitivity matrix
    {
    case 1:// Gaussian model
      switch(*like)
	{
	case 1:// Conditional likelihood:
	  if(*spt) Sens_Cond_Gauss_st(cormod,data,eps,flagcor,flagnuis,nuis,np,
				      npar,nparc,parcor,score,sensmat);
	  else Sens_Cond_Gauss(cormod,data,eps,flagcor,flagnuis,nuis,np,
			       npar,nparc,parcor,score,sensmat);
	  break;
	case 3: // Marginal likelihood:
	  switch(*type)
	    {
	    case 1: // Sensitivity for the difference likelihood case
	      if(*spt) Sens_Diff_Gauss_st(cormod,data,eps,flagcor,flagnuis,nuis,np,
				          npar,nparc,parcor,score,sensmat);
	      else Sens_Diff_Gauss(cormod,data,eps,flagcor,flagnuis,nuis,np,
				   npar,nparc,parcor,score,sensmat);
	      break;
	    case 2: // Sensitivity for the pairwise likelihood case
	      if(*spt) Sens_Pair_Gauss_st(cormod,data,eps,flagcor,flagnuis,nuis,np,
					  npar,nparc,parcor,score,sensmat);
	      else Sens_Pair_Gauss(cormod,data,eps,flagcor,flagnuis,nuis,np,
				   npar,nparc,parcor,score,sensmat);
	      break;
	    }
	  break;
	}
      break;
    case 2:// Binary gaussian random field:
      break;
    case 3:// Brown Resnick random field:
    for(i=0; i<(*ncoord-1);i++){
    for(j=(i+1); j<*ncoord;j++){
             if(mlags[i][j]<*maxdist){
	vario=VarioFct(cormod,mlags[i][j],parcor,0);
	GradVarioFct(vario,cormod,eps,flagcor,gradcor,mlags[i][j],parcor,0);
	for(n=0;n<*nrep;n++){
	  Grad_Brow_Resn(vario,flagnuis,gradcor,grad,npar,nuis,
			 data[(n+i* *nrep)],data[(n+j* *nrep)]);
	  for(d=0;d<*npar;d++) score[d]=score[d]+grad[d];}
	  b++;}}}
	  *np=b;
      break;
    case 4:// Extremal Gaussian random field:
    for(i=0; i<(*ncoord-1);i++){
    for(j=(i+1); j<*ncoord;j++){
                 if(mlags[i][j]<*maxdist){
	rho=CorFct(cormod,mlags[i][j],0,parcor);// Compute the correlation function
	//Compute the gradient for a given correlation model
	GradCorrFct(rho,cormod,eps,flagcor,gradcor,mlags[i][j],0,parcor);
	for(n=0;n<*nrep;n++){
	  Grad_Ext_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,
			 data[(n+i* *nrep)],data[(n+j* *nrep)]);
	  for(d=0;d<*npar;d++) score[d]=score[d]+grad[d];}
	  b++;}}}
	  *np=b;
    case 5:// Extremal-t random field:
    for(i=0; i<(*ncoord-1);i++){
    for(j=(i+1); j<*ncoord;j++){
                 if(mlags[i][j]<*maxdist){
	rho=CorFct(cormod,mlags[i][j],0,parcor);// Compute the correlation function
	//Compute the gradient for a given correlation model
	GradCorrFct(rho,cormod,eps,flagcor,gradcor,mlags[i][j],0,parcor);
	for(n=0;n<*nrep;n++){
	  Grad_Ext_T(rho,flagnuis,gradcor,grad,npar,nuis,
		     data[(n+i* *nrep)],data[(n+j* *nrep)]);
	  for(d=0;d<*npar;d++) score[d]=score[d]+grad[d];}
	  b++;}}}
	  *np=b;
      break;
    }

  return;
}
// Compute the Sensitivity matrix for the composite Gaussian difference likelihood:
void Sens_Diff_Gauss(int *cormod, double *data,double *eps, int *flagcor, int *flagnuis,
		     double *nuis, int *np,int *npar, int *nparc, double *parcor, double *score,
		     double *sensmat)
{
  // Initialization variables:
  int b=0,i=0,j=0,d=0,n=0;
  double rho=0.0, *grad, *gradcor, *gradient;

  gradcor=(double *) R_alloc(*nparc, sizeof(double));// Correlation gradient
  gradient=(double *) R_alloc(*npar, sizeof(double));// Overall gradient
  grad=(double *) R_alloc(*npar, sizeof(double));// ijth component of the gradient

     for(i=0; i<(*ncoord-1);i++){
    for(j=(i+1); j<*ncoord;j++){
             if(mlags[i][j]<*maxdist){
    //Compute the correlation function for the elements i,j
    rho=CorFct(cormod,mlags[i][j],0,parcor);
    //Compute the gradient for the given correlation
    GradCorrFct(rho,cormod,eps,flagcor,gradcor,mlags[i][j],0,parcor);
    //Compute the gradient of the variogram for the given pair
    Grad_Diff_Vario(rho,flagnuis,gradcor,gradient,npar,nuis);
    for(n=0;n<*nrep;n++){
      //Compute the gradient of the Gaussian difference composite
      Grad_Diff_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,
		      data[(n+i* *nrep)],data[(n+j* *nrep)]);
      //ADD TO THE SENSITIVITY MATRIX THE CONTRIBUTE OF THE GIVEN PAIR
      Sens_Diff_Gauss_ij(gradient,npar,sensmat);
      for(d=0;d<*npar;d++) score[d]=score[d]+grad[d];}
      b++;}}}
      *np=b;
  return;
}
// Compute the Sensitivity matrix for the space time composite Gaussian difference likelihood:
void Sens_Diff_Gauss_st(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis, double *nuis,int *np,
		        int *npar, int *nparc, double *parcor, double *score, double *sensmat)
{
  // Initialization variables:
  int b=0,d=0,i=0,j=0,n=0,t=0,v=0;
  double rho=0.0, *grad, *gradcor, *gradient;

  gradcor=(double *) R_alloc(*nparc, sizeof(double));// Correlation gradient
  gradient=(double *) R_alloc(*npar, sizeof(double));// Overall gradient
  grad=(double *) R_alloc(*npar, sizeof(double));// ijth component of the gradient
  for(i=0;i<*ncoord;i++){
    for(t=0;t<*ntime;t++){
      for(j=i;j<*ncoord;j++){
	if(i==j){// marginal temporal log-likelihood:
	  for(v=t+1;v<*ntime;v++){
	    if(mlagt[t][v]<=*maxtime){
	      // Compute the temporal correlation function for given pair:
	      rho=CorFct(cormod,0,mlagt[t][v],parcor);
	      // Compute the gradient for the given temporal correlation
	      GradCorrFct(rho,cormod,eps,flagcor,gradcor,0,mlagt[t][v],parcor);
	      // Compute the gradient of the variogram for the given pair (temporal):
	      Grad_Diff_Vario(rho,flagnuis,gradcor,gradient,npar,nuis);
	      for(n=0;n<*nrep;n++){
		//Compute the gradient of the Gaussian difference composite
		Grad_Diff_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,
				data[(t+*ntime*i)+n* *ntime* *ncoord],
				data[(v+*ntime*j)+n* *ntime* *ncoord]);
		//ADD TO THE SENSITIVITY MATRIX THE CONTRIBUTE OF THE GIVEN PAIR (temporal):
		Sens_Diff_Gauss_ij(gradient,npar,sensmat);
		for(d=0;d<*npar;d++) score[d]=score[d]+grad[d];}
	      b++;}}}
	else{// spatial-temporal log-likelihood:
	  for(v=0;v<*ntime;v++){
	    if(mlagt[t][v]<=*maxtime && mlags[i][j]<=*maxdist){
	      // Compute the spatial-temporal correlation function for a given pair
	      rho=CorFct(cormod,mlags[i][j],mlagt[t][v],parcor);
	      // Compute the gradient for the given correlation
	      GradCorrFct(rho,cormod,eps,flagcor,gradcor,mlags[i][j],mlagt[t][v],parcor);
	      // Compute the gradient of the variogram for the given pair
	      Grad_Diff_Vario(rho,flagnuis,gradcor,gradient,npar,nuis);
	      for(n=0;n<*nrep;n++){
		//Compute the gradient of the Gaussian difference composite
		Grad_Diff_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,
				data[(t+*ntime*i)+n* *ntime* *ncoord],
				data[(v+*ntime*j)+n* *ntime* *ncoord]);
		//ADD TO THE SENSITIVITY MATRIX THE CONTRIBUTE OF THE GIVEN PAIR
		Sens_Diff_Gauss_ij(gradient,npar,sensmat);
		for(d=0;d<*npar;d++) score[d]=score[d]+grad[d];}
	      b++;}}}}}}
	      *np=b;
  return;
}
// Compute the Sensitivity matrix for a single Gaussian difference likelihood component:
void Sens_Diff_Gauss_ij(double *gradient, int *npar, double *sensmat)
{
  // Initialization variables:
  int h=0,i=0,j=0;

  for(i=0;i<*npar;i++)
    for(j=i;j<*npar;j++){
      sensmat[h]=sensmat[h]+0.5*gradient[i]*gradient[j];
      h++;}
  return;
}

// Compute the Sensitivity matrix for the pairwise composite Gaussian likelihood:
void Sens_Pair_Gauss(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
		     double *nuis, int *np, int *npar, int *nparc, double *parcor,double *score, double *sensmat)
{
  // Initialization variables:
  int b=0,i=0,j=0,l=0,d=0,n=0,nsens=0;
  double rho=0.0, *grad, *gradcor, *sens;

  nsens=*npar * (*npar+1)/2;
  gradcor=(double *) R_alloc(*nparc, sizeof(double));// Correlation gradient
  grad=(double *) R_alloc(*npar, sizeof(double));// ijth component of the gradient
  sens=(double *) R_alloc(nsens, sizeof(double));// One sensitive contribute

    for(i=0; i<(*ncoord-1);i++){
    for(j=(i+1); j<*ncoord;j++){
    //Compute the correlation function for the elements i,j
      if(mlags[i][j]<*maxdist){
    rho=CorFct(cormod,mlags[i][j],0,parcor);
    //Compute the gradient for the given correlation
    GradCorrFct(rho,cormod,eps,flagcor,gradcor,mlags[i][j],0,parcor);
    for(n=0;n<*nrep;n++){
      //Compute the gradient of the log pairwise likelihood
      Grad_Pair_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,
		      data[(n+i* *nrep)],data[(n+j* *nrep)]);
      //ADD TO THE SENSITIVITY MATRIX THE CONTRIBUTE OF THE GIVEN PAIR
      Sens_Pair_Gauss_ij(rho,flagnuis,gradcor,npar,nparc,nuis,sens);
      for(l=0;l<nsens;l++) sensmat[l]=sensmat[l]-sens[l];
      for(d=0;d<*npar;d++) score[d]=score[d]+grad[d];}
      b++;}}}
      *np=b;
  return;
}
// Compute the Sensitivity matrix for the space time pairwise composite Gaussian likelihood:
void Sens_Pair_Gauss_st(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis, double *nuis, int *np,
			int *npar, int *nparc, double *parcor, double *score, double *sensmat)
{
  // Initialization variables:
  int  b=0,d=0,i=0,l=0,n=0,nsens=0,j=0,t=0,v=0;
  double rho=0.0, *grad, *gradcor, *sens;

  nsens=*npar * (*npar+1)/2;
  gradcor=(double *) R_alloc(*nparc, sizeof(double));// Correlation gradient
  grad=(double *) R_alloc(*npar, sizeof(double));// ijth component of the gradient
  sens=(double *) R_alloc(nsens, sizeof(double));// One sensitive contribute
  for(i=0;i<*ncoord;i++){
    for(t=0;t<*ntime;t++){
      for(j=i;j<*ncoord;j++){
	if(i==j){// marginal temporal log-likelihood:
	  for(v=t+1;v<*ntime;v++){
	    if(mlagt[t][v]<=*maxtime){
	      // Compute the temporal correlation function for for a given pair:
	      rho=CorFct(cormod,0,mlagt[t][v],parcor);
	      // Compute the gradient for the temporal correlation function:
	      GradCorrFct(rho,cormod,eps,flagcor,gradcor,0,mlagt[t][v],parcor);
	      for(n=0;n<*nrep;n++){
		//Compute the gradient of the log pairwise likelihood
		Grad_Pair_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,
				data[(t+*ntime*i)+n* *ntime* *ncoord],
				data[(v+*ntime*j)+n* *ntime* *ncoord]);
		//ADD TO THE SENSITIVITY MATRIX THE CONTRIBUTE OF THE GIVEN PAIR
		Sens_Pair_Gauss_ij(rho,flagnuis,gradcor,npar,nparc,nuis,sens);
		for(l=0;l<nsens;l++) sensmat[l]=sensmat[l]-sens[l];
		for(d=0;d<*npar;d++) score[d]=score[d]+grad[d];}
	      b++;}}}
	else {// marginal spatial-temporal log-likelihood:
	  for(v=0;v<*ntime;v++){
	    if(mlagt[t][v]<=*maxtime && mlags[i][j]<=*maxdist){
	      // Compute the spatial-temporal correlation function for a given pair:
	      rho=CorFct(cormod,mlags[i][j],mlagt[t][v],parcor);
	      // Compute the gradient of the patial-temporal correlation function:
	      GradCorrFct(rho,cormod,eps,flagcor,gradcor,mlags[i][j],mlagt[t][v],parcor);
	      for(n=0;n<*nrep;n++){
		//Compute the gradient of the log pairwise likelihood
		Grad_Pair_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,
				data[(t+*ntime*i)+n* *ntime* *ncoord],
				data[(v+*ntime*j)+n* *ntime* *ncoord]);
		//ADD TO THE SENSITIVITY MATRIX THE CONTRIBUTE OF THE GIVEN PAIR
		Sens_Pair_Gauss_ij(rho,flagnuis,gradcor,npar,nparc,nuis,sens);
		for(l=0;l<nsens;l++) sensmat[l]=sensmat[l]-sens[l];
		for(d=0;d<*npar;d++) score[d]=score[d]+grad[d];}
		b++;
		}}}
          }}}
          *np=b;
  return;
}
// Compute the Sensitivity matrix for the Gaussian difference likelihood:
void Sens_Pair_Gauss_ij(double rho, int *flag, double *gradcor, int *npar,
			int *nparc, double *par, double *sensmat)
{
  // Initialization variables:
  double nugget=par[1], sill=par[2];
  double C=0.0, ER=0.0, EL=0.0, ES=0.0, EM=0.0, EQ=0.0;
  double a=0.0, b=0.0, c=0.0, d=0.0, e=0.0;
  double f=0.0, q=0.0, t=0.0;
  int h=0, i=0, l=0, j=0, p=0;
  // Defines useful quantities:
  a=nugget+sill;
  b=sill*rho;
  c=pow(a,2)-pow(b,2);
  d=a/c;
  e=b/c;
  f=a+b;
  q=d-rho*e;
  t=1-pow(rho,2);
  // Expected values:
  ER=2*a;
  EL=b;
  ES=0.5*ER/c;
  EM=EL/c;
  EQ=d*ER-2*e*EL-1;
  //---- START COMPUTATION THE SENSITIVITY MATRIX-----//
  // Derivatives of the PLLik respect with the mean
  if(flag[0]==1)
    {
      sensmat[p]=-2/f; // Second derivative
       i++; p++; l++;
      // Mixed derivatives:
      if(flag[1]==1){ sensmat[p]=0; p++; l++; }// nugget
      if(flag[2]==1){ sensmat[p]=0; p++; l++; }// sill
      for(j=l;j<*npar;j++){ sensmat[p]=0; p++;}// correlation parameters
    }
  // Derivative of the PLLik respect with the nugget
  if(flag[1]==1)
    {
      l = i; i++;
      sensmat[p]=EQ*(1/c-4*pow(d,2))+(2*d*ER-1)/c;
      p++; l++;
      if(flag[2]==1) // sill
	{
	  sensmat[p]=2*(ES*(3*d-rho*e)-EM*(d*rho+e)-
			d*(d-rho*e)*(2*EQ+1))-1/c;
	  p++; l++;
	}
      h=0;
      for(j=l;j<*npar;j++) // correlation parameters
	{
	  sensmat[p]=2*sill*(d*e*(2*EQ+1)-e*ES-d*EM)*gradcor[h];
	  h++; p++;
	}
    }
  // Derivative of the PLLik respect with the sill
  if(flag[2]==1)
    {
      l = i; i++;
      sensmat[p]=-2*pow(q,2)*(2*EQ+1)+
	2*(ES*(2*q+d*t)-EM*(2*rho*q+e*t))-t/c;
      p++; l++;
      h=0;
      for(j=l;j<*npar;j++) // correlation parameters
	{
	  sensmat[p]=(2*e*(sill*(q*(2*EQ+1)-ES)-EQ+e*b)+
		      EM*(1-2*sill*q))*gradcor[h];
	  h++; p++;
	}
    }
  // Derivative of the PLLik respect with the correlation parameters*/
  C=-pow(sill/c,2)*(2*pow(b,2)+c);

  for(j=0;j<*nparc;j++)
    for(l=j;l<*nparc;l++)
      {
	sensmat[p]=C*gradcor[j]*gradcor[l];
	p++;
      }
  //---- END COMPUTATION OF THE SENSITIVITY MATRIX-----//
  return;
}

// Compute the Sensitivity matrix for the conditional composite Gaussian likelihood:
void Sens_Cond_Gauss(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis, double *nuis, int *np,
		     int *npar, int *nparc, double *parcor, double *score, double *sensmat)
{
  // Initialization variables:
  int b=0,d=0,i=0,j=0,l=0,n=0,nsens=0;
  double rho=0.0, *grad, *gradcor, *sens;

  nsens=*npar * (*npar+1)/2;
  gradcor=(double *) R_alloc(*nparc, sizeof(double));// Correlation gradient
  grad=(double *) R_alloc(*npar, sizeof(double));// ijth component of the gradient
  sens=(double *) R_alloc(nsens, sizeof(double));// One sensitive contribute

      for(i=0; i<(*ncoord-1);i++){
    for(j=(i+1); j<*ncoord;j++){
               if(mlags[i][j]<*maxdist){
    //Compute the correlation function for the elements i,j
    rho=CorFct(cormod,mlags[i][j],0,parcor);
    //Compute the gradient for the given correlation
    GradCorrFct(rho,cormod,eps,flagcor,gradcor,mlags[i][j],0,parcor);
    for(n=0;n<*nrep;n++){
      //Compute the gradient of the conditional likelihood:
      Grad_Cond_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,
		      data[(n+i* *nrep)],data[(n+j* *nrep)]);
      //ADD TO THE SENSITIVITY MATRIX THE CONTRIBUTE OF THE GIVEN PAIR
      Sens_Cond_Gauss_ij(rho,flagnuis,gradcor,npar,nparc,nuis,sens);
      for(l=0;l<nsens;l++) sensmat[l]=sensmat[l]-sens[l];
      for(d=0;d<*npar;d++) score[d]=score[d]+grad[d];}
      b++;}}}
      *np=b;
  return;
}
// Compute the Sensitivity matrix for the space time conditional composite Gaussian likelihood:
void Sens_Cond_Gauss_st(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis, double *nuis, int *np,
			int *npar, int *nparc, double *parcor, double *score, double *sensmat)
{
  // Initialization variables:
  int b=0,d=0,i=0,l=0,n=0,nsens=0,j=0,t=0,v=0;
  double rho=0.0, *grad, *gradcor, *sens;

  nsens=*npar * (*npar+1)/2;
  gradcor=(double *) R_alloc(*nparc, sizeof(double));// Correlation gradient
  grad=(double *) R_alloc(*npar, sizeof(double));// ijth component of the gradient
  sens=(double *) R_alloc(nsens, sizeof(double));// One sensitive contribute

    for(i=0;i<*ncoord;i++){
      for(t=0;t<*ntime;t++){
	for(j=i;j<*ncoord;j++){
	  if(i==j){// marginal temporal log-likelihood:
	    for(v=t+1;v<*ntime;v++){
	      if(mlagt[t][v]<=*maxtime){
		//Compute the temporal correlation function for a given pair:
		rho=CorFct(cormod,0,mlagt[t][v],parcor);
		//Compute the gradient of the temporal correlation function for a given pair:
		GradCorrFct(rho,cormod,eps,flagcor,gradcor,0,mlagt[t][v],parcor);
		for(n=0;n<*nrep;n++){
		  //Compute the gradient of the conditional likelihood:
		  Grad_Cond_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,
				  data[(t+*ntime*i)+n* *ntime* *ncoord],
				  data[(v+*ntime*j)+n* *ntime* *ncoord]);
		  //ADD TO THE SENSITIVITY MATRIX THE CONTRIBUTE OF THE GIVEN PAIR
		  Sens_Cond_Gauss_ij(rho,flagnuis,gradcor,npar,nparc,nuis,sens);
		  for(l=0;l<nsens;l++) sensmat[l]=sensmat[l]-sens[l];
		  for(d=0;d<*npar;d++) score[d]=score[d]+grad[d];}
		  b++;}}}
	else{// marginal spatial-temporal log-likelihood:
	  for(v=0;v<*ntime;v++){
	    if(mlagt[t][v]<=*maxtime && mlags[i][j]<=*maxdist){
	      // Compute the spatial-temporal correlation function for a given pair:
	      rho=CorFct(cormod,mlags[i][j],mlagt[t][v],parcor);
	      // Compute the gradient of the spatial-temporal correlation  for a given pair:
	      GradCorrFct(rho,cormod,eps,flagcor,gradcor,mlags[i][j],mlagt[t][v],parcor);
	      for(n=0;n<*nrep;n++){
		//Compute the gradient of the conditional likelihood:
		Grad_Cond_Gauss(rho,flagnuis,gradcor,grad,npar,nuis,
				data[(t+*ntime*i)+n* *ntime* *ncoord],
				data[(v+*ntime*j)+n* *ntime* *ncoord]);
		//ADD TO THE SENSITIVITY MATRIX THE CONTRIBUTE OF THE GIVEN PAIR
		Sens_Cond_Gauss_ij(rho,flagnuis,gradcor,npar,nparc,nuis,sens);
		for(l=0;l<nsens;l++) sensmat[l]=sensmat[l]-sens[l];
		for(d=0;d<*npar;d++) score[d]=score[d]+grad[d];}
	      b++;}}}}}}
	      *np=b;
    return;
}
// Compute the Sensitivity matrix for a single Gaussian conditional likelihood:
void Sens_Cond_Gauss_ij(double rho, int *flag, double *gradcor, int *npar,
			int *nparc, double *par, double *sensmat)
{
 // Initialization variables:
 double nugget=par[1], sill=par[2];
 double C=0.0, ER=0.0, EL=0.0, ES=0.0, EM=0.0, EQ=0.0;
 double a=0.0, b=0.0, c=0.0, d=0.0;
 double e=0.0, f=0.0, q=0.0, t=0.0;
 int h=0, i=0, l=0, j=0, p=0;
 // Defines useful quantities:
 a=nugget+sill;
 b=sill*rho;
 c=pow(a,2)-pow(b,2);
 d=a/c;
 e=b/c;
 f=a+b;
 q=d-rho*e;
 t=1-pow(rho,2);
 // Expected values:
 ER=2*a;
 EL=b;
 ES=0.5*ER/c;
 EM=EL/c;
 EQ=d*ER-2*e*EL-1;
 //---- START COMPUTATION THE SENSITIVITY MATRIX-----//
 // Derivatives of the PLLik respect with the mean
 if(flag[0]==1)
   {
     sensmat[p]=-4/f+2/a; // Second derivative
      i++; p++; l++;
     // Mixed derivatives:
     if(flag[1]==1){ sensmat[p]=0; p++; l++; }// nugget
     if(flag[2]==1){ sensmat[p]=0; p++; l++; }// sill
     for(j=l;j<*npar;j++){ sensmat[p]=0;p++; }// correlation parameters
   }
 // Derivative of the PLLik respect with the nugget
 if(flag[1]==1)
   {
     l = i; i++;
     sensmat[p]=2*(EQ*(1/c-4*pow(d,2))+(2*d*ER-1)/c)-1/pow(a,2);
     p++; l++;
     if(flag[2]==1)// sill
	{
	  sensmat[p]=2*(2*(ES*(3*d-rho*e)-EM*(d*rho+e)-
			   d*(d-rho*e)*(2*EQ+1))-1/c)-1/pow(a,2);
	  p++; l++;
	}
     h=0;
     for(j=l;j<*npar;j++) // correlation parameters
	{
	  sensmat[p]=4*sill*(d*e*(2*EQ+1)-e*ES-d*EM)*gradcor[h];
	  h++; p++;
	}
   }
 // Derivative of the PLLik respect with the sill
 if(flag[2]==1)
   {
     l = i; i++;
     sensmat[p]=2*(-2*pow(q,2)*(2*EQ+1)+
		   2*(ES*(2*q+d*t)-EM*(2*rho*q+e*t))-t/c)+1/pow(a,2);
     p++; l++;
     h=0;
     for(j=l;j<*npar;j++) // correlation parameters
	{
	  sensmat[p]=2*((2*e*(sill*(q*(2*EQ+1)-ES)-EQ+e*b)+
			 EM*(1-2*sill*q))*gradcor[h]);
	  h++; p++;
	}
   }
 // Derivative of the PLLik respect with the correlation parameters
 C=-pow(sill/c,2)*(2*pow(b,2)+c);
 for(j=0;j<*nparc;j++)
   for(l=j;l<*nparc;l++)
     {
	sensmat[p]=2*C*gradcor[j]*gradcor[l];
	p++;
     }

 //---- END COMPUTATION OF THE SENSITIVITY MATRIX-----//
 return;
}

void Vari_SubSamp(double *coordx, double *coordy, int *cormod, double *data,
		  int *dist,double *eps, int *flagcor, int *flagnuis, int *grid,
            int *like, int *model, int *npar, int *nparc, double *nuis, int *np,
		  double *parcor, double *thr, int *type, double *varimat,
		  double *winc, double *winstp)
{
  double rho=0.0, lag=0.0, *gradcor,*gradient, psm=0, ps, *rangex, *rangey;
  double *ecoordx,*ecoordy,*scoordx,*scoordy,*sdata,*sumgrad,*subvari,*xgrid,*ygrid;
  double delta=0.0, deltax=0.0, deltay=0.0, dimwinx=0.0, dimwiny=0.0, vario;
  double winstx=winstp[0], winsty=winstp[0];
  int *npts, numintx=0, numinty=0;
  int h=0,i=0,l=0,m=0,n=0,nsub=0,nvari=0,nwpair=0,p=0,q=0,j=0;

  nvari=*npar * (*npar+1)/2;

  gradcor=(double *) R_alloc(*nparc, sizeof(double));
  gradient=(double *) R_alloc(*npar, sizeof(double));
  sumgrad=(double *) R_alloc(*npar, sizeof(double));
  subvari=(double *) R_alloc(nvari, sizeof(double));

  npts=(int *) R_alloc(1, sizeof(int));

  rangex=(double *) R_alloc(2, sizeof(double));
  rangey=(double *) R_alloc(2, sizeof(double));

  scoordx=(double *) R_alloc(*ncoord, sizeof(double));
  scoordy=(double *) R_alloc(*ncoord, sizeof(double));
  sdata=(double *) R_alloc(*ncoord, sizeof(double));

  if(*grid){
    ecoordx=(double *) R_alloc(*ncoord, sizeof(double));
    ecoordy=(double *) R_alloc(*ncoord, sizeof(double));
    for(i=0;i<*ncoordx;i++)
      for(j=0;j<*ncoordy;j++){
	ecoordx[h]=coordx[i];
	ecoordy[h]=coordy[j];
	h++;}
    coordx=ecoordx;
    coordy=ecoordy;}

  Range(coordx,rangex,ncoord);// range of the x-coordinate
  Range(coordy,rangey,ncoord);// range of the y-coordinate
  // set the sub-sampling window based on prototype unit window (R_0)
  // and scaling factor (lambda_n)
  deltax=rangex[1]-rangex[0];// R_n = lambda_n * R_0
  deltay=rangey[1]-rangey[0];
  if(!winc[0]){
    delta=fmin(deltax,deltay);
    winc[0]=(delta/sqrt(delta))/2;}
  if(!winstp[0]) winstp[0]=0.5;
  dimwinx=winc[0] * sqrt(deltax);// sub-window x length depends on a constant: deafault??
  dimwiny=winc[0] * sqrt(deltay);// sub-window y length depends on a constant: deafault??
  winstx=*winstp * dimwinx;     // x step is a  proportion of sub-window x length (deafult is 0.5)
  winsty=*winstp * dimwiny;     // y step is a  proportion of sub-window y length (deafult is 0.5)
  numintx=floor((deltax-dimwinx)/winstx+1);   //number of overlapping sub-windows
  numinty=floor((deltay-dimwiny)/winsty+1);   //number of overlapping sub-windows

  xgrid=(double *) R_alloc(numintx, sizeof(double));
  ygrid=(double *) R_alloc(numinty, sizeof(double));

  SeqStep(rangex,numintx,winstx,xgrid);
  SeqStep(rangey,numinty,winsty,ygrid);

  for(n=0;n<*nrep;n++){
    for(h=0;h<nvari;h++) subvari[h]=0;//initialize the  variance on the subwindows
    nsub=0;
    for(i=0;i<numintx;i++){
      for(j=0;j<numinty;j++){
	*npts=0;
	for(h=0;h<*npar;h++){ sumgrad[h]=0;gradient[h]=0;} //initialize the gradient of the sub-window
	for(h=0;h<*nparc;h++)gradcor[h]=0;
	SetSampling(coordx,coordy,data,n,npts,scoordx,scoordy,
		    sdata,xgrid[i]+dimwinx,xgrid[i],ygrid[j]+dimwiny,ygrid[j]);//  create data and coordinates of the sub-windows
	if(*npts>4){
	  nsub++;
	  nwpair=0;//initialize the number of pairs in the window
	  switch(*model){
	  case 1: // Gaussian random field:
	    for(l=0;l<(*npts-1);l++){
	      for(m=(l+1);m<*npts;m++){
           if(*dist==0) lag=hypot(scoordx[l]-scoordx[m],scoordy[l]-scoordy[m]);
		 if(*dist==1) lag=Dist_chordal(scoordx[l],scoordy[l],scoordx[m],scoordy[m]);
           if(*dist==2) lag=Dist_geodesic(scoordx[l],scoordy[l],scoordx[m],scoordy[m]);
                if(lag<*maxdist){
		  nwpair++;
		  rho=CorFct(cormod,lag,0,parcor);
		  GradCorrFct(rho,cormod,eps,flagcor,gradcor,lag,0,parcor);
               switch(*like){//select the type of composite likelihood
	       case 1:// Conditional likelihood:
		 switch(*type){
		 case 2:
		   Grad_Cond_Gauss(rho,flagnuis,gradcor,gradient,npar,nuis,sdata[l],sdata[m]);
		   break;}
		 break;
               case 3: // Marginal likelihood:
		 switch(*type){
		 case 1: // Gradient of the log difference likelihood
		   Grad_Diff_Gauss(rho,flagnuis,gradcor,gradient,npar,nuis,sdata[l],sdata[m]);
		   break;
		 case 2:// Gradient of the log pairwise likelihood
		   Grad_Pair_Gauss(rho,flagnuis,gradcor,gradient,npar,nuis,sdata[l],sdata[m]);
		   break;}
		 break;}
	       for(h=0;h<*npar;h++) {sumgrad[h]=sumgrad[h]+gradient[h]; }// sum the gradient in the subwindow
		}}}
	    break;
	  case 2: // Binary-Gaussian random field:
	    psm=pnorm((nuis[0]-*thr)/sqrt(nuis[2]+nuis[1]),0,1,1,0);
	    for(l=0;l<(*npts-1);l++){
	      for(m=(l+1);m<*npts;m++){
           if(*dist==0) lag=hypot(scoordx[l]-scoordx[m],scoordy[l]-scoordy[m]);
		 if(*dist==1) lag=Dist_chordal(scoordx[l],scoordy[l],scoordx[m],scoordy[m]);
           if(*dist==2) lag=Dist_geodesic(scoordx[l],scoordy[l],scoordx[m],scoordy[m]);
		if(lag<*maxdist){
		  nwpair++;
		  rho=CorFct(cormod,lag,0,parcor);
		  GradCorrFct(rho,cormod,eps,flagcor,gradcor,lag,0,parcor);
		  ps=pbnorm(cormod,lag,0,nuis,parcor,*thr);
               switch(*like){//select the type of composite likelihood
               case 1:// Conditional likelihood:
		 switch(*type){
		 case 2:// Gradient of the log binnary Gauss conditional lik
		   Grad_Cond_Bin(rho,ps,psm,flagnuis,gradcor,gradient,npar,nuis,thr,sdata[l],sdata[m]);
		   break;}
		 break;
               case 3:// Marginal likelihood:
		 switch(*type){
		 case 1:// Gradient of the log difference likelihood
		   Grad_Diff_Bin(rho,ps,psm,flagnuis,gradcor,gradient,npar,nuis,thr,sdata[l],sdata[m]);
		   //Insert the gradient of binnary Gauss Diff
		   break;
		 case 2:// Gradient of the log pairwise likelihood
		   Grad_Pair_Bin(rho,ps,psm,flagnuis,gradcor,gradient,npar,nuis,thr,sdata[l],sdata[m]);
		   break;}
		 break;}
               for(h=0;h<*npar;h++) sumgrad[h]=sumgrad[h]+gradient[h];// sum the gradient in the subwindow
		}}}
	    break;
	  case 3:// Brown-Resnick random field:
	    for(l=0;l<(*npts-1);l++){
	      for(m=(l+1);m<*npts;m++){
           if(*dist==0) lag=hypot(scoordx[l]-scoordx[m],scoordy[l]-scoordy[m]);
		 if(*dist==1) lag=Dist_chordal(scoordx[l],scoordy[l],scoordx[m],scoordy[m]);
           if(*dist==2) lag=Dist_geodesic(scoordx[l],scoordy[l],scoordx[m],scoordy[m]);
		if(lag<*maxdist){
		  nwpair++;
		  vario=VarioFct(cormod,lag,parcor,0);
		  GradVarioFct(vario,cormod,eps,flagcor,gradcor,lag,parcor,0);
		  switch(*like){//select the type of composite likelihood
		  case 3:// Marginal likelihood:
		    switch(*type){
		    case 2:// Gradient of the log pairwise likelihood
		      Grad_Brow_Resn(vario,flagnuis,gradcor,gradient,npar,nuis,sdata[l],sdata[m]);
		      break;}
		    break;}
		  for(h=0;h<*npar;h++) sumgrad[h]=sumgrad[h]+gradient[h];   // sum the gradient in the subwindow
		}}}
	    break;
	  case 4:// Extremal Gaussian random field:
	    for(l=0;l<(*npts-1);l++){
	      for(m=(l+1);m<*npts;m++){
           if(*dist==0) lag=hypot(scoordx[l]-scoordx[m],scoordy[l]-scoordy[m]);
		 if(*dist==1) lag=Dist_chordal(scoordx[l],scoordy[l],scoordx[m],scoordy[m]);
           if(*dist==2) lag=Dist_geodesic(scoordx[l],scoordy[l],scoordx[m],scoordy[m]);
		if(lag<*maxdist){
		  nwpair++;
		  rho=CorFct(cormod,lag,0,parcor);
		  GradCorrFct(rho,cormod,eps,flagcor,gradcor,lag,0,parcor);
		  Grad_Ext_Gauss(rho,flagnuis,gradcor,gradient,npar,nuis, sdata[l],sdata[m]);
		  for(h=0;h<*npar;h++) sumgrad[h]=sumgrad[h]+gradient[h];  // sum the gradient in the subwindow
		}}}
	    break;
	  case 5:// Extremal t random field:
	    for(l=0;l<(*npts-1);l++){
	      for(m=(l+1);m<*npts;m++){
		 if(*dist==0) lag=hypot(scoordx[l]-scoordx[m],scoordy[l]-scoordy[m]);
		 if(*dist==1) lag=Dist_chordal(scoordx[l],scoordy[l],scoordx[m],scoordy[m]);
           if(*dist==2) lag=Dist_geodesic(scoordx[l],scoordy[l],scoordx[m],scoordy[m]);
		if(lag<*maxdist){
		  nwpair++;
		  rho=CorFct(cormod,lag,0,parcor);
		  GradCorrFct(rho,cormod,eps,flagcor,gradcor,lag,0,parcor);
		  switch(*like){//select the type of composite likelihood
		  case 3:// Marginal likelihood:
		    switch(*type){
		    case 2:// Gradient of the log pairwise likelihood
		      Grad_Ext_T(rho,flagnuis,gradcor,gradient,npar,nuis,sdata[l],sdata[m]);
		      break;}
		    break;}
		  for(h=0;h<*npar;h++) sumgrad[h]=sumgrad[h]+gradient[h];// sum the gradient in the subwindow
		}}}
	    break;}
	  h=0;
	  for(p=0;p<*npar;p++)//update the sub-variance in the subwindow
	    for(q=p;q<*npar;q++){subvari[h]=subvari[h]+sumgrad[p]*sumgrad[q]/nwpair;h++;}
	}}}}
  for(h=0;h<nvari;h++) varimat[h]= np[0] * subvari[h]/(nsub* *nrep);//update variability matrix
  return;
}
// Computes the variability matrix based on the sub-sampling method (spatial-temporal):
void Vari_SubSamp_st(int *cormod, double *data, int *dist, double *eps, int *flagcor,
		           int *flagnuis,int *like,int *npar, int *nparc, double *nuis,int *np, double *parcor,
		           int *type,double *varimat, double *winc, double *winstp)
{
  double beta, rho=0.0, lagt=0.0, *gradcor, *gradient;
  double *sdata, *sumgrad, *subvari, *sublagt;
  double step=*minimtime;
  int nsub=0, nstime=0, nvari=0, nwpair=0;
  int h=0, i=0, k=0, n=0, j=0, p=0, q=0, t=0, v=0;
  nvari=*npar * (*npar+1)/2;
  gradcor=(double *) R_alloc(*nparc, sizeof(double));
  gradient=(double *) R_alloc(*npar, sizeof(double));
  sumgrad=(double *) R_alloc(*npar, sizeof(double));
  subvari=(double *) R_alloc(nvari, sizeof(double));

  //default sub window temporal length
  if(!(*winc)){
    beta=CorFct(cormod,0,1,parcor);
    *winc=R_pow(2*beta/(1-R_pow(beta,2)),(double) 2/3)*R_pow(*ntime*1.5,(double)1/3);
    if(*winc<4*step) *winc=2*step;// if the length is too small
    if(*winc>=*ntime) *winc=*ntime-step;} // if the length is too big
  //set the spatial-temporal windows:
  int wint = (int) *winc;
  if(*winstp==0) *winstp=wint; //defualt for the forward step:the minimum distance
  //else *winstp=*winstp*wint;   //otherwise a proportion of the temporal window
  sublagt=(double *) R_alloc(wint, sizeof(double));
  sublagt[0]=step;
  sdata=(double *) R_alloc(*ncoord*wint, sizeof(double));
  //set the temporal distances for the sub-sample:
  for(i=0;i<wint;i++){sublagt[i+1]=sublagt[i]+step;nstime++;}
  nsub=floor(((*ntime-wint)/winstp[0]+1));

//start the sub-sampling procedure:
  for(n=0;n<*nrep;n++){
    for(h=0;h<nvari;h++)subvari[h]=0;//initialize the  variance on the subwindows
      for(k=0;k<nsub;k++){//loop for the number of sub-sampling:
      // set the sub-sample of the data:
	SetSampling_st(data,sdata,ncoord,ntime,wint,k,n,nrep);
	nwpair=0;
	for(h=0;h<*npar;h++){sumgrad[h]=0;gradient[h]=0;} //initialize the gradient of the sub-window
	for(h=0;h<*nparc;h++) gradcor[h]=0;
	for(i=0;i<*ncoord;i++){
	  for(t=0;t<nstime;t++){
	    for(j=i;j<*ncoord;j++){
	      if(i==j){
	      for(v=t+1;v<nstime;v++){
	         lagt=fabs(sublagt[t]-sublagt[v]);
	         if(lagt<=*maxtime){
		   nwpair++;
		   rho=CorFct(cormod,0,lagt,parcor);GradCorrFct(rho,cormod,eps,flagcor,gradcor,0,lagt,parcor);
		   switch(*like){
		   case 1:// Conditional likelihood:
		     switch(*type){
		     case 2:
		      Grad_Cond_Gauss(rho,flagnuis,gradcor,gradient,npar,nuis,
				      sdata[(t+nstime*i)+n*nstime* *ncoord],
				      sdata[(v+nstime*j)+n*nstime* *ncoord]);
		      break;}
		     break;
		   case 3: // Marginal likelihood:
		     switch(*type){
		     case 1:// Gradient of the log difference likelihood
		       Grad_Diff_Gauss(rho,flagnuis,gradcor,gradient,npar,nuis,
				       sdata[(t+nstime*i)+n*nstime* *ncoord],
				       sdata[(v+nstime*j)+n*nstime* *ncoord]);
		       break;
		     case 2:// Gradient of the log pairwise likelihood
		       Grad_Pair_Gauss(rho,flagnuis,gradcor,gradient,npar,nuis,
				       sdata[(t+nstime*i)+n*nstime* *ncoord],
				       sdata[(v+nstime*j)+n*nstime* *ncoord]);
		       break;}
		     break;}
		   //sum the gradiend in the subwindow
		   for(h=0;h<*npar;h++) sumgrad[h]=sumgrad[h]+gradient[h];}}}
	      else{
		for(v=t+1;v<nstime;v++){
		  lagt=fabs(sublagt[t]-sublagt[v]);
		  if(lagt<=*maxtime && mlags[i][j]<=*maxdist){
		    nwpair++;
		    rho=CorFct(cormod,mlags[i][j],lagt,parcor);
		    GradCorrFct(rho,cormod,eps,flagcor,gradcor,mlags[i][j],lagt,parcor);
		    switch(*like){
		    case 1:// Conditional likelihood:
		      switch(*type){
		      case 2:
			Grad_Cond_Gauss(rho,flagnuis,gradcor,gradient,npar,nuis,
					sdata[(t+nstime*i)+n*nstime* *ncoord],
					sdata[(v+nstime*j)+n*nstime* *ncoord]);
			break;}
		      break;
		    case 3:// Marginal likelihood:
		      switch(*type){
		      case 1:// Gradient of the log difference likelihood
			Grad_Diff_Gauss(rho,flagnuis,gradcor,gradient,npar,nuis,
					sdata[(t+nstime*i)+n*nstime* *ncoord],
					sdata[(v+nstime*j)+n*nstime* *ncoord]);
			break;
		      case 2:// Gradient of the log pairwise likelihood
			Grad_Pair_Gauss(rho,flagnuis,gradcor,gradient,npar,nuis,
					sdata[(t+nstime*i)+n*nstime* *ncoord],
					sdata[(v+nstime*j)+n*nstime* *ncoord]);
			break;}
		      break;}
		     // sum the gradient in the subwindow
		    for(h=0;h<*npar;h++) sumgrad[h]=sumgrad[h]+gradient[h];}}}}}}
	h=0;
	for(p=0;p<*npar;p++)//update the sub-variance in the subwindow
	  for(q=p;q<*npar;q++){subvari[h]=subvari[h]+sumgrad[p]*sumgrad[q]/nwpair;h++;}}}
  for(h=0;h<nvari;h++) varimat[h]= np[0]*subvari[h]/(nsub * *nrep);//update variability matrix
  return;
}
