#include "header.h"

// Composite conditional log-likelihood for the spatial Gaussian model:
void Comp_Cond_Gauss(int *cormod, double *data, double *nuis, double *par, double *thr, double *res)
{
  int i=0, h=0, j=0, n=0;
  double s1=0.0, s12=0.0;
  double det=0.0, u=0.0, u2=0.0, v=0.0, v2=0.0;
  // Checks the validity of the the nuisance and correlation parameters (nugget, sill and corr):
  if(nuis[1]<0 || nuis[2]<=0 || CheckCor(cormod,par)==-2){
    *res=LOW; return;}
  // Set nuisance parameters:
  s1=nuis[1]+nuis[2];//set nugget + sill
  // Computes the log-likelihood:
  for(i=0; i<(*ncoord-1);i++)
    for(j=(i+1); j<*ncoord;j++){
      // Pairwise distances
      if(lags[h]<=*maxdist){
	s12=nuis[2]*CorFct(cormod, lags[h], 0, par); //sill * corr
	det=pow(s1,2)-pow(s12,2);
	for(n=0;n<*nrep;n++){
	  u=data[(n+i * *nrep)]-nuis[0]; //data[si] - mean
	  v=data[(n+j * *nrep)]-nuis[0]; //data[sj] - mean
	  u2=pow(u,2);
	  v2=pow(v,2);
	  *res+= -log(2*M_PI)-log(det)+log(s1)+
	    (u2+v2)*(0.5/s1-s1/det)+2*s12*u*v/det;}}
      h++;}
  // Checks the return values
  if(!R_FINITE(*res))
    *res = LOW;
  return;
}
// Composite conditional log-likelihood for the spatial-temporal Gaussian model:
void Comp_Cond_Gauss_st(int *cormod, double *data, double *nuis, double *par, double *thr, double *res)
{
  int i=0, j=0, n=0, t=0, v=0;
  double s1=0.0, s12=0.0;
  double det=0.0, u=0.0, u2=0.0, w=0.0, w2=0.0;
  // Checks the validity of the nuisance and correlation parameters (nugget, sill and corr):
  if(nuis[1]<0 || nuis[2]<=0 || CheckCor(cormod,par)==-2){
    *res=LOW; return;}
  // Set nuisance parameters:
  s1=nuis[1]+nuis[2];
 // Computes the log-likelihood:
 for(i=0;i<*ncoord;i++){
   for(t=0;t<*ntime;t++){
     for(j=i;j<*ncoord;j++){
       if(i==j){// marginal temporal log-likelihood:
	 for(v=t+1;v<*ntime;v++){
	   if(mlagt[t][v]<=*maxtime){
	     s12=nuis[2]*CorFct(cormod,0,mlagt[t][v],par);
	     det=pow(s1,2)-pow(s12,2);
	     for(n=0;n<*nrep;n++){
	       u=data[(t+*ntime*i)+n* *nrep]-nuis[0];
	       w=data[(v+*ntime*j)+n* *nrep]-nuis[0];
	       u2=pow(u,2);
	       w2=pow(w,2);
	       *res+= -log(2*M_PI)-log(det)+log(s1)+
		 (u2+w2)*(0.5/s1-s1/det)+2*s12*u*w/det;}}}}
       else{
	 for(v=0;v<*ntime;v++){// spatial-temporal log-likelihood:
	   if(mlagt[t][v]<=*maxtime && mlags[i][j]<=*maxdist){
	     s12=nuis[2]*CorFct(cormod,mlags[i][j],mlagt[t][v],par);
	     det=pow(s1,2)-pow(s12,2);
	     for(n=0;n<*nrep;n++){
	       u=data[(t+*ntime*i)+n* *nrep]-nuis[0];
	       w=data[(v+*ntime*j)+n* *nrep]-nuis[0];
	       u2=pow(u,2);
	       w2=pow(w,2);
	       *res+= -log(2*M_PI)-log(det)+log(s1)+
		 (u2+w2)*(0.5/s1-s1/det)+2*s12*u*w/det;}}}}
     }
   }
 }
  if(!R_FINITE(*res))
    *res = LOW;
  return;
}
// Composite conditional log-likelihood for the spatial Binary Gaussian model:
void Comp_Cond_BinGauss(int *cormod, double *data, double *nuis, double *par, double *thr, double *res)
{
  int i=0, h=0, j=0, n=0;
  double u=0.0, v=0.0;
  double psm=0.0;//probability of marginal success
  double psj=0.0;//probability of joint success
  // Checks the validity of the nuisance and correlation parameters (nugget, sill and corr):
  if(nuis[1]<0 || nuis[2]<=0 || nuis[2]>1 || CheckCor(cormod,par)==-2){
    *res=LOW; return;}
  nuis[1]=1-nuis[2];// define the nugget
  //compute the composite-likelihood:
  for(i=0; i<(*ncoord-1);i++)
    for(j=(i+1); j<*ncoord;j++){
      if(lags[h]<=*maxdist){
	psj=pbnorm(cormod,lags[h],0,nuis,par,*thr);
	psm=pnorm((nuis[0]-*thr)/sqrt(nuis[2]+nuis[1]),0,1,1,0);
	for(n=0;n<*nrep;n++){
	  u=data[(n+i * *nrep)];v=data[(n+j * *nrep)];
	  *res+=2*(((u*v)*log(psj)+(1-u)*(1-v)*log(1-2*psm+psj)+(u*(1-v)+(1-u)*v)*log(psm-psj)))-
	    ((u+v)*log(psm)+log(1-psm)*(2-u-v));}}
      h++;}
  return;
}
// Composite conditional log-likelihood for the binary spatial-temporal Gaussian model:
void Comp_Cond_BinGauss_st( int *cormod, double *data, double *nuis, double *par,double *thr, double *res)
{
  int i=0, j=0, n=0, t=0, v=0;
  double u=0.0, w=0.0;
  double psm=0.0;//probability of marginal success
  double psj=0.0;//probability of joint success
  // Checks the validity of the nuisance and correlation parameters (nugget, sill and corr)
  if(nuis[1]<0 || nuis[2]<=0 || nuis[2]>1 || CheckCor(cormod,par)==-2){
    *res=LOW; return;}
  nuis[1]=1-nuis[2];// define the nugget
  // Computes the log-likelihood:
  for(i=0;i<*ncoord;i++){
    for(t=0;t<*ntime;t++){
      for(j=i;j<*ncoord;j++){
	if(i==j){
	  for(v=t+1;v<*ntime;v++){
	    if(mlagt[t][v]<=*maxtime){
	      psj=pbnorm(cormod,0,mlagt[t][v],nuis,par,*thr);
	      psm=pnorm((nuis[0]-*thr)/sqrt(nuis[2]+nuis[1]),0,1,1,0);
	      for(n=0;n<*nrep;n++){
		u=data[(t+*ntime*i)+n* *nrep];w=data[(v+*ntime*j)+n* *nrep];
		*res+=2*((u*w)*log(psj)+(1-u)*(1-w)*log(1-2*psm+psj)+(u*(1-w)+(1-u)*w)*log(psm-psj))-
		  ((u+w)*log(psm)+log(1-psm)*(2-u-w));}}}}
	else{
	  for(v=0;v<*ntime;v++){
	    if(mlagt[t][v]<=*maxtime && mlags[i][j]<=*maxdist){
	      psj=pbnorm(cormod,mlags[i][j],mlagt[t][v],nuis,par,*thr);
	      psm=pnorm((nuis[0]-*thr)/sqrt(nuis[2]+nuis[1]),0,1,1,0);
	      for(n=0;n<*nrep;n++){
		u=data[(t+*ntime*i)+n* *nrep];
		w=data[(v+*ntime*j)+n* *nrep];
		*res+=2*((u*w)*log(psj)+(1-u)*(1-w)*log(1-2*psm+psj)+
			 (u*(1-w)+(1-u)*w)*log(psm-psj))-
		  ((u+w)*log(psm)+log(1-psm)*(2-u-w));}}}}}}}
 if(!R_FINITE(*res))
    *res = LOW;
 return;
}
// Composite marginal (difference) log-likelihood for the spatial Gaussian model:
void Comp_Diff_Gauss(int *cormod, double *data, double *nuis, double *par, double *thr, double *res)
{
  int i=0, h=0, j=0, n=0;
  double vario=0.0;
  // Checks the validity of the nuisance and correlation parameters (nugget, sill and corr):
  if(nuis[1]<0 || nuis[2]<=0 || CheckCor(cormod,par)==-2){
    *res=LOW; return;}
  // Computes the log-likelihood:
  for(i=0; i<(*ncoord-1);i++)
    for(j=(i+1);j<*ncoord;j++){
      // Pairwise distances
      if(lags[h]<=*maxdist){
	// Variogram: nugget+sill*(1-corr)
	vario=Variogram(cormod,lags[h],0,nuis,par);
	for(n=0;n<*nrep;n++)
	  *res+= -0.5*(log(2*M_PI)+log(vario)+
		       pow(data[(n+i * *nrep)]-
			   data[(n+j * *nrep)],2)/(2*vario));}
      h++;}
  // Checks the return values
  if(!R_FINITE(*res))
    *res = LOW;
  return;
}
// Composite marginal (difference) log-likelihood for the spatial-temporal Gaussian model:
void Comp_Diff_Gauss_st(int *cormod, double *data, double *nuis, double *par, double *thr, double *res)
{
  int i=0, j=0, n=0, t=0, v=0;
  double vario=0.0;
  //Checks the validity of the nuisance and correlation parameters (nugget, sill and corr):
  if(nuis[1]<0 || nuis[2]<=0 || CheckCor(cormod,par)==-2){
    *res=LOW; return;}
  // Computes the log-likelihood:
  for(i=0;i<*ncoord;i++){
    for(t=0;t<*ntime;t++){
      for(j=i;j<*ncoord;j++){
	if(i==j){// marginal temporal log-likelihood:
	  for(v=t+1;v<*ntime;v++){
	    if(mlagt[t][v]<=*maxtime){
	      vario=Variogram(cormod,0,mlagt[t][v],nuis,par);
	      for(n=0;n<*nrep;n++)
		*res+= -0.5*(log(2*M_PI)+log(vario)+
			     pow(data[(t+*ntime*i)]-
				 data[(v+*ntime*j)],2)/(2*vario));}}}
	else{
	  for(v=0;v<*ntime;v++){
	    if(mlagt[t][v]<=*maxtime && mlags[i][j]<=*maxdist){
	      vario=Variogram(cormod,mlags[i][j],mlagt[t][v],nuis,par);
	      for(n=0;n<*nrep;n++)
		*res+= -0.5*(log(2*M_PI)+log(vario)+
			     pow(data[(t+*ntime*i)]-
				 data[(v+*ntime*j)],2)/(2*vario));}}}
      }
    }
  }
 if(!R_FINITE(*res))
    *res = LOW;
}
// Composite marginal (difference) log-likelihood for the binary spatial Gaussian model:
void Comp_Diff_BinGauss( int *cormod, double *data, double *nuis, double *par, double *thr,double *res)
{
  int i=0, h=0, j=0, n=0;
  double diff=0.0;
  double psm=0.0;//probability of marginal success
  double psj=0.0;//probability of joint success
  // Checks the validity of the nuisance and correlation parameters (nugget, sill and corr):
  if(nuis[1]<0 || nuis[2]<=0 || nuis[2]>1 || CheckCor(cormod,par)==-2){
    *res=LOW; return;}
  nuis[1]=1-nuis[2];// define the nugget
  // Computes the compostite log-likelihood:
  for(i=0; i<(*ncoord-1);i++)
    for(j=(i+1); j<*ncoord;j++){
      if(lags[h]<=*maxdist){
	psj=pbnorm(cormod,lags[h],0,nuis,par,*thr);
	psm=pnorm((nuis[0]-*thr)/sqrt(nuis[2]+nuis[1]),0,1,1,0);
	for(n=0;n<*nrep;n++){
	  diff=data[(n+i * *nrep)]-data[(n+j * *nrep)];
	  *res+= (1-R_pow(diff,2))*log(1-2*(psm-psj))+R_pow(diff,2)*log(psm-psj);}}
      h++;}
  return;
}
// Composite marginal (difference) log-likelihood for binary the spatial-temporal Gaussian model:
void Comp_Diff_BinGauss_st(int *cormod, double *data, double *nuis, double *par,double *thr, double *res)
{
  int i=0, j=0, n=0, t=0, v=0;
  double diff=0.0;
  double psm=0.0;//probability of marginal success
  double psj=0.0;//probability of joint success
  // Checks the validity of the nuisance and correlation parameters (nugget, sill and corr)
  if(nuis[1]<0 || nuis[2]<=0 || nuis[2]>1 || CheckCor(cormod,par)==-2){
    *res=LOW; return;}
  nuis[1]=1-nuis[2];// define the nugget
  // Computes the log-likelihood:
  for(i=0;i<*ncoord;i++){
    for(t=0;t<*ntime;t++){
      for(j=i;j<*ncoord;j++){
	if(i==j){
	  for(v=t+1;v<*ntime;v++){
	    if(mlagt[t][v]<=*maxtime){
	      psj=pbnorm(cormod,0,mlagt[t][v],nuis,par,*thr);
	      psm=pnorm((nuis[0]-*thr)/sqrt(nuis[2]+nuis[1]),0,1,1,0);
	      for(n=0;n<*nrep;n++){
		diff=data[(t+*ntime*i)+n* *nrep]-data[(v+*ntime*j)+n* *nrep];
		*res+= (1-R_pow(diff,2))*log(1-2*(psm-psj))+R_pow(diff,2)*log(psm-psj);}}}}
	else{
	  for(v=0;v<*ntime;v++){
	    if(mlagt[t][v]<=*maxtime && mlags[i][j]<=*maxdist){
	      psj=pbnorm(cormod,mlags[i][j],mlagt[t][v],nuis,par,*thr);
	      psm=pnorm((nuis[0]-*thr)/sqrt(nuis[2]+nuis[1]),0,1,1,0);
	      for(n=0;n<*nrep;n++){
		diff=data[(t+*ntime*i)+n* *nrep]-data[(v+*ntime*j)+n* *nrep];
		*res+=(1-R_pow(diff,2))*log(1-2*(psm-psj))+R_pow(diff,2)*log(psm-psj);}}}}
}}}
 if(!R_FINITE(*res))
    *res = LOW;
 return;
}
// Composite marginal (pariwise) log-likelihood for the spatial Gaussian model:
void Comp_Pair_Gauss(int *cormod, double *data, double *nuis, double *par, double *thr, double *res)
{
  int i=0, h=0, j=0, n=0;
  double s1=0.0, s12=0.0;
  double det=0.0, u=0.0, u2=0.0, v=0.0, v2=0.0;
  // Checks the validity of the nuisance and correlation parameters (nugget, sill and corr):
  if(nuis[1]<0 || nuis[2]<=0 || CheckCor(cormod,par)==-2){
    *res=LOW; return;}
  // Set nuisance parameters:
  s1=nuis[1]+nuis[2];//set nugget + sill
  // Computes the log-likelihood:
  for(i=0;i<(*ncoord-1);i++)
    for(j=(i+1); j<*ncoord;j++){
      // Pairwise distances
      if(lags[h]<=*maxdist){
	s12=nuis[2]*CorFct(cormod,lags[h],0,par); //sill * corr
	det=pow(s1,2)-pow(s12,2);
	for(n=0;n<*nrep;n++){
	  u=data[(n+i * *nrep)]-nuis[0]; //data[si] - mean
	  v=data[(n+j * *nrep)]-nuis[0]; //data[sj] - mean
	  u2=R_pow(u,2);
	  v2=R_pow(v,2);
	  *res+= -0.5*(2*log(2*M_PI)+log(det)+
		       (s1*(u2+v2)-2*s12*u*v)/det);}}
      h++;}
  // Checks the return values
  if(!R_FINITE(*res))
    *res = LOW;
  return;
}
// Composite marginal (pariwise) log-likelihood for the spatial-temporal Gaussian model:
void Comp_Pair_Gauss_st(int *cormod, double *data, double *nuis, double *par, double *thr, double *res)
{
  int i=0, j=0, n=0, t=0, v=0;
  double s1=0.0, s12=0.0;
  double det=0.0, u=0.0, u2=0.0, w=0.0, w2=0.0;
  // Checks the validity of the nuisance and correlation parameters (nugget, sill and corr)
  if(nuis[1]<0 || nuis[2]<=0 || CheckCor(cormod,par)==-2){
    *res=LOW; return;}
  // Set nuisance parameters:
  s1=nuis[1]+nuis[2];
  // Computes the log-likelihood:
  for(i=0;i<*ncoord;i++){
    for(t=0;t<*ntime;t++){
      for(j=i;j<*ncoord;j++){
	if(i==j){
	  for(v=t+1;v<*ntime;v++){
	    if(mlagt[t][v]<=*maxtime){
	      s12=nuis[2]*CorFct(cormod,0, mlagt[t][v],par);
	      det=pow(s1,2)-pow(s12,2);
	      for(n=0;n<*nrep;n++){
		u=data[(t+*ntime*i)+n* *nrep]-nuis[0];
		w=data[(v+*ntime*j)+n* *nrep]-nuis[0];
		u2=pow(u,2);
		w2=pow(w,2);
		*res+= -0.5*(2*log(2*M_PI)+log(det)+
			     (s1*(u2+w2)-2*s12*u*w)/det);}}}}
	else{
	  for(v=0;v<*ntime;v++){
	    if(mlagt[t][v]<=*maxtime && mlags[i][j]<=*maxdist){
	      s12=nuis[2]*CorFct(cormod,mlags[i][j], mlagt[t][v],par);
	      det=pow(s1,2)-pow(s12,2);
	      for(n=0;n<*nrep;n++){
		u=data[(t+*ntime*i)+n* *nrep]-nuis[0];
		w=data[(v+*ntime*j)+n* *nrep]-nuis[0];
		u2=pow(u,2);
		w2=pow(w,2);
		*res+= -0.5*(2*log(2*M_PI)+log(det)+
			     (s1*(u2+w2)-2*s12*u*w)/det);}}}}
      }
    }
  }
 if(!R_FINITE(*res))
    *res = LOW;
}
// Composite marginal pairwise log-likelihood for the binary spatial Gaussian model:
void Comp_Pair_BinGauss( int *cormod, double *data, double *nuis, double *par, double *thr,double *res)
{
  int i=0, h=0, j=0, n=0;
  double u=0.0, v=0.0;
  double psm=0.0;//probability of marginal success
  double psj=0.0;//probability of joint success
  // Checks the validity of the nuisance and correlation parameters (nugget, sill and corr):
  if(nuis[1]<0 || nuis[2]<=0 || nuis[2]>1 || CheckCor(cormod,par)==-2){
    *res=LOW; return;}
  nuis[1]=1-nuis[2];// define the nugget
  //compute the composite log-likelihood:
  for(i=0; i<(*ncoord-1);i++)
    for(j=(i+1); j<*ncoord;j++){
      if(lags[h]<=*maxdist){
	psj=pbnorm(cormod,lags[h],0,nuis,par,*thr);
	psm=pnorm((nuis[0]-*thr)/sqrt(nuis[2]+nuis[1]),0,1,1,0);
	for(n=0;n<*nrep;n++){
	  u=data[(n+i * *nrep)];v=data[(n+j * *nrep)];
	  *res+=((u*v)*log(psj)+(1-u)*(1-v)*log(1-2*psm+psj)+(u*(1-v)+(1-u)*v)*log(psm-psj));}}
      h++;}
  return;
}
// Composite marginal (pariwise) log-likelihood for the binary spatial-temporal Gaussian model:
void Comp_Pair_BinGauss_st(int *cormod, double *data, double *nuis, double *par, double *thr,double *res)
{
  int i=0, j=0, n=0, t=0, v=0;
  double u=0.0, w=0.0;
  double psm=0.0;//probability of marginal success
  double psj=0.0;//probability of joint success
  // Checks the validity of the nuisance and correlation parameters (nugget, sill and corr)
  if(nuis[1]<0 || nuis[2]<=0 || nuis[2]>1 || CheckCor(cormod,par)==-2){
    *res=LOW; return;}
  nuis[1]=1-nuis[2];// define the nugget
  // Computes the log-likelihood:
  for(i=0;i<*ncoord;i++){
    for(t=0;t<*ntime;t++){
      for(j=i;j<*ncoord;j++){
	if(i==j){
	  for(v=t+1;v<*ntime;v++){
	    if(mlagt[t][v]<=*maxtime){
	      psj=pbnorm(cormod,0,mlagt[t][v],nuis,par,*thr);
	      psm=pnorm((nuis[0]-*thr)/sqrt(nuis[2]+nuis[1]),0,1,1,0);
	      for(n=0;n<*nrep;n++){
		u=data[(t+*ntime*i)+n* *nrep];w=data[(v+*ntime*j)+n* *nrep];
		*res+=((u*w)*log(psj)+(1-u)*(1-w)*log(1-2*psm+psj)+(u*(1-w)+(1-u)*w)*log(psm-psj));}}}}
	else{
	  for(v=0;v<*ntime;v++){
	    if(mlagt[t][v]<=*maxtime && mlags[i][j]<=*maxdist){
	      psj=pbnorm(cormod,mlags[i][j],mlagt[t][v],nuis,par,*thr);
	      psm=pnorm((nuis[0]-*thr)/sqrt(nuis[2]+nuis[1]),0,1,1,0);
	      for(n=0;n<*nrep;n++){
		u=data[(t+*ntime*i)+n* *nrep]-nuis[0];w=data[(v+*ntime*j)+n* *nrep]-nuis[0];
		*res+=((u*w)*log(psj)+(1-u)*(1-w)*log(1-2*psm+psj)+(u*(1-w)+(1-u)*w)*log(psm-psj));}}}}
}}}
 if(!R_FINITE(*res))
    *res = LOW;
 return;
}
// Composite pariwise log-likelihood for the max-stable spatial extremal Gaussian model:
void Comp_Ext_Gauss(int *cormod, double *data, double *nuis, double *par, double *thr, double *res)
{
  int i=0, h=0, j=0, n=0;
  double a=0.0, rho=0.0, d2V=0.0, dxV=0.0, dyV=0.0;
  double x=0.0, x2=0.0, xy2=0, y=0.0, y2=0.0, V=0.0;
  // Checks the validity of the nuisance and correlation parameters (nugget, sill and corr):
  if(nuis[0]<=0 || nuis[0]>1 || CheckCor(cormod,par)==-2){
    *res=LOW; return;}
  // Computes the log-likelihood:
  for(i=0;i<(*ncoord-1);i++)
    for(j=(i+1); j<*ncoord;j++){
      // Pairwise distances
      if(lags[h]<=*maxdist){
	rho=nuis[0]*CorFct(cormod,lags[h],0,par);//rho=sill*corr
	if(rho>.99999996){
	  for(n=0;n<*nrep;n++){
	    x=data[(n+i * *nrep)]; //data[si]
	    y=data[(n+j * *nrep)]; //data[sj]
	    if(x>=y) *res+=-2*log(y)-1/y;
	    else *res+=-2*log(x)-1/x;}}
	else{
	  for(n=0;n<*nrep;n++){
	    x=data[(n+i * *nrep)]; //data[si]
	    y=data[(n+j * *nrep)]; //data[sj]
            x2=x*x;//pow(u,2);
	    y2=y*y;//pow(v,2);
	    xy2=2*x*y;
	    a=sqrt(x2+y2-xy2*rho);//sqrt of the quadratic form
	    V=-(x+y+a)/xy2;
	    dxV=-0.5*(x*rho-a-y)/(x2*a);
	    dyV=-0.5*(y*rho-a-x)/(y2*a);
	    d2V=0.5*(1-rho*rho)/(a*a*a);
	    *res+=V+log(d2V+dxV*dyV);}}}
      h++;}
  // Checks the return values
  if(!R_FINITE(*res))
    *res = LOW;
  return;
}
// Composite pariwise log-likelihood for the max-stable spatial extremal-t model:
void Comp_Ext_T(int *cormod, double *data, double *nuis, double *par, double *thr, double *res)
{
  int i=0, h=0, j=0, n=0;
  double a=0.0, ac=0.0, aci=0.0, c=0.0, df=nuis[0], df1=nuis[0]+1;
  double d1tx=0.0, d1ty=0.0, dtx=0.0, dty=0.0, d2V=0.0, dxV=0.0, dyV=0.0;
  double opdf=1+1/nuis[0], ptx=0.0, pty=0.0, rho=0.0, x=0.0, x2=0.0;
  double x2d=0.0, xyd=0.0, y=0.0, y2=0.0, y2d=0.0, V=0.0, w=0.0, z=0.0;
  // Checks the validity of the nuisance and correlation parameters (nugget, sill and corr):
  if(df<=0 || nuis[1]<=0 || nuis[1]>1 || CheckCor(cormod,par)==-2){
    *res=LOW; return;}
  // Computes the log-likelihood:
  for(i=0;i<(*ncoord-1);i++)
    for(j=(i+1); j<*ncoord;j++){
      // Pairwise distances
      if(lags[h]<=*maxdist){
	rho=nuis[1]*CorFct(cormod,lags[h],0,par);//rho=sill*corr
	a=sqrt(df1/(1-pow(rho, 2)));
	for(n=0;n<*nrep;n++){
	  x=data[(n+i * *nrep)]; //data[si]
	  y=data[(n+j * *nrep)]; //data[sj]
	  c=pow(y/x,1/df);
	  ac=a*c;
	  aci=a/c;
	  z=(c-rho)*a;
	  w=(1/c-rho)*a;
	  x2=pow(x,2);
	  y2=pow(y,2);
	  x2d=x2*df;
	  y2d=y2*df;
	  xyd=x*y*df;
	  ptx=pt(z,df1,1,0);
	  pty=pt(w,df1,1,0);
	  dtx=dt(z,df1,0);
	  dty=dt(w,df1,0);
	  d1tx=d1x_dt(z,df1);
	  d1ty=d1x_dt(w,df1);
	  //defines the log-likelihood components:
	  V=-ptx/x-pty/y;
	  dxV=ptx/x2+dtx*ac/x2d-dty*aci/xyd;
	  dyV=pty/y2+dty*aci/y2d-dtx*ac/xyd;
	  d2V=ac*(dtx*opdf+d1tx*ac/df)/x2d/y+
	    aci*(dty*opdf+d1ty*aci/df)/y2d/x;
	  *res+=V+log(d2V+dxV*dyV);}}//log-likelihood
      h++;}
  // Checks the return values
  if(!R_FINITE(*res))
    *res = LOW;
  return;
}
// Pairwise log-likelihood for the max-stable spatial Brown-Resnick model:
void Comp_Brow_Resn(int *cormod, double *data, double *nuis, double *par, double *thr, double *res)
{
  int i=0, j=0, h=0, n=0;
  double a=0.0, ao2=0.0, ax2=0.0, axy=0.0, ay2=0.0, d2V=0.0;
  double dx=0.0, dxV=0.0, dy=0.0, dyV=0.0, lyx=0.0, px=0.0;
  double py=0.0, vario=0.0, V=0.0, x=0.0, y=0.0, x2=0.0;
  double y2=0.0, w=0.0, z=0.0;
  // Checks the validity of the variogram parameters:
  if(CheckCor(cormod,par)==-2){
    *res=LOW; return;}
  // Computes the log-likelihood:
  for(i=0;i<(*ncoord-1);i++)
    for(j=(i+1); j<*ncoord;j++){
      // Pairwise distances
      if(lags[h]<=*maxdist){
	vario=VarioFct(cormod,lags[h],par);//vario~log(n)(1-rho(n))
	a=sqrt(vario);// Husler-Reiss coefficient (lambda)
	ao2=0.5*a;
	for(n=0;n<*nrep;n++){
	  // Defines useful componets:
	  x=data[(n+i * *nrep)]; //data[si];
	  y=data[(n+j * *nrep)]; //data[sj];
	  axy=a*x*y;
	  x2=pow(x,2);
	  y2=pow(y,2);
	  ax2=a*x2;
	  ay2=a*y2;
	  lyx=log(y/x)/a;
	  z=ao2+lyx;
	  w=ao2-lyx;
	  px=pnorm(z,0,1,1,0);
	  py=pnorm(w,0,1,1,0);
	  dx=dnorm(z,0,1,0);
	  dy=dnorm(w,0,1,0);
	  // Defines the likelihood components:
	  V=-px/x-py/y;
	  dxV=px/x2+dx/ax2-dy/axy;
	  dyV=py/y2+dy/ay2-dx/axy;
	  d2V=(w*dx*y+z*dy*x)/(ax2*ay2);
	  // defines the pairwise loglikelihood:
	  *res+=V+log(d2V+dxV*dyV);}}
      h++;}
  // Checks the return values
  if(!R_FINITE(*res))
    *res = LOW;
  return;
}

