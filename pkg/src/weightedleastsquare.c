#include "header.h"

// binned spatial variogram:
void Binned_Variogram(double *bins, double *data, int *lbins, double *moms, int *nbins)
{
  int h=0, i=0, j=0, n=0, p=0;
  double step=0.0;
  //Set the binnes step:
  step=(*maxdist-*minimdista)/(*nbins-1);
  bins[0]= *minimdista;
  //define bins:
  for(h=1;h<*nbins;h++)
    bins[h]=bins[h-1]+step;
  //Computes the binned moments:
  for(i=0;i<(*ncoord-1);i++)
    for(j=(i+1);j<*ncoord;j++){
      if(lags[p]<=*maxdist){
	for(h=0;h<(*nbins-1);h++)
	  if((bins[h]<=lags[p]) && (lags[p]<bins[h+1])){
	    for(n=0;n<*nrep;n++){
	      moms[h]+=0.5*pow(data[n+i * *nrep]-data[n+j * *nrep],2);
	      lbins[h]+=1;}}}
      p++;}
  return;
}
// binned spatial-temporal variogram:
void Binned_Variogram_st(double *bins, double *bint, double *data, int *lbins, int *lbinst,
			 int *lbint, double *moms,double *momst, double *momt, int *nbins, int *nbint)
{
int h=0, i=0, j=0, n=0;
  int q=0, t=0, u=0, v=0;
  double step=0.0;
  //defines the spatial bins:
  step=(*maxdist-*minimdista)/(*nbins-1);
  bins[0]= *minimdista;
  for(h=1;h<*nbins;h++)
    bins[h]=bins[h-1]+step;
  //defines the temporal bins:
  bint[0]=*minimtime;
  for(h=1;h<*nbint;h++)
    bint[h]=bint[h-1]+bint[0];
  //computes the empirical variogram:
  for(i=0;i<*ncoord;i++)
    for(t=0;t<*ntime;t++)
      for(j=i;j<*ncoord;j++){
	if(i==j)//computes the marignal temporal variogram:
	  for(v=(t+1);v<*ntime;v++){
	    if(mlagt[t][v]<=*maxtime)
	      for(u=0;u<*nbint;u++)
		if(bint[u]==mlagt[t][v])
		  for(n=0;n<*nrep;n++){
		    momt[u]+=0.5*pow(data[(t+i * *ntime)+n* *nrep]-data[(v+i * *ntime)+n* *nrep], 2);
		    lbint[u]+=1;}}
	else{
	  for(v=0;v<*ntime;v++){
	    if(t==v){// computes the marginal spatial variogram:
	      if(mlags[i][j]<=*maxdist)
		for(h=0;h<(*nbins-1);h++)
		  if((bins[h]<=mlags[i][j]) && (mlags[i][j]<bins[h+1]))
		    for(n=0;n<*nrep;n++){
		      moms[h]+=0.5*pow(data[(t+i * *ntime)+n* *nrep]-data[(t+j * *ntime)+n* *nrep], 2);
		      lbins[h]+=1;}}
	    if(v>t){// computes the spatial-temporal variogram:
	      if(mlags[i][j]<=*maxdist && mlagt[t][v]<=*maxtime){
		q=0;
		for(h=0;h<(*nbins-1);h++)
		  for(u=0;u<*nbint;u++){
		    if((bins[h]<=mlags[i][j]) && (mlags[i][j]<bins[h+1]) && (bint[u]==mlagt[t][v]))
		      for(n=0;n<*nrep;n++){
			momst[q]+=0.5*pow(data[(t+i * *ntime)+n* *nrep]-data[(v+j * *ntime)+n* *nrep], 2);
			lbinst[q]+=1;}
		    q++;}}}}}}
  return;
}
/*
void Binned_Variogram_st(double *bins, double *bint, double *data, int *lbins,
			 int *lbinst, int *lbint, double *moms, double *momst,
			 double *momt, int *nbins, int *nbint)
{
  int h=0, i=0, j=0, n=0, ndata=*ncoord * *ntime;
  int q=0, t=0, u=0, v=0;
  double step=0.0;
  //defines the spatial bins:
  step=(*maxdist-*minimdista)/(*nbins-1);
  bins[0]= *minimdista;
  for(h=1;h<*nbins;h++)
    bins[h]=bins[h-1]+step;
  //defines the temporal bins:
  bint[0]=*minimtime;
  for(h=1;h<*nbint;h++)
    bint[h]=bint[h-1]+bint[0];
  //loops for the cross sptial-temporal differences
  for(t=0;t<*ntime;t++) // first temporal component
    {
      for(v=t;v<*ntime;v++) // second temporal component
	{
	  for(i=0;i<*ncoord;i++) // first spatial component
	    {
	      for(j=0;j<*ncoord;j++) // second spatial component
		{
		  // START marginal spatial variogram
		  if(t==v && j>i)
		    {
		      if(mlags[i][j]<=*maxdist)
			  {
			    for(h=0;h<(*nbins-1);h++)
			      if((bins[h]<=mlags[i][j]) && (mlags[i][j]<bins[h+1]))
				{
				  for(n=0;n<*nrep;n++)
				    {
				      moms[h]+=0.5*pow(data[(t+i * *ntime)+n* *nrep]-
						       data[(t+j * *ntime)+n* *nrep], 2);
				      lbins[h]+=1;
				    }
				}
			  }
		    }
		  // END marginal spatial variogram
                  // START marginal temporal variogram
		  if(v>t && i==j)
		    {
		      if(mlagt[t][v]<=*maxtime)
			{
			  for(u=0;u<*nbint;u++)
			    if(bint[u]==mlagt[t][v])
			      {
				for(n=0;n<*nrep;n++)
				  {
				    momt[u]+=0.5*pow(data[(t+i * *ntime)+n* *nrep]-
						     data[(v+i * *ntime)+n* *nrep], 2);
				    lbint[u]+=1;
				  }
			      }
			}
		    }
		  // END marignal temporal variogram
		  // START cross spatial-temporal variogram
		  if(v>t && j!=i)
		    {
		      if(mlags[i][j]<=*maxdist && mlagt[t][v]<=*maxtime)
			{
			  q=0;
			  for(h=0;h<(*nbins-1);h++)
			    for(u=0;u<*nbint;u++)
			      {
				if((bins[h]<=mlags[i][j]) && (mlags[i][j]<bins[h+1]) && (bint[u]==mlagt[t][v]))
				  {
				    for(n=0;n<*nrep;n++)
				      {
					momst[q]+=0.5*pow(data[(t+i * *ntime)+n* *nrep]-
							  data[(v+j * *ntime)+n* *nrep], 2);
					lbinst[q]+=1;
				      }
				  }
				q++;
			      }
			}
		    }
		  // END cross spatial-temporal variogram
		}
	    }
	}
    }
  return;
  }*/
// binned madogram:
void Binned_Madogram(double *bins, double *data, int *lbins, double *moms, int *nbins)
{
  int h=0, i=0, j=0, n=0, p=0;
  double step=0.0;
  //Set the binnes step:
  step=(*maxdist-*minimdista)/(*nbins-1);
  bins[0]=*minimdista;
  //define bins:
  for(h=1;h<*nbins;h++)
    bins[h]=bins[h-1]+step;
  //Computes the binned moments:
  for(i=0;i<(*ncoord-1);i++)
    for(j=(i+1);j<*ncoord;j++){
      if(lags[p]<=*maxdist){
	for(h=0;h<(*nbins-1);h++)
	  if((bins[h]<=lags[p]) && (lags[p]<bins[h+1])){
	    for(n=0; n<*nrep;n++){
	      moms[h]+=0.5*fabs(data[n+i * *nrep]-data[n+j * *nrep]);
	      lbins[h]+=1;}}}
      p++;}
  return;
}
// variogram cloud:
void Cloud_Variogram(double *bins, double *data, int *lbins, double *moms, int *nbins)
{
  int h=0, i=0, j=0, n=0;
 //Computes the cloud moments:
  for(i=0;i<(*ncoord-1);i++)
    for(j=(i+1);j<*ncoord;j++){
      bins[h]=lags[h];
      for(n=0;n<*nrep;n++)
	moms[h]+=0.5*pow(data[n+i * *nrep]-data[n+j * *nrep],2);
      lbins[h]=*nrep;
      h++;}
  return;
}
// madogram cloud:
void Cloud_Madogram(double *bins, double *data, int *lbins, double *moms, int *nbins)
{
  int h=0, i=0, j=0, n=0;
 //Computes the cloud momemts:
  for(i=0;i<(*ncoord-1);i++)
    for(j=(i+1);j<*ncoord;j++){
      bins[h]=lags[h];
      for(n=0;n<*nrep;n++)
	moms[h]+=0.5*fabs(data[n+i * *nrep]-data[n+j * *nrep]);
      lbins[h]=*nrep;
      h++;}
  return;
}
// Least square method for Gaussian model:
/*void LeastSquare_G(double *bins, double *bint, int *cormod, double *lbins, double *moms,
		   int *nbins, int *nbint, double *nuis, double *par, double *res)
{
  int h=0;
  double vario=0.0, varhat=0.0;
  //Checks the nuisance parameters (nugget, sill and corr):
  if(nuis[1]<0 || nuis[2]<=0 || CheckCor(cormod,1,0,par)==-2) {
    *res=LOW; return;}
  // Computes the least squares:
  for(h=0;h<(*nbins-1);h++)
    if(lbins[h]){
      vario=moms[h]/lbins[h];// Computes the empirical variogram
	//variogram=nugget+sill*(1-corr)
      varhat=Variogram(cormod,0.5*(bins[h]+bins[h+1]),0,nuis,par);
      if(varhat)
	*res=*res-pow(vario-varhat,2)*lbins[h];}// Computes the least squares
  return;
  }*/
// Least square method for Gaussian spatial-temporal random field:
void LeastSquare_G(double *bins, double *bint, int *cormod, double *lbins, double *moms,
		      int *nbins, int *nbint, double *nuis, double *par, double *res)
{
  int h=0, i=0, u=0;
  double vario=0.0, varhat=0.0;
  //Checks the nuisance parameters (nugget, sill and corr):
  if(nuis[1]<0 || nuis[2]<=0 || CheckCor(cormod,par)==-2) {
    *res=LOW; return;}
  // Computes the least squares:
  for(u=0;u<*nbint;u++)
    for(h=0;h<*nbins;h++){
      vario=moms[i]/lbins[i];// Computes the empirical variogram
      varhat=Variogram(cormod,0.5*(bins[h]+bins[h+1]),bint[u],nuis,par);
      *res=*res-pow(vario-varhat,2)*lbins[i];// Computes the least squares
      i++;}
  return;
}
// Weighted least square method for Gaussian spatial random field:
/*void WLeastSquare_G(double *bins, double *bint, int *cormod, double *lbins, double *moms,
		    int *nbins, int *nbint, double *nuis, double *par, double *res)
{
  int h=0;
  double vario=0.0, varhat=0.0;
  //Checks the nuisance parameters (nugget, sill and corr):
  if(nuis[1]<0 || nuis[2]<=0 || CheckCor(cormod,1,0,par)==-2) {
    *res=LOW; return;}
  // Computes the least squares:
  for(h=0;h<(*nbins-1);h++)
    if(lbins[h]){
      vario=moms[h]/lbins[h];// Computes the empirical variogram
	//variogram=nugget+sill*(1-corr)
      varhat=Variogram(cormod,0.5*(bins[h]+bins[h+1]),0,nuis,par);
      if(varhat)
	*res=*res-pow(vario/varhat-1,2)*lbins[h];}// Computes the weighted least squares
  return;
  }*/
// Weighted least square method for Gaussian spatial-temporal random field:
void WLeastSquare_G(double *bins, double *bint, int *cormod, double *lbins, double *moms,
		       int *nbins, int *nbint, double *nuis, double *par, double *res)
{
  int h=0, i=0, u=0;
  double vario=0.0, varhat=0.0;
  //Checks the nuisance parameters (nugget, sill and corr):
  if(nuis[1]<0 || nuis[2]<=0 || CheckCor(cormod,par)==-2) {
    *res=LOW; return;}
  // Computes the least squares:
  for(u=0;u<*nbint;u++)
    for(h=0;h<*nbins;h++){
      vario=moms[i]/lbins[i];// Computes the empirical variogram
      varhat=Variogram(cormod,0.5*(bins[h]+bins[h+1]),bint[u],nuis,par);
      *res=*res-pow(vario/varhat-1,2)*lbins[i];// Computes the weighted least squares
      i++;}
  return;
}
// Least square method for max-stable extremal-Gaussian model:
void LeastSquare_MEG(double *bins, double *bint, int *cormod, double *lbins, double *moms,
		     int *nbins, int *nbint, double *nuis, double *par, double *res)
{
  int h=0;
  double vario=0.0, extcoeff=0.0, extcfhat=0.0;
  //Checks the nuisance and correlation parameters (sill, corr):
  if(nuis[0]<=0 || nuis[0]>1 || CheckCor(cormod,par)==-2) {
    *res=LOW; return;}
  // Computes the least squares:
  for(h=0;h<(*nbins-1);h++)
    if(lbins[h]){
      vario=moms[h]/lbins[h];// Computes the variogram
      extcoeff=(1+2*vario)/(1-2*vario);// Computes the extremal coefficient
      // Computes the extremal coefficient for the Gaussian extremal model:
      extcfhat=1+sqrt(0.5*(1-nuis[0]*CorFct(cormod,0.5*(bins[h]+bins[h+1]),0,par)));
      *res=*res-pow(extcoeff-extcfhat,2);}// Computes the least squares
  return;
}

// Weighted least square method for max-stable extremal-Gaussian model:
void WLeastSquare_MEG(double *bins, double *bint, int *cormod, double *lbins, double *moms,
		      int *nbins, int *nbint, double *nuis, double *par, double *res)
{
  int h=0;
  double variogram=0.0, vario=0.0;

  for(h=0;h<(*nbins-1);h++)
    if(lbins[h])
      {
	variogram=moms[h]/lbins[h];
	vario=nuis[0]+nuis[1]*//nugget+sill*(1-corr)
	  (1-CorFct(cormod, (bins[h]+bins[h + 1])/2,0,par));
	if(vario)
	  *res=*res-pow(variogram/vario-1, 2)*lbins[h];
      }
  return;
}

// Least square method for max-stable Brown-Resnick model:
void LeastSquare_MBR(double *bins, double *bint, int *cormod, double *lbins, double *moms,
		     int *nbins, int *nbint, double *nuis, double *par, double *res)
{
  int h=0;
  double vario=0.0, extcoeff=0.0, extcfhat=0.0;
  //Checks the nuisance and correlation parameters (sill, corr):
  if(VarioFct(cormod,1,par)==-2) {
    *res=LOW; return;}
  // Computes the least squares:
  for(h = 0; h < (*nbins - 1); h++)
    if(lbins[h]){
      vario=moms[h]/lbins[h];// Computes the variogram
      extcoeff=(1+2*vario)/(1-2*vario);// Computes the extremal coefficient
      // Computes the extremal coefficient for the Brown-Renick model:
      extcfhat=2*pnorm(0.5*sqrt(VarioFct(cormod,0.5*(bins[h]+bins[h+1]),par)),0,1,1,0);
      *res=*res-pow(extcoeff-extcfhat,2);}// Computes the least squares
  return;
}

// Least square method for max-stable extremal-t model:
void LeastSquare_MET(double *bins, double *bint, int *cormod, double *lbins, double *moms,
		     int *nbins, int *nbint, double *nuis, double *par, double *res)
{
  int h=0;
  double df1=0.0, rho=0.0, vario=0.0, extcoeff=0.0, extcfhat=0.0;
  //Checks the nuisance and correlation parameters (df, sill, corr):
  if(nuis[0]<=0 || nuis[1]<=0 || nuis[1]>1 || CheckCor(cormod,par)==-2) {
    *res=LOW; return;}
  df1=nuis[0]+1;// Set the df of the t cdf:
  // Computes the least squares:
  for(h=0;h<(*nbins-1);h++)
    if(lbins[h]){
      vario=moms[h]/lbins[h];// Computes the variogram
      extcoeff=(1+2*vario)/(1-2*vario);// Computes the extremal coefficient
      rho=nuis[1]*CorFct(cormod,0.5*(bins[h]+bins[h+1]),0,par);
      // Computes the extremal coefficient for the extremal-t model:
      extcfhat=2*pt(sqrt(df1*(1-rho)/(1+rho)),df1,1,0);
      *res=*res-pow(extcoeff-extcfhat,2);}// Computes the least squares
  return;
}


