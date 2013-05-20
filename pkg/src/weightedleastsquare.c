#include "header.h"




// binned spatial lorelogram:
void Binned_Lorelogram(double *bins, double *data, int *lbins, double *moms,int *nbins)
{
  int h=0, i=0, j=0, n=0;
  double step=0.0,*n11,*n10,*n01,*n00;

  n11=(double *) calloc(*nbins-1, sizeof(double));
  n10=(double *) calloc(*nbins-1, sizeof(double));
  n01=(double *) calloc(*nbins-1, sizeof(double));
  n00=(double *) calloc(*nbins-1, sizeof(double));

  //Set the binnes step:
  step=(*maxdist-*minimdista)/(*nbins-1);
  bins[0]= *minimdista;
  //define bins:
  for(h=1;h<*nbins;h++)
    bins[h]=bins[h-1]+step;
  //Computes the binned statistics:
  for(i=0;i<(*ncoord-1);i++)
    for(j=(i+1);j<*ncoord;j++){
      if(mlags[i][j]<=*maxdist){
	for(h=0;h<(*nbins-1);h++)
	  if((bins[h]<=mlags[i][j]) && (mlags[i][j]<bins[h+1])){
	    for(n=0;n<*nrep;n++){
	      if(data[n+i * *nrep] && data[n+j * *nrep]) n11[h]++;
	      if(data[n+i * *nrep] && !data[n+j * *nrep]) n10[h]++;
	      if(!data[n+i * *nrep] && data[n+j * *nrep]) n01[h]++;
	      if(!data[n+i * *nrep] && !data[n+j * *nrep]) n00[h]++;
	      }}}}
// computing log odds ration in each bin
 for(h=0;h<(*nbins-1);h++){
   if(n11[h]&&n10[h]&&n01[h]&&n00[h]){
     moms[h]=log((n11[h]*n00[h])/(n01[h]*n10[h]));
     lbins[h]=1;}
   else{
     moms[h]=1;
     lbins[h]=0;}}
  return;
}
















// binned spatial-temporal variogram:
void Binned_Lorelogram_st(double *bins, double *bint, double *data, int *lbins, int *lbinst,
			 int *lbint, double *moms,double *momst, double *momt, int *nbins, int *nbint)
{
  int h=0, i=0, j=0, n=0;
  int q=0, t=0, u=0, v=0;
  double step=0.0,*n11s,*n10s,*n01s,*n00s,*n11t;
  double *n10t,*n01t,*n00t,*n11,*n10,*n01,*n00;
  n11s=(double *) calloc((*nbins-1), sizeof(double));
  n10s=(double *) calloc((*nbins-1), sizeof(double));
  n01s=(double *) calloc((*nbins-1), sizeof(double));
  n00s=(double *) calloc((*nbins-1), sizeof(double));
  n11t=(double *) calloc((*nbint), sizeof(double));
  n10t=(double *) calloc((*nbint), sizeof(double));
  n01t=(double *) calloc((*nbint), sizeof(double));
  n00t=(double *) calloc((*nbint), sizeof(double));
  n11=(double *) calloc((*nbins-1)*(*nbint), sizeof(double));
  n10=(double *) calloc((*nbins-1)*(*nbint), sizeof(double));
  n01=(double *) calloc((*nbins-1)*(*nbint), sizeof(double));
  n00=(double *) calloc((*nbins-1)*(*nbint), sizeof(double));
   //defines the spatial bins:
  step=*maxdist/(*nbins-1);
  bins[0]=0;
  for(h=1;h<*nbins;h++)
    bins[h]=bins[h-1]+step;
  //defines the temporal bins:
  bint[0]=0;
  for(u=1;u<*nbint;u++)
    bint[u]=bint[u-1]+minimtime[0];
  //computes the empirical variogram:
  for(i=0;i<*ncoord;i++){
    for(t=0;t<*ntime;t++){
      for(j=i;j<*ncoord;j++){
	if(i==j){//computes the marignal temporal lorelogram:
	  for(v=(t+1);v<*ntime;v++){
	    if(mlagt[t][v]<=*maxtime){
	      for(u=0;u<*nbint;u++){
		if(bint[u]==mlagt[t][v]){

		  for(n=0;n<*nrep;n++){
		    if(data[(t+i * *ntime)+n* *nrep]&&data[(v+i * *ntime)+n* *nrep]) n11t[u]++;
		    if(data[(t+i * *ntime)+n* *nrep]&&!data[(v+i * *ntime)+n* *nrep]) n10t[u]++;
		    if(!data[(t+i * *ntime)+n* *nrep]&&data[(v+i * *ntime)+n* *nrep]) n01t[u]++;
		    if(!data[(t+i * *ntime)+n* *nrep]&&!data[(v+i * *ntime)+n* *nrep]) n00t[u]++;
		    }}}}}}
	else{
	  for(v=0;v<*ntime;v++){
	    if(t==v){// computes the marginal spatial lorelogram:
	      if(mlags[i][j]<=*maxdist){
		for(h=0;h<(*nbins-1);h++){
		  if((bins[h]<=mlags[i][j]) && (mlags[i][j]<bins[h+1])){
		    for(n=0;n<*nrep;n++){
		      if(data[(t+i * *ntime)+n* *nrep]&&data[(t+j * *ntime)+n* *nrep]) n11s[h]++;
		      if(data[(t+i * *ntime)+n* *nrep]&&!data[(t+j * *ntime)+n* *nrep]) n10s[h]++;
		      if(!data[(t+i * *ntime)+n* *nrep]&&data[(t+j * *ntime)+n* *nrep]) n01s[h]++;
		      if(!data[(t+i * *ntime)+n* *nrep]&&!data[(t+j * *ntime)+n* *nrep]) n00s[h]++;}}}}}
	    else{// computes the spatial-temporal lorelogram:
	      if(mlags[i][j]<=*maxdist && mlagt[t][v]<=*maxtime){
		q=0;
		for(h=0;h<(*nbins-1);h++){
		  for(u=0;u<*nbint;u++){
		    if((bins[h]<=mlags[i][j]) && (mlags[i][j]<bins[h+1]) && (bint[u]==mlagt[t][v])){
		      for(n=0;n<*nrep;n++){
			if(data[(t+i * *ntime)+n* *nrep]&&data[(v+j * *ntime)+n* *nrep])  n11[q]++;
			if(data[(t+i * *ntime)+n* *nrep]&&!data[(v+j * *ntime)+n* *nrep])  n10[q]++;
			if(!data[(t+i * *ntime)+n* *nrep]&&data[(v+j * *ntime)+n* *nrep])  n01[q]++;
			if(!data[(t+i * *ntime)+n* *nrep]&&!data[(v+j * *ntime)+n* *nrep])  n00[q]++;
			}}
			q++;}}}}}}}
    }}
  // computing  space time log odds ratio in each spacetime bin
  q=0;
  for(h=0;h<(*nbins-1);h++){
   if(!n10s[h]||!n01s[h]||!n11s[h]||!n00s[h]||((n11s[h]*n00s[h])==(n10s[h]*n01s[h])))
     moms[h]=0;
   else moms[h]=log((n11s[h]*n00s[h])/(n01s[h]*n10s[h]));}
   for(u=0;u<(*nbint);u++){
     if(!n11t[u]||!n01t[u]||!n10t[u]||!n00t[u]||((n11t[u]*n00t[u])==(n10t[u]*n01t[u])))
       momt[u]=0;
     else momt[u]=log((n11t[u]*n00t[u])/(n01t[u]*n10t[u]));}
   for(q=0;q<((*nbins-1)*(*nbint));q++){
     if(!n11[q]||!n01[q]||!n10[q]||!n00[q]||((n11[q]*n00[q])==(n10[q]*n01[q])))
       momst[q]=0;
     else momst[q]=log((n11[q]*n00[q])/(n01[q]*n10[q]));}
  return;
}










// binned spatial variogram:
void Binned_Variogram(double *bins, double *data, int *lbins, double *moms, int *nbins)
{
  int h=0, i=0, j=0, n=0;
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
      if(mlags[i][j]<=*maxdist){
	for(h=0;h<(*nbins-1);h++)
	  if((bins[h]<=mlags[i][j]) && (mlags[i][j]<bins[h+1])){
	    for(n=0;n<*nrep;n++){
	      moms[h]+=0.5*pow(data[n+i * *nrep]-data[n+j * *nrep],2);
	      lbins[h]+=1;}}}}
  return;
}

// binned spatial variogram:
void Binned_Variogram_2(double *bins, double *data, int *lbins, double *moms, int *nbins)
{
  int h=0, i=0, j=0, n=0,p=0;
  double step=0.0;
  //Set the binnes step:
  step=(*maxdist-*minimdista)/(*nbins-1);
  bins[0]= *minimdista;
  //define bins:
  for(h=1;h<*nbins;h++)
    bins[h]=bins[h-1]+step;
  //Computes the binned moments:
  for(i=0;i<(*ncoord-1);i++){
    for(j=(i+1);j<*ncoord;j++){
      if(lags[p]<=*maxdist){
	for(h=0;h<(*nbins-1);h++)
	  if((bins[h]<=lags[p]) && (lags[p]<bins[h+1])){
	    for(n=0;n<*nrep;n++){
	      moms[h]+=0.5*pow(data[n+i * *nrep]-data[n+j * *nrep],2);
	      lbins[h]+=1;}}
	      p++;}}}
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
  step=*maxdist/(*nbins-1);
  bins[0]=0;
  for(h=1;h<*nbins;h++)
    bins[h]=bins[h-1]+step;
  //defines the temporal bins:
  bint[0]=0;
  for(h=1;h<*nbint;h++)
    bint[h]=bint[h-1]+minimtime[0];
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
	    else{// computes the spatial-temporal variogram:
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

// binned spatial-temporal variogram:
void Binned_Variogram_st2(double *bins, double *bint, double *data, int *lbins, int *lbinst,
			 int *lbint, double *moms,double *momst, double *momt, int *nbins, int *nbint)
{
int h=0, i=0, j=0, n=0, p=0;
  int q=0, t=0, u=0, v=0;
  double step=0.0;
  //defines the spatial bins:
  step=*maxdist/(*nbins-1);
  bins[0]=0;
  for(h=1;h<*nbins;h++)
    bins[h]=bins[h-1]+step;
  //defines the temporal bins:
  bint[0]=0;
  for(h=1;h<*nbint;h++)
    bint[h]=bint[h-1]+minimtime[0];
  //computes the empirical variogram:
  for(i=0;i<*ncoord;i++)
    for(t=0;t<*ntime;t++)
      for(j=i;j<*ncoord;j++){
	if(i==j)//computes the marignal temporal variogram:
	  for(v=(t+1);v<*ntime;v++){
	    if(lagt[p]<=*maxtime){
	      for(u=0;u<*nbint;u++){
		if(bint[u]==lagt[p]){
		  for(n=0;n<*nrep;n++){
		    momt[u]+=0.5*pow(data[(t+i * *ntime)+n* *nrep]-data[(v+i * *ntime)+n* *nrep], 2);
		    lbint[u]+=1;}}}p++;}}
	else{
	  for(v=0;v<*ntime;v++){
	    if(t==v){// computes the marginal spatial variogram:
	      if(lags[p]<=*maxdist){
		for(h=0;h<(*nbins-1);h++){
		  if((bins[h]<=lags[p]) && (lags[p]<bins[h+1])){
		    for(n=0;n<*nrep;n++){
		      moms[h]+=0.5*pow(data[(t+i * *ntime)+n* *nrep]-data[(t+j * *ntime)+n* *nrep], 2);
		      lbins[h]+=1;}}}p++;}}
	    else{// computes the spatial-temporal variogram:
	      if(lags[p]<=*maxdist && lagt[p]<=*maxtime){
		q=0;
		for(h=0;h<(*nbins-1);h++){
		  for(u=0;u<*nbint;u++){
		    if((bins[h]<=lags[p]) && (lags[p]<bins[h+1]) && (bint[u]==lagt[p])){
		      for(n=0;n<*nrep;n++){
			momst[q]+=0.5*pow(data[(t+i * *ntime)+n* *nrep]-data[(v+j * *ntime)+n* *nrep], 2);
			lbinst[q]+=1;}}
		    q++;}}
		    p++;}}}}}
  return;
}

// binned spatial madogram:
void Binned_Madogram(double *bins, double *data, int *lbins, double *moms, int *nbins)
{
  int h=0, i=0, j=0, n=0;
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
      if(mlags[i][j]<=*maxdist){
	for(h=0;h<(*nbins-1);h++)
	  if((bins[h]<=mlags[i][j]) && (mlags[i][j]<bins[h+1])){
	    for(n=0; n<*nrep;n++){
	      moms[h]+=0.5*fabs(data[n+i * *nrep]-data[n+j * *nrep]);
	      lbins[h]+=1;}}}
     }
  return;
}






// binned spatial-temporal variogram:
void Binned_Madogram_st(double *bins, double *bint, double *data, int *lbins, int *lbinst,
			 int *lbint, double *moms,double *momst, double *momt, int *nbins, int *nbint)
{
int h=0, i=0, j=0, n=0;
  int q=0, t=0, u=0, v=0;
  double step=0.0;
  //defines the spatial bins:
  step=*maxdist/(*nbins-1);
  bins[0]=0;
  for(h=1;h<*nbins;h++)
    bins[h]=bins[h-1]+step;
  //defines the temporal bins:
  bint[0]=0;
  for(h=1;h<*nbint;h++)
    bint[h]=bint[h-1]+minimtime[0];
  //computes the empirical variogram:
  for(i=0;i<*ncoord;i++)
    for(t=0;t<*ntime;t++)
      for(j=i;j<*ncoord;j++){
	if(i==j){//computes the marignal temporal variogram:
	  for(v=(t+1);v<*ntime;v++){
	    if(mlagt[t][v]<=*maxtime)
	      for(u=0;u<*nbint;u++)
		if(bint[u]==mlagt[t][v])
		  for(n=0;n<*nrep;n++){
		    momt[u]+=0.5*fabs(data[(t+i * *ntime)+n* *nrep]-data[(v+i * *ntime)+n* *nrep]);
		    lbint[u]+=1;}}}
	else{
	  for(v=0;v<*ntime;v++){
	    if(t==v){// computes the marginal spatial variogram:
	      if(mlags[i][j]<=*maxdist)
		for(h=0;h<(*nbins-1);h++)
		  if((bins[h]<=mlags[i][j]) && (mlags[i][j]<bins[h+1]))
		    for(n=0;n<*nrep;n++){
		      moms[h]+=0.5*fabs(data[(t+i * *ntime)+n* *nrep]-data[(t+j * *ntime)+n* *nrep]);
		      lbins[h]+=1;}}
	    else{// computes the spatial-temporal variogram:
	      if(mlags[i][j]<=*maxdist && mlagt[t][v]<=*maxtime){
		q=0;
		for(h=0;h<(*nbins-1);h++)
		  for(u=0;u<*nbint;u++){
		    if((bins[h]<=mlags[i][j]) && (mlags[i][j]<bins[h+1]) && (bint[u]==mlagt[t][v]))
		      for(n=0;n<*nrep;n++){
			momst[q]+=0.5*fabs(data[(t+i * *ntime)+n* *nrep]-data[(v+j * *ntime)+n* *nrep]);
			lbinst[q]+=1;}
		    q++;}}}}}}
  return;
}

// variogram cloud:
void Cloud_Variogram(double *bins, double *data, int *lbins, double *moms, int *nbins)
{
  int  h=0,i=0, j=0, n=0;
 //Computes the cloud moments:
  for(i=0;i<(*ncoord-1);i++)
    for(j=(i+1);j<*ncoord;j++){
      bins[h]=mlags[i][j];
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
       bins[h]=mlags[i][j];
      for(n=0;n<*nrep;n++)
	moms[h]+=0.5*fabs(data[n+i * *nrep]-data[n+j * *nrep]);
      lbins[h]=*nrep;
      h++;}
  return;
}
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
    for(h=0;h<(*nbins-1);h++){
      vario=moms[i]/lbins[i];// Computes the empirical variogram
      varhat=Variogram(cormod,0.5*(bins[h]+bins[h+1]),bint[u],nuis,par);
      *res=*res-pow(varhat-vario,2);// Computes the least squares
      i++;}
  return;
}
// Weighted least square method for Gaussian spatial-temporal random field:
void WLeastSquare_G(double *bins, double *bint, int *cormod, double *lbins, double *moms,
		       int *nbins, int *nbint, double *nuis, double *par, double *res)
{
  int h=0,i=0,u=0;
  double vario=0.0,varhat=0.0;
  //Checks the nuisance parameters (nugget, sill and corr):
  if(nuis[1]<0 || nuis[2]<=0 || CheckCor(cormod,par)==-2) {
    *res=LOW; return;}
  // Computes the least squares:
  for(u=0;u<*nbint;u++)
    for(h=0;h<(*nbins-1);h++){
      vario=moms[i]/lbins[i];// Computes the empirical variogram
      varhat=Variogram(cormod,0.5*(bins[h]+bins[h+1]),bint[u],nuis,par);
      if(vario) *res=*res-pow(varhat-vario,2)*(lbins[i]/pow(vario,2));// Computes the weighted least squares
      i++;}
  return;
}
// Least square method for max-stable extremal-Gaussian spatio-temporal random field:
void LeastSquare_MEG(double *bins, double *bint, int *cormod, double *lbins, double *moms,
		     int *nbins, int *nbint, double *nuis, double *par, double *res)
{
  int h=0,i=0,u=0;
  double vario=0.0,extcoeff=0.0,extcfhat=0.0;
  //Checks the nuisance and correlation parameters (sill, corr):
  if(nuis[0]<=0 || nuis[0]>1 || CheckCor(cormod,par)==-2) {
    *res=LOW; return;}
  // Computes the least squares:
  for(u=0;u<*nbint;u++)
    for(h=0;h<(*nbins-1);h++){
      vario=moms[i]/lbins[i];// Computes the madogram
      extcoeff=(1+2*vario)/(1-2*vario);// Computes the extremal coefficient
      // Computes the extremal coefficient for the Gaussian extremal model:
      extcfhat=1+sqrt(0.5*(1-nuis[0]*CorFct(cormod,0.5*(bins[h]+bins[h+1]),bint[u],par)));
      *res=*res-pow(extcoeff-extcfhat,2);// Computes the least squares
      i++;}
  return;
}
// Weighted least square method for max-stable extremal-Gaussian spatio-temporal random field:
void WLeastSquare_MEG(double *bins, double *bint, int *cormod, double *lbins, double *moms,
		      int *nbins, int *nbint, double *nuis, double *par, double *res)
{
  int h=0,i=0,u=0;
  double extcfhat=0.0,extcoeff=0.0,vario=0.0;
  //Checks the nuisance and correlation parameters (sill, corr):
  if(nuis[0]<=0 || nuis[0]>1 || CheckCor(cormod,par)==-2) {
    *res=LOW; return;}
 // Computes the least squares:
  for(u=0;u<*nbint;u++)
    for(h=0;h<(*nbins-1);h++){
      vario=moms[i]/lbins[i];// Computes the madogram
      extcoeff=(1+2*vario)/(1-2*vario);// Computes the extremal coefficient
      // Computes the extremal coefficient for the Gaussian extremal model:
      extcfhat=1+sqrt(0.5*(1-nuis[0]*CorFct(cormod,0.5*(bins[h]+bins[h+1]),bint[u],par)));
      *res=*res-pow(extcoeff/extcfhat-1,2)*lbins[h];
      i++;}
  return;
}
// Least square method for max-stable Brown-Resnick spatio-temporal random field:
void LeastSquare_MBR(double *bins, double *bint, int *cormod, double *lbins, double *moms,
		     int *nbins, int *nbint, double *nuis, double *par, double *res)
{
  int h=0,i=0,u=0;
  double vario=0.0,extcoeff=0.0,extcfhat=0.0;
  //Checks the nuisance and correlation parameters (sill, corr):
  if(nuis[0]<=0 || nuis[0]>1 || CheckCor(cormod,par)==-2) {
    *res=LOW; return;}
  // Computes the least squares:
  for(u=0;u<*nbint;u++)
    for(h=0;h<(*nbins-1);h++){
      vario=moms[i]/lbins[i];// Computes the madogram
      extcoeff=(1+2*vario)/(1-2*vario);// Computes the extremal coefficient
      // Computes the extremal coefficient for the Gaussian extremal model:
      extcfhat=2*pnorm(0.5*sqrt(VarioFct(cormod,0.5*(bins[h]+bins[h+1]),par,bint[u])),0,1,1,0);
      *res=*res-pow(extcoeff-extcfhat,2);// Computes the least squares
      i++;}
  return;
}
// Least square method for max-stable extremal-t spatio-temporal random field:
void LeastSquare_MET(double *bins, double *bint, int *cormod, double *lbins, double *moms,
		     int *nbins, int *nbint, double *nuis, double *par, double *res)
{
  int h=0,i=0,u=0;
  double df1=0.0,rho=0.0,vario=0.0,extcoeff=0.0,extcfhat=0.0;
  //Checks the nuisance and correlation parameters (sill, corr):
  if(nuis[0]<=0 || nuis[1]<=0 || nuis[1]>1 || CheckCor(cormod,par)==-2) {
    *res=LOW; return;}
  df1=nuis[0]+1;// Set the df of the t cdf
  // Computes the least squares:
  for(u=0;u<*nbint;u++)
    for(h=0;h<(*nbins-1);h++){
      vario=moms[i]/lbins[i];// Computes the madogram
      extcoeff=(1+2*vario)/(1-2*vario);// Computes the extremal coefficient
      rho=nuis[1]*CorFct(cormod,0.5*(bins[h]+bins[h+1]),bint[u],par);
      // Computes the extremal coefficient for the extremal-t model model:
      extcfhat=2*pt(sqrt(df1*(1-rho)/(1+rho)),df1,1,0);
      *res=*res-pow(extcoeff-extcfhat,2);// Computes the least squares
      i++;}
  return;
}

