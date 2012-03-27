#include "header.h"
#define REARTH 6378.388

void RangeDist(double *max, double *min)
{
  *max=*maximdista;
  *min=*minimdista;
  return;
}
// Computes the spatial distances on irregular and regural grid:
void Space_Dist(double *coordx,double *coordy,int grid,int *ia,int *idx,
		int *ja,double tapran,int tap,int type)
{
  int i=0,h=0,j=0,m=0,n=0;
  int icount=0,jcount=0;
  double dij=0.0,*tlags;

  if(tap){
    switch(type){
    case 0:// Euclidean distances:
      ia[0]=1;
      if(grid){// spatial grid
	for(i=0;i<*ncoordx;i++)
	  for(j=0;j<*ncoordy;j++){
	    jcount=0;
	    for(m=0;m<*ncoordx;m++)
	      for(n=0;n<*ncoordy;n++){
		dij=hypot(coordx[m]-coordx[i], coordy[n]-coordy[j]);
		*maximdista=fmax(*maximdista, dij);
		if(dij) *minimdista=fmin(*minimdista, dij);
		if(dij<=tapran){
		  lags[h]=dij;
		  ja[h]=jcount+1;
		  idx[h]=icount*(*ncoord)+jcount+1;
		  ia[icount+1]=ia[icount+1]+1;
		  h++;}
		jcount++;}
		icount++;}
	for(i=0;i<*ncoord;i++)
	  ia[i+1]=ia[i+1]+ia[i];}
      else{
	for(i=0;i<*ncoord;i++)
	  for(j=0;j<*ncoord;j++){
	    dij=hypot(coordx[i]-coordx[j],coordy[i]-coordy[j]);
	    *maximdista=fmax(*maximdista, dij);
	    if(dij) *minimdista=fmin(*minimdista,dij);
	    if(dij<= tapran){
	      lags[h]=dij;
	      ja[h]=j+1;
	      idx[h]=i*(*ncoord)+j+1;
	      ia[i+1]=ia[i+1]+1;
	      h++;}}
	for(i=0;i<*ncoord;i++)
	  ia[i+1]=ia[i+1]+ia[i];}
      break;
    case 1:// Geodesic distances:
        ia[0]=1;
      if(grid){// spatial grid
	for(i=0;i<*ncoordx;i++)
	  for(j=0;j<*ncoordy;j++){
	    jcount=0;
	    for(m=0;m<*ncoordx;m++)
	      for(n=0;n<*ncoordy;n++){
		dij=Dist_geodesic(coordx[m],coordy[n],coordx[i],coordy[j]);
		*maximdista=fmax(*maximdista,dij);
		if(dij) *minimdista=fmin(*minimdista,dij);
		if(dij<= tapran){
		  lags[h]=dij;
		  ja[h]=jcount+1;
		  idx[h]=icount*(*ncoord)+jcount+1;
		  ia[icount+1]=ia[icount+1]+1;
		  h++;}
		jcount++;}
		icount++;}
	for(i=0;i<*ncoord;i++)
	  ia[i+1]=ia[i+1]+ia[i];}
      else{
	for(i=0;i<*ncoord;i++)
	  for(j=0;j<*ncoord;j++){
	    dij=Dist_geodesic(coordx[i],coordy[i],coordx[j],coordy[j]);
	    *maximdista=fmax(*maximdista,dij);
	    if(dij) *minimdista=fmin(*minimdista,dij);
	    if(dij<= tapran){
	      lags[h]=dij;
	      ja[h]=j+1;
	      idx[h]=i*(*ncoord)+j+1;
	      ia[i+1]=ia[i+1]+1;
	      h++;}}
	for(i=0;i<*ncoord;i++)
	  ia[i+1]=ia[i+1]+ia[i];}
      break;}
    *npairs=h;
    tlags=(double *) malloc(*npairs*sizeof(double));
    for(i=0;i<*npairs;i++)
      tlags[i]=lags[i];
    lags=tlags;}
  else
    switch(type){
    case 0:// Euclidean distances:
      if(grid){// in case of a equispaced grid of coordinates:
	  for(i=0;i<*ncoordx;i++)
	    for(j=0;j<*ncoordy;j++){
		for(m=i;m<*ncoordx;m++){
		    if(m==i) n=(j+1);
		    else n=0;
		    for(n=n;n<*ncoordy;n++){
		      lags[h]=hypot(coordx[m]-coordx[i], coordy[n]-coordy[j]);
		      *maximdista=fmax(*maximdista, lags[h]);
		      *minimdista=fmin(*minimdista, lags[h]);
		      h++;}}}}
      else{// in case of an irregular grid of coordinates:
	  for(i=0;i<(*ncoord-1);i++)
	    for(j=(i+1);j<*ncoord;j++){
	      lags[h]=hypot(coordx[i]-coordx[j],coordy[i]-coordy[j]);
		*maximdista=fmax(*maximdista, lags[h]);
		*minimdista=fmin(*minimdista, lags[h]);
		h++;}}
      break;
    case 1:// Geodesic distances:
      if(grid){// in case of a equispaced grid of coordinates:
	for(i=0;i<*ncoordx;i++)
	  for(j=0;j<*ncoordy;j++){
	    for(m=i;m<*ncoordx;m++){
	      if(m==i) n=(j+1);
	      else n=0;
	      for(n=n;n<*ncoordy;n++){
		lags[h]=Dist_geodesic(coordx[m],coordy[n],coordx[i],coordy[j]);
		*maximdista=fmax(*maximdista, lags[h]);
		*minimdista=fmin(*minimdista, lags[h]);
		h++;}}}}
      else{// in case of a irregular grid of coordinates:
	  for(i=0;i<(*ncoord-1);i++)
	    for(j=(i+1);j<*ncoord;j++){
		lags[h]=Dist_geodesic(coordx[i],coordy[i],coordx[j],coordy[j]);
		*maximdista=fmax(*maximdista, lags[h]);
		*minimdista=fmin(*minimdista, lags[h]);
		h++;}}
      break;}
  return;
}
/*
void Space_Dist(double *coordx, double *coordy, int *grid, int *tap, double tapran, int *type)
{
  int i=0,icount=0,h=0,j=0,jcount=0,m=0,n=0;
  double dij=0.0;

  if(*tap)
    switch(*type){
    case 0:// Euclidean distances:
      ia[0]=1;
      if(*grid){
	for(i=0;i<*ncoordx;i++)
	  for(j=0;j<*ncoordy;j++){
	    jcount=0;
	    if(j>0) icount++;
	    for(m=0;m<*ncoordx;m++)
	      for(n=0;n<*ncoordy;n++){
		dij=hypot(coordx[m]-coordx[i], coordy[n]-coordy[j]);
		if(dij<= tapran){
		  lags[h]=dij;
		  ja[h]=jcount+1;
		  idx[h]=icount*(*ncoord)+jcount+1;
		  ia[icount+1]=ia[icount+1]+1;
		  *maximdista=fmax(*maximdista, lags[h]);
		  *minimdista=fmin(*minimdista, lags[h]);
		  h++;}
		jcount++;}}
      for(i=0;i<*ncoord;i++)
	  ia[i+1]=ia[i+1]+ia[i];}
      else{
	for(i=0;i<*ncoord;i++)
	  for(j=0;j<*ncoord;j++){
	    dij=hypot(coordx[i]-coordx[j],coordy[i]-coordy[j]);
	    if(dij<= tapran){
	      lags[h]=dij;
	      ja[h]=j+1;
	      idx[h]=i*(*ncoord)+j+1;
	      ia[i+1]=ia[i+1]+1;
	      *maximdista=fmax(*maximdista, lags[h]);
	      *minimdista=fmin(*minimdista, lags[h]);
	      h++;}}
	for (i=0;i<*ncoord;i++)
	  ia[i+1]=ia[i+1]+ia[i];}
      break;
    case 1:// Geodesic distances:
      ia[0]=1;
      for(i=0;i<*ncoord;i++)
	for(j=0;j<*ncoord;j++){
	  dij=Dist_geodesic(coordx[i],coordy[i],coordx[j],coordy[j]);
	  if(dij<= tapran){
	    lags[h]=dij;
	    ja[h]=j+1;
	    idx[h]=i*(*ncoord)+j+1;
	    ia[i+1]=ia[i+1]+1;
	    *maximdista=fmax(*maximdista, lags[h]);
	    *minimdista=fmin(*minimdista, lags[h]);
	    h++;}}
      for (i=0;i<*ncoord;i++)
	ia[i+1]=ia[i+1]+ia[i];
      break;}
  else
    switch(*type){
    case 0:// Euclidean distances:
      if(*grid){// in case of a equispaced grid of coordinates:
	  for(i=0;i<*ncoordx;i++)
	    for(j=0;j<*ncoordy;j++){
		for(m=i;m<*ncoordx;m++){
		    if(m==i) n=(j+1);
		    else n=0;
		    for(n=n;n<*ncoordy;n++){
		      lags[h]=hypot(coordx[m]-coordx[i], coordy[n]-coordy[j]);
		      *maximdista=fmax(*maximdista, lags[h]);
		      *minimdista=fmin(*minimdista, lags[h]);
		      h++;}}}}
      else{// in case of an irregular grid of coordinates:
	  for(i=0;i<(*ncoord-1);i++)
	    for(j=(i+1);j<*ncoord;j++){
	      lags[h]=hypot(coordx[i]-coordx[j],coordy[i]-coordy[j]);
		*maximdista=fmax(*maximdista, lags[h]);
		*minimdista=fmin(*minimdista, lags[h]);
		h++;}}
      break;
    case 1:// Geodesic distances:
      if(*grid){// in case of a equispaced grid of coordinates:
	for(i=0;i<*ncoordx;i++)
	  for(j=0;j<*ncoordy;j++){
	    for(m=i;m<*ncoordx;m++){
	      if(m==i) n=(j+1);
	      else n=0;
	      for(n=n;n<*ncoordy;n++){
		lags[h]=Dist_geodesic(coordx[m],coordy[n],coordx[i],coordy[j]);
		*maximdista=fmax(*maximdista, lags[h]);
		*minimdista=fmin(*minimdista, lags[h]);
		h++;}}}}
      else{// in case of a irregular grid of coordinates:
	  for(i=0;i<(*ncoord-1);i++)
	    for(j=(i+1);j<*ncoord;j++){
		lags[h]=Dist_geodesic(coordx[i],coordy[i],coordx[j],coordy[j]);
		*maximdista=fmax(*maximdista, lags[h]);
		*minimdista=fmin(*minimdista, lags[h]);
		h++;}}
      break;}
  return;
  }*/
// Computes the spatial-temporal distances on regular and irregular grid:
void SpaceTime_Dist(double *coordx, double *coordy, double *coordt, int *grid, int *type)
{
  int h=0, i=0, j=0, k=0, m=0, n=0, t=0, v=0;

  switch(*type)
    {
    case 0:// Euclidean distances:
      // Computes the spatial distances:
      if(*grid){// in case of a equispaced grid of coordinates:
	for(i=0;i<*ncoordx;i++)
	  for(j=0;j<*ncoordy;j++){
	    k=i* *ncoordx+j+1;
	    for(m=i;m<*ncoordx;m++){
	      if(m==i) n=(j+1);
	      else n=0;
	      for(n=n;n<*ncoordy;n++){
		mlags[h][k]=hypot(coordx[m]-coordx[i],coordy[n]-coordy[j]);
		*maximdista=fmax(*maximdista,mlags[h][k]);
		*minimdista=fmin(*minimdista,mlags[h][k]);
		k++;}} h++;}}
      else{// in case of an irregular grid of coordinates:
	for(i=0;i<*ncoord;i++)
	  for(j=(i+1);j<*ncoord;j++){
	    mlags[i][j]=hypot(coordx[i]-coordx[j],coordy[i]-coordy[j]);
	    *maximdista=fmax(*maximdista,mlags[i][j]);
	    *minimdista=fmin(*minimdista,mlags[i][j]);}}
      break;
    case 1:// Geodesic distances:
      // Computes the spatial distances:
      if(*grid){// in case of an equispaced grid of coordinates:
	for(i=0;i<*ncoordx;i++)
	  for(j=0;j<*ncoordy;j++){
	    k=i* *ncoordx+j+1;
	    for(m=i;m<*ncoordx;m++){
	      if(m==i) n=(j+1);
	      else n=0;
	      for(n=n;n<*ncoordy;n++){
		mlags[h][k]=Dist_geodesic(coordx[m],coordy[n],coordx[i],coordy[j]);
		*maximdista=fmax(*maximdista,mlags[h][k]);
		*minimdista=fmin(*minimdista,mlags[h][k]);
		k++;}} h++;}}
      else{// in case of an irregular grid of coordinates:
	for(i=0;i<*ncoord;i++)
	  for(j=(i+1);j<*ncoord;j++){
	    mlags[i][j]=Dist_geodesic(coordx[i],coordy[i],coordx[j],coordy[j]);
	    *maximdista=fmax(*maximdista,mlags[i][j]);
	    *minimdista=fmin(*minimdista,mlags[i][j]);}}
      break;}
  // Computes the temporal distances:
  for(t=0;t<*ntime;t++)
    for(v=(t+1);v<*ntime;v++){
      mlagt[t][v]=fabs(coordt[t]-coordt[v]);
      mlagt[v][t]=mlagt[t][v];
      *maximtime=fmax(*maximtime,mlagt[t][v]);
      *minimtime=fmin(*minimtime,mlagt[t][v]);}
  return;
}
// Computes the Geodesic distance between to coordinates:
double Dist_geodesic(double lonx, double latx, double lony, double laty)
/*
This function
compute geodetic inter-site distance between two
sites, given the latitude and longitude of both sites.
Input:

   lonx longitude site 1
   latx  latitude site 1
   lony longitude site 2
   laty  latitude site 2

Output:

   val the distance in km
*/
{
  double ax, bx, ay, by, val;
  val = 0.0;
  // this is a trick
  if (lonx == lony)
    {
      if(latx == laty)
	return val;
    }
  ax = (latx)*M_PI/180;
  bx = (lonx)*M_PI/180;
  ay = (laty)*M_PI/180;
  by = (lony)*M_PI/180;

  val = sin(ax) * sin(ay) + cos(ax) * cos(ay) * cos(fabs(bx - by));
  val = acos(val) *  REARTH;

  return val;
}

// Computes the spatial distances on irregular and regural grid:
void GeoDist(double *coordx, double *coordy, int *ncoord, double *res)
{
  int i=0, h=0, j=0;

  for(i=0;i<(*ncoord-1);i++)
    for(j=(i+1);j<*ncoord;j++){
      res[h]=Dist_geodesic(coordx[i],coordy[i],coordx[j],coordy[j]);
      h++;}

  return;
}

void ComputeMaxima(double *df, double *maxima, int *model,
		   int *nblock, int *nsite, double *sim)
{
  int i=0, k=0;
  double an=0.0, bn=0.0, chi2=1.0, ln=0, n=0.0;

  // Set the number of blocks:
  n=(double) *nblock;
  //Brown-Resnick (Husler-Reiss) field:
  if(*model==3){
    ln=log(n);
       bn=sqrt(2*ln)-(0.5*log(ln)+log(2*sqrt(M_PI)))/
	 sqrt(2*ln);
       an=1/bn;
       // First loop: number of blocks
       for(i=0;i<n;i++){
	 // Second loop: number of sites
	 for(k=0;k< *nsite;k++){
	   // Compute the componentwise maxima
	   maxima[k]=fmax(maxima[k],sim[k+i * *nsite]);
	   if(i==(n-1))
	     maxima[k]= exp((maxima[k]-bn)/an);}}}
  //Extremal-t field:
  if(*model==5){
    an=pow(n*gammafn((*df+1)/2)*pow(*df,*df/2-1)/gammafn(*df/2)/sqrt(M_PI),1/ *df);
    // First loop: number of blocks
      for(i=0;i<n;i++){
	chi2=sqrt(rchisq(*df) / *df);
	// Second loop: number of sites
	for(k=0;k< *nsite;k++){
	  // Compute the componentwise maxima
	  maxima[k]=fmax(maxima[k],sim[k+i * *nsite]/chi2);
	  if(i==(n-1))
	    maxima[k]=pow(maxima[k]/an, *df);}}}
  return;
}
// check if 'val1' is equal to 'val2'  when both are double
int is_equal(double val1, double val2)
{
  return fabs(val1-val2)<MAXERR;
}

double Maxima(double *x, int *size)
{
  double res=0.0;
  int i=0;

  res = x[0];

  for(i = 1; i < *size; i++)
    res = fmax(res, x[i]);

  return res;
}

double Minima(double *x, int *size)
{
  double res=0.0;
  int i=0;

  res = x[0];

  for(i = 1; i < *size; i++)
    res = fmin(res, x[i]);

  return res;
}

void Range(double *x, double *ran, int *size)
{
  int i=0;

  ran[0] = x[0];
  ran[1] = x[0];

  for(i = 1; i < *size; i++)
    {
      ran[0] = fmin(ran[0], x[i]);
      ran[1] = fmax(ran[1], x[i]);
    }

  return;
}
// define a sequence of points from x[0] to x[1] of length 'len'
void Seq(double *x, int len, double *res)
{
  double delta=0.0;
  int i=0;

  res[0] = x[0];
  delta = (x[1] - x[0]) / (len - 1);

  for(i = 1; i < len; i++)
    res[i] = res[i - 1] + delta;

  return;
}
// define a sequence of 'len' points of of 'delta' steps from the starting point x[0]
void SeqStep(double *x, int len, double step, double *res)
{
  int i=0;
  res[0]=x[0];
  for(i=0;i<len;i++) res[i+1]=res[i]+step;
  return;
}
// Determine (for the sub-sampling procedure) the sub-coordinates and
// the sub-data given spatial data, coordinates and an interval:
void SetSampling(double *coordx, double *coordy, double *data, int n, int *npts,
		 double *scoordx, double *scoordy, double *sdata, double xmax,
		 double xmin, double ymax, double ymin)
{
  int i=0, j=0;
  for(i=0;i<*ncoord;i++)
    if((xmin<coordx[i]||is_equal(xmin,coordx[i]))&&
       (xmax>coordx[i]||is_equal(xmax,coordx[i]))&&
       (ymin<coordy[i]||is_equal(ymin,coordy[i]))&&
       (ymax>coordy[i]||is_equal(ymax,coordy[i]))){
	scoordx[j]=coordx[i];
	scoordy[j]=coordy[i];
	sdata[j]=data[n* *nrep+i];
	j++;}
  *npts = j;
  return;
}
// Determine (for the sub-sampling procedure) the sub-coordinates and
// the sub-data given spatial data, coordinates and an interval:
void SetSampling_st(double *data,double *sdata,int *ncoord,int *ntime,
		    int wint,int k,int n,int *nrep)
{
  int i=0,j=0,p=0;
  for(i=0;i<(*ncoord);i++)
    for(j=(k+(*ntime*i))+n* *nrep;j<(k+wint+(*ntime*i))+n* *nrep;j++)
      {sdata[p]=data[j];p++;}
  return;
}
// Set the global variables for the spatial and spatial-temporal fitting:
void SetGlobalVar(double *coordx,double *coordy,double *coordt,int *grid,int *ia,
		  int *idx,int *ismal,int *ja,int *nsite,int *nsitex,int *nsitey,
		  int *npair,int *replic,int *spatim,double *srange,int *times,
		  double *trange,int *tap,int *type,int *weighted)
{
  //Spatial settings:
  maxdist=(double *) malloc(1*sizeof(double));//spatial threshould
  if(maxdist==NULL) {*ismal=0; return;}
  maxtime=(double *) malloc(1*sizeof(double));//temporal threshould
  if(maxtime==NULL) {*ismal=0; return;}
  maximdista=(double *) malloc(1*sizeof(double)); //maximum spatial distance
  if(maximdista==NULL) {*ismal=0; return;}
  *maximdista=0;
  minimdista=(double *) malloc(1*sizeof(double)); //minimum spatial distance
  if(minimdista==NULL) {*ismal=0; return;}
  *minimdista=1.0e15;
  ncoord=(int *) malloc(1 * sizeof(int));//number of total spatial coordinates
  if(ncoord==NULL) {*ismal=0; return;}
  *ncoord=*nsite;
  ncoordx=(int *) malloc(1 * sizeof(int));//number of the first spatial coordinates
  if(ncoordx==NULL) {*ismal=0; return;}
  *ncoordx=*nsitex;
  ncoordy=(int *) malloc(1 * sizeof(int));//number of the second spatial coordinates
  if(ncoordy==NULL) {*ismal=0; return;}
  *ncoordy=*nsitey;
  npairs=(int *) malloc(1 * sizeof(int));//number of spatial pairs
  if(npairs==NULL) {*ismal=0; return;}
  //Temporal settings:
  maximtime=(double *) malloc(1*sizeof(double)); //maximum temporal distance
  if(maximtime==NULL) {*ismal=0; return;}
  minimtime=(double *) malloc(1*sizeof(double)); //minimum temporal distance
  if(minimtime==NULL) {*ismal=0; return;}
  ntime=(int *) malloc(1*sizeof(int));//number of times
  if(ntime==NULL) {*ismal=0; return;}
  *ntime=*times;
  if (!*tap) npairs[0]=(ncoord[0]*ntime[0])*((ncoord[0]*ntime[0])-1)/2;
  else npairs[0]=ncoord[0]*ntime[0]*ncoord[0]*ntime[0];
  //Random field replications:
  nrep=(int *) malloc(1 * sizeof(int));//number of iid replicates of the random field
  if(nrep==NULL) {*ismal=0; return;}
  *nrep=*replic;
  // Computes the spatial or spatial-temporal distances
  // and the minima and maxima of these distances:
  if(*spatim){// Computes the vectors of spatia-temporal distances:
      int i=0;
      *maximtime=0;// set the initial maximum time
      *minimtime=1.0e15;// set the initial minimum time
      // allocates the vector of temporal distances:
      mlagt=malloc(*ntime*sizeof(double *));
      if(mlagt==NULL) {*ismal=0; return;}
      for(i=0;i<*ntime;i++){
	  mlagt[i]=malloc(*ntime*sizeof(double));
	  if(mlagt[i]==NULL) {*ismal=0; return;}
	  mlagt[i][i]=0;}
      // allocates the matrix of spatial distances:
      mlags=malloc(*ncoord*sizeof(double *));
      if(mlags==NULL) {*ismal=0; return;}
      for(i=0;i<*ncoord;i++){
	  mlags[i]=malloc(*ncoord*sizeof(double));
	  if(mlags[i]==NULL) {*ismal=0; return;}
	  mlags[i][i]=0;}
      // set the vector of temporal distances and
      // the matrix of spatial distances:
      SpaceTime_Dist(coordx, coordy, coordt, grid, type);
      // Set the range of the temporal intervals:
      trange[0]=*minimtime;
      *maxtime=*maximtime;
      if(trange[1]!=0) {if(trange[1]>*minimtime && trange[1]<*maximtime)
		      *maxtime=trange[1];}
      else trange[1]=*maxtime;}
  else{// Computes the vector of spatial distances:
      // allocates the vectors of spatial distances
    lags=(double *) malloc(*npairs*sizeof(double));
    if(lags==NULL){*ismal=0; return;}
    Space_Dist(coordx,coordy,grid[0],ia,idx,ja,srange[1],tap[0],type[0]);}
  // Set the range of the spatial distances:
  srange[0]=*minimdista;
  *maxdist=*maximdista;
  if(srange[1]!=0) {if(srange[1]>*minimdista && srange[1]<*maximdista)
		  *maxdist=srange[1];}
  else srange[1]=*maxdist;
  npair[0]=npairs[0];
  return;
}
/*
void indx(double *ind, int *n)
{
  int i=0, j=0, h=0;

  for(i=0;i<*n;i++)
    for(j=0;j<i;j++)
      {
	if(j==0) ind[h]=i-1;
	else ind[h]=ind[h-1]+(*n-j-1);
	h++;
      }
  return;
}
*/
void DeleteGlobalVar()
{
  // Delete all the global variables:
  free(maxdist); maxdist=NULL;
  free(maxtime); maxtime=NULL;
  free(lags); lags=NULL;
  //free(lagt); lagt=NULL;
  free(maximdista); maximdista=NULL;
  free(maximtime); maximtime=NULL;
  free(minimdista); minimdista=NULL;
  free(minimtime); minimtime=NULL;
  free(mlags); mlags=NULL;
  free(mlagt); mlagt=NULL;
  free(ncoord); ncoord=NULL;
  free(ncoordx); ncoord=NULL;
  free(ncoordy); ncoord=NULL;
  free(npairs); npairs=NULL;
  free(npairt); npairt=NULL;
  free(nrep); nrep=NULL;
  free(ntime); ntime=NULL;
  return;
}
