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
		int *ismal,int *ja,double thres,int type)

{
  int i=0,h=0,j=0,k=0, m=0, n=0, z=0;
  double dij=0.0;
  if(*istap){   // tapering case
    switch(type){
    case 0:// Euclidean distances:
      ia[0]=1;
  /******************************************************************************/
      if(grid){// spatial grid
     int icount=0,jcount=0;
	for(i=0;i<*ncoordx;i++){
	  for(j=0;j<*ncoordy;j++){
	    jcount=0;
	    for(m=0;m<*ncoordx;m++)
	      for(n=0;n<*ncoordy;n++){
		dij=hypot(coordx[m]-coordx[i], coordy[n]-coordy[j]);
		*maximdista=fmax(*maximdista, dij);
		if(dij) *minimdista=fmin(*minimdista, dij);
		if(dij<=thres){
		  lags[h]=dij;
		  ja[h]=jcount+1;
		  idx[h]=icount*(*ncoord)+jcount+1;
		  ia[icount+1]=ia[icount+1]+1;
		  h++;}
		jcount++;}
		icount++;}}
	for(i=0;i<*ncoord;i++)
	  ia[i+1]=ia[i+1]+ia[i];}
      else{   //no spatial grid
	for(i=0;i<*ncoord;i++)
	  for(j=0;j<*ncoord;j++){
	    dij=hypot(coordx[i]-coordx[j],coordy[i]-coordy[j]);
	    *maximdista=fmax(*maximdista, dij);
	    if(dij) *minimdista=fmin(*minimdista,dij);
	    if(dij<= thres){
	      lags[h]=dij;
	      ja[h]=j+1;
	      idx[h]=i*(*ncoord)+j+1;
	      ia[i+1]=ia[i+1]+1;
	      h++;}}
	for(i=0;i<*ncoord;i++)
	  ia[i+1]=ia[i+1]+ia[i];}
      break;
/******************************************************************************/
/******************************************************************************/
    case 1:// Chordal distances:
        ia[0]=1;
      if(grid){// spatial grid
           int icount=0,jcount=0;
	for(i=0;i<*ncoordx;i++)
	  for(j=0;j<*ncoordy;j++){
	    jcount=0;
	    for(m=0;m<*ncoordx;m++)
	      for(n=0;n<*ncoordy;n++){
		dij=Dist_chordal(coordx[m],coordy[n],coordx[i],coordy[j]);
		*maximdista=fmax(*maximdista,dij);
		if(dij) *minimdista=fmin(*minimdista,dij);
		if(dij<= thres){
		  lags[h]=dij;
		  ja[h]=jcount+1;
		  idx[h]=icount*(*ncoord)+jcount+1;
		  ia[icount+1]=ia[icount+1]+1;
		  h++;}
		jcount++;}
		icount++;}
	for(i=0;i<*ncoord;i++)
	  ia[i+1]=ia[i+1]+ia[i];}
      else{ //no spatial grid
	for(i=0;i<*ncoord;i++)
	  for(j=0;j<*ncoord;j++){
	    dij=Dist_chordal(coordx[i],coordy[i],coordx[j],coordy[j]);
	    *maximdista=fmax(*maximdista,dij);
	    if(dij) *minimdista=fmin(*minimdista,dij);
	    if(dij<= thres){
	      lags[h]=dij;
	      ja[h]=j+1;
	      idx[h]=i*(*ncoord)+j+1;
	      ia[i+1]=ia[i+1]+1;
	      h++;}}
	for(i=0;i<*ncoord;i++)
	  ia[i+1]=ia[i+1]+ia[i];}
      break;
/******************************************************************************/
/******************************************************************************/
    case 2:// Geodesic distances:
        ia[0]=1;
      if(grid){// spatial grid
           int icount=0,jcount=0;
	for(i=0;i<*ncoordx;i++)
	  for(j=0;j<*ncoordy;j++){
	    jcount=0;
	    for(m=0;m<*ncoordx;m++)
	      for(n=0;n<*ncoordy;n++){
           dij=Dist_geodesic(coordx[m],coordy[n],coordx[i],coordy[j]);
		*maximdista=fmax(*maximdista,dij);
		if(dij) *minimdista=fmin(*minimdista,dij);
		if(dij<= thres){
		 lags[h]=dij;
		  ja[h]=jcount+1;
		  idx[h]=icount*(*ncoord)+jcount+1;
		  ia[icount+1]=ia[icount+1]+1;
		  h++;}
		jcount++;}
		icount++;}
	for(i=0;i<*ncoord;i++)
	  ia[i+1]=ia[i+1]+ia[i];}
      else{ // no spatial grid
	for(i=0;i<*ncoord;i++)
	  for(j=0;j<*ncoord;j++){
          dij=Dist_geodesic(coordx[i],coordy[i],coordx[j],coordy[j]);
	    *maximdista=fmax(*maximdista,dij);
	    if(dij) *minimdista=fmin(*minimdista,dij);
	    if(dij<= thres){
	      lags[h]=dij;
	      ja[h]=j+1;
	      idx[h]=i*(*ncoord)+j+1;
	      ia[i+1]=ia[i+1]+1;
	      h++;}}
	for(i=0;i<*ncoord;i++)
	  ia[i+1]=ia[i+1]+ia[i];}
      break;
      }
/******************************************************************************/
/******************************************************************************/
    double *tlags;   //saving  distance for tapering
    npairs[0]=h;
    tlags=(double *)  Calloc(h,double);
    for(i=0;i<h;i++) tlags[i]=lags[i];
    Free(lags);
    lags=(double *)  Calloc(h,double);
    //lags=tlags;    // lags is the vector of spatial distances
    for(i=0;i<h;i++) lags[i]=tlags[i];
    Free(tlags);
    }  //end tapering

  else{  //no tapering: classical case
  switch(type){
/******************************************************************************/
/******************************************************************************/
      case 0:// Euclidean distances:
	if(grid){// in case of a equispaced grid of coordinates:
	  for(i=0;i<*ncoordx;i++){
	    for(j=0;j<*ncoordy;j++){
	         k=i* *ncoordx+j+1;

	      for(m=i;m<*ncoordx;m++){
		if(m==i) n=(j+1);
		else n=0;
		for(z=n;z<*ncoordy;z++){
		  mlags[h][k]=hypot(coordx[m]-coordx[i], coordy[z]-coordy[j]);
		  //Rprintf("%f %d %d\n",mlags[h][k],h,k);
		  *maximdista=fmax(*maximdista, mlags[h][k]);
		  *minimdista=fmin(*minimdista, mlags[h][k]);
		  k++;}};h++;}}
		  }
	else{// in case of an irregular grid of coordinates:
	  for(i=0;i<(*ncoord-1);i++){
	    for(j=(i+1);j<*ncoord;j++){
	      mlags[i][j]=hypot(coordx[i]-coordx[j],coordy[i]-coordy[j]);
	      *maximdista=fmax(*maximdista,mlags[i][j]);
	      *minimdista=fmin(*minimdista,mlags[i][j]);
	      h++;}}}
	break;
/******************************************************************************/
/******************************************************************************/
      case 1:// Chordal distances:
	if(grid){// in case of a equispaced grid of coordinates:
	  for(i=0;i<*ncoordx;i++){
	    for(j=0;j<*ncoordy;j++){
	         k=i* *ncoordx+j+1;
	      for(m=i;m<*ncoordx;m++){
		if(m==i) n=(j+1);
		else n=0;
		for(z=n;z<*ncoordy;z++){
            mlags[h][k]=Dist_chordal(coordx[m],coordy[z],coordx[i],coordy[j]);
		*maximdista=fmax(*maximdista, mlags[h][k]);
		  *minimdista=fmin(*minimdista, mlags[h][k]);
		  k++;}};h++;}}
		  }
	else{// in case of a irregular grid of coordinates:
	  for(i=0;i<(*ncoord-1);i++)
	    for(j=(i+1);j<*ncoord;j++){
          mlags[i][j]=Dist_chordal(coordx[i],coordy[i],coordx[j],coordy[j]);
	      *maximdista=fmax(*maximdista, mlags[i][j]);
	      *minimdista=fmin(*minimdista, mlags[i][j]);
	      h++;}}
	break;
/******************************************************************************/
/******************************************************************************/
	case 2:// Geodesic distances:
	if(grid){// in case of a equispaced grid of coordinates:
	    for(i=0;i<*ncoordx;i++){
	    for(j=0;j<*ncoordy;j++){
	         k=i* *ncoordx+j+1;
	      for(m=i;m<*ncoordx;m++){
		if(m==i) n=(j+1);
		else n=0;
		for(z=n;z<*ncoordy;z++){
            mlags[h][k]=Dist_geodesic(coordx[m],coordy[z],coordx[i],coordy[j]);
             //Rprintf("%f %d %d\n",mlags[h][k],h,k);
		  *maximdista=fmax(*maximdista, mlags[h][k]);
		  *minimdista=fmin(*minimdista, mlags[h][k]);
		  k++;}};h++;}}
		  }
	else{// in case of a irregular grid of coordinates:
	  for(i=0;i<(*ncoord-1);i++)
	    for(j=(i+1);j<*ncoord;j++){
           mlags[i][j]=Dist_geodesic(coordx[i],coordy[i],coordx[j],coordy[j]);
           //Rprintf("%f\n",mlags[i][j]);
	      *maximdista=fmax(*maximdista, mlags[i][j]);
	      *minimdista=fmin(*minimdista, mlags[i][j]);
	      h++;}}
	break;
	/******************************************************************************/
/******************************************************************************/
	}}

	return;
}



// Computes the spatial-temporal distances on regular and irregular grid:
void SpaceTime_Dist(double *coordx,double *coordy,double *coordt,int *grid,int *ia,int *idx,int *ismal,int *ja,
                    int *tapmodel,double thres,double thret,int type)
{
  int i=0,cc=0,j=0,h=0,k=0,m=0,n=0,t=0,v=0, z=0;

  if (*istap) {       // tapering case
  double *thre,*c_supp;
  c_supp=(double *) R_alloc(2, sizeof(double));
  thre=(double *)   R_alloc(2, sizeof(double));
  thre[0]=thres;thre[1]=thret;
  double dij=0.0,dtv=0.0;
  int count=0;
  if(*grid){   // it doesn't work
    int icount=0,jcount=0;
    switch(type){
    case 0:  // euclidean distance   no funciona!!!

          ia[0] = 1;
          for(i=0;i<*ncoord;i++){
          for(t=0;t<*ntime;t++){
          cc=0;
          for(j=0;j<*ncoord;j++){
          jcount=0;
            for(v=0;v<*ntime;v++){
               dtv=fabs(coordt[t]-coordt[v]);
               *maximtime=fmax(*maximtime, dtv);
               if(dtv) *minimtime=fmin(*minimtime, dtv);

	      for(m=0;m<*ncoordx;m++){
	      for(n=0;n<*ncoordy;n++){
		dij=hypot(coordx[m]-coordx[i], coordy[n]-coordy[j]);
		*maximdista=fmax(*maximdista, dij);
		if(dij) *minimdista=fmin(*minimdista, dij);
                Comp_supp(c_supp,tapmodel, dij, dtv,thre);
              if(dij<=c_supp[0]&&dtv<=c_supp[1]){
                               lags[count]=dij;
                               lagt[count]=dtv;
                               idx[count] =(icount * (*ntime) * (*ntime) * *ncoord) +  (t*  *ncoord *  *ntime) +  (1+v+ *ntime * jcount);
                               ja[count]=1+v+(*ntime) * jcount;
                               cc=cc+1;
                               count = count +1 ;}
                }
                jcount++;}}
                icount++;}
                ia[k+1]=ia[k]+cc;
                k=k+1;}}
          break;
    case 1:  // chordal distance
    break;
    case 2:  // geod distance
    break;
    }}
     else{
  switch(type){
  case 0:  // euclidean distance
   ia[0] = 1;
          for(i=0;i<*ncoord;i++){
          for(t=0;t<*ntime;t++){
          cc=0;
          for(j=0;j<*ncoord;j++){
          dij=hypot(coordx[i]-coordx[j],coordy[i]-coordy[j]);
          *maximdista=fmax(*maximdista, dij);
          if(dij) *minimdista=fmin(*minimdista, dij);
          for(v=0;v<*ntime;v++){
               dtv=fabs(coordt[t]-coordt[v]);
               *maximtime=fmax(*maximtime, dtv);
               if(dtv) *minimtime=fmin(*minimtime, dtv);
                Comp_supp(c_supp,tapmodel, dij, dtv,thre);
              if(dij<=c_supp[0]&&dtv<=c_supp[1]){
                               lags[count]=dij;
                               lagt[count]=dtv;
                               idx[count] =(i * (*ntime) * (*ntime) * *ncoord) +  (t*  *ncoord *  *ntime) +  (1+v+ *ntime * j);
                               ja[count]=1+v+(*ntime) * j;
                               cc=cc+1;
                               count = count +1 ;
                }}}
                ia[k+1]=ia[k]+cc;
                k=k+1;}}
          break;
  case 1:
         ia[0] = 1;
          for(i=0;i<*ncoord;i++){
          for(t=0;t<*ntime;t++){
          cc=0;
          for(j=0;j<*ncoord;j++){
          dij=Dist_chordal(coordx[i],coordy[i],coordx[j],coordy[j]);
          *maximdista=fmax(*maximdista, dij);
          if(dij) *minimdista=fmin(*minimdista, dij);
          for(v=0;v<*ntime;v++){
               dtv=fabs(coordt[t]-coordt[v]);
               *maximtime=fmax(*maximtime, dtv);
               if(dtv) *minimtime=fmin(*minimtime, dtv);
                Comp_supp(c_supp,tapmodel, dij, dtv,thre);
              if(dij<=c_supp[0]&&dtv<=c_supp[1]){
                               lags[count]=dij;
                               lagt[count]=dtv;
                               idx[count] =(i * (*ntime) * (*ntime) * *ncoord) +  (t*  *ncoord *  *ntime) +  (1+v+ *ntime * j);
                               ja[count]=1+v+(*ntime) * j;
                               cc=cc+1;
                               count = count +1 ;}
                }}
                ia[k+1]=ia[k]+cc;
                k=k+1;}}
          break;
  case 2:
   ia[0] = 1;
          for(i=0;i<*ncoord;i++){
          for(t=0;t<*ntime;t++){
          cc=0;
          for(j=0;j<*ncoord;j++){
          dij=Dist_geodesic(coordx[i],coordy[i],coordx[j],coordy[j]);
          *maximdista=fmax(*maximdista, dij);
          if(dij) *minimdista=fmin(*minimdista, dij);
          for(v=0;v<*ntime;v++){
               dtv=fabs(coordt[t]-coordt[v]);
               *maximtime=fmax(*maximtime, dtv);
               if(dtv) *minimtime=fmin(*minimtime, dtv);
                Comp_supp(c_supp,tapmodel, dij, dtv,thre);
              if(dij<=c_supp[0]&&dtv<=c_supp[1]){
                               lags[count]=dij;
                               lagt[count]=dtv;
                               idx[count] =(i * (*ntime) * (*ntime) * *ncoord) +  (t*  *ncoord *  *ntime) +  (1+v+ *ntime * j);
                               ja[count]=1+v+(*ntime) * j;
                               cc=cc+1;
                               count = count +1 ;}
                }}
                ia[k+1]=ia[k]+cc;
                k=k+1;}}
          break;
         }}
  npairs[0]=count;
  double *tlags, *tlagt;
  tlags=(double *) Calloc(count,double);
  tlagt=(double *) Calloc(count,double);
  for(i=0;i<count;i++) { tlags[i]=lags[i];tlagt[i]=lagt[i];}
  Free(lags);Free(lagt);
  lags=(double *)  Calloc(count,double);
  lagt=(double *)  Calloc(count,double);
  for(i=0;i<count;i++) {lags[i]=tlags[i];lagt[i]=tlagt[i];}
  Free(tlags);Free(tlagt);
  }

  else {   // no tapering


  switch(type){
	  case 0:// Euclidean distances:
      // Computes the spatial distances:
      if(*grid){// in case of a equispaced grid of coordinates:
        for(i=0;i<*ncoordx;i++){
	    for(j=0;j<*ncoordy;j++){
	         k=i* *ncoordx+j+1;
	      for(m=i;m<*ncoordx;m++){
		if(m==i) n=(j+1);
		else n=0;
		for(z=n;z<*ncoordy;z++){
		  mlags[h][k]=hypot(coordx[m]-coordx[i], coordy[z]-coordy[j]);
		  //Rprintf("fgdfg  %d %d %f \n",h,k,mlags[h][k]);
		  *maximdista=fmax(*maximdista, mlags[h][k]);
		  *minimdista=fmin(*minimdista, mlags[h][k]);
		  k++;}};h++;}}
		  }
      else{// in case of an irregular grid of coordinates:
	for(i=0;i<*ncoord;i++)
	  for(j=(i+1);j<*ncoord;j++){
	    mlags[i][j]=hypot(coordx[i]-coordx[j],coordy[i]-coordy[j]);
	     //Rprintf("fgdfg  %d %d %f \n",i,j,mlags[i][j]);
	    *maximdista=fmax(*maximdista,mlags[i][j]);
	    *minimdista=fmin(*minimdista,mlags[i][j]);}}
      break;

      case 1:// chordal distances:
      // Computes the spatial distances:
      if(*grid){// in case of a equispaced grid of coordinates:
	  for(i=0;i<*ncoordx;i++){
	    for(j=0;j<*ncoordy;j++){
	         k=i* *ncoordx+j+1;
	      for(m=i;m<*ncoordx;m++){
		if(m==i) n=(j+1);
		else n=0;
		for(z=n;z<*ncoordy;z++){
            mlags[h][k]=Dist_chordal(coordx[m],coordy[z],coordx[i],coordy[j]);
		*maximdista=fmax(*maximdista, mlags[h][k]);
		  *minimdista=fmin(*minimdista, mlags[h][k]);
		  k++;}};h++;}}
		}
      else{// in case of an irregular grid of coordinates:
	for(i=0;i<*ncoord;i++)
	  for(j=(i+1);j<*ncoord;j++){
	    mlags[i][j]=Dist_chordal(coordx[i],coordy[i],coordx[j],coordy[j]);
          Dist_chordal(coordx[i],coordy[i],coordx[j],coordy[j]);
	    *maximdista=fmax(*maximdista,mlags[i][j]);
	    *minimdista=fmin(*minimdista,mlags[i][j]);}}
      break;


    case 2:// Geodesic distances:
      // Computes the spatial distances:
      if(*grid){// in case of an equispaced grid of coordinates:
        for(i=0;i<*ncoordx;i++){
	    for(j=0;j<*ncoordy;j++){
	         k=i* *ncoordx+j+1;
	      for(m=i;m<*ncoordx;m++){
		if(m==i) n=(j+1);
		else n=0;
		for(z=n;z<*ncoordy;n++){
            mlags[h][k]=Dist_geodesic(coordx[m],coordy[z],coordx[i],coordy[j]);
		  *maximdista=fmax(*maximdista, mlags[h][k]);
		  *minimdista=fmin(*minimdista, mlags[h][k]);
		  k++;}};h++;}}
		}
      else{// in case of an irregular grid of coordinates:
	for(i=0;i<*ncoord-1;i++)
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
  }

  return;
}
double Dist_chordal(double loni, double lati, double lonj, double latj)
 {
   double ai, bi, aj, bj, val=0.0;
   if (loni == lonj && lati == latj) return val;
   ai = (lati)*M_PI/180;
   bi = (loni)*M_PI/180;
   aj = (latj)*M_PI/180;
   bj = (lonj)*M_PI/180;
 val=REARTH  *sqrt(R_pow(cos(ai) * cos(bi)-cos(aj)  *cos(bj) ,2) +
                    R_pow(cos(ai) * sin(bi)-cos(aj) * sin(bj) ,2)+
                         R_pow(sin(ai)-sin(aj) ,2));
 return val;
 }

// Computes the Geodesic distance between to coordinates:
double Dist_geodesic(double loni, double lati, double lonj, double latj)
{
  double ai, bi, aj, bj, val;
  val = 0.0;
 if (loni == lonj && lati == latj) return val;
  ai = (lati)*M_PI/180;
  bi = (loni)*M_PI/180;
  aj = (latj)*M_PI/180;
  bj = (lonj)*M_PI/180;
  val = sin(ai) * sin(aj) + cos(ai) * cos(aj) * cos(bi - bj);
  val = acos(val) *  REARTH;
  return val;
}
// Computes the spatial distances on irregular and regural grid:
void GeoDist(double *coordx, double *coordy, int *ncoord, double *res,int *type_dist)
{
  int i=0, h=0, j=0;

  for(i=0;i<(*ncoord-1);i++)
    for(j=(i+1);j<*ncoord;j++){
      if(*type_dist==1) res[h]=Dist_chordal(coordx[i],coordy[i],coordx[j],coordy[j]);
      if(*type_dist==2) res[h]=Dist_geodesic(coordx[i],coordy[i],coordx[j],coordy[j]);
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
// define a sequence of 'len' points of 'delta' steps from the starting point x[0]
void SeqStep(double *x, int len, double step, double *res)
{
  int i=0;
  res[0]=x[0];
  for(i=1;i<len;i++) res[i]=res[i-1]+step;
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
		  int *npair,int *replic,double *srange, double *sep, int *times,double *trange,
		  int *tap,int *tapmodel,int *type,int *weighted)
{
  //Spatial settings: //Spatial settings:
  maxdist=(double *) Calloc(1,double);//spatial threshold
  if(maxdist==NULL) {*ismal=0; return;}
  maximdista=(double *) Calloc(1,double);//spatial threshold
  if(maximdista==NULL) {*ismal=0; return;}
  *maximdista=0;
  minimdista=(double *) Calloc(1,double); //minimum spatial distance
  if(minimdista==NULL) {*ismal=0; return;}
  *minimdista=-LOW;
  ncoord=(int *) Calloc(1,int);//number of total spatial coordinates
  if(ncoord==NULL) {*ismal=0; return;}
  *ncoord=*nsite;
  ncoordx=(int *) Calloc(1,int);//number of the first spatial coordinates
  if(ncoordx==NULL) {*ismal=0; return;}
  *ncoordx=*nsitex;
  ncoordy=(int *) Calloc(1,int);//number of the second spatial coordinates
  if(ncoordy==NULL) {*ismal=0; return;}
  *ncoordy=*nsitey;
  npairs=(int *) Calloc(1,int);//effective number of pairs
  if(npairs==NULL) {*ismal=0; return;}
  isst=(int *) Calloc(1,int);//is a spatio-temporal random field?
  if(isst==NULL) {*ismal=0; return;}
  if(*times>1) *isst=1; else *isst=0;//set the spatio-temporal flag
  istap=(int *) Calloc(1,int);//is tapering?
  if(istap==NULL) {*ismal=0; return;}
  *istap=*tap;//set the tapering flag
  //set the total number of pairs
  //Random field replications:
  nrep=(int *) Calloc(1,int);//number of iid replicates of the random field
  if(nrep==NULL) {*ismal=0; return;}
  *nrep=*replic;
  tapsep=(double *) Calloc(1,double);//temporal threshold
  if(tapsep==NULL) {*ismal=0; return;}
  *tapsep=sep[0];

  // Computes the spatial or spatial-temporal distances
  // and the minima and maxima of these distances:
   if(!*isst) {// spatial case
        if(*istap) {// tapering case
           npairs[0]=ncoord[0]*ncoord[0];
           lags=(double *) Calloc(*npairs,double);
           if(lags==NULL){*ismal=0; return;}}
      else {
           int i=0;
           npairs[0]=ncoord[0]*(ncoord[0]-1)*0.5;
           mlags= (double **) Calloc(ncoord[0],double *);
           if(mlags==NULL) {*ismal=0; return;}
           for(i=0;i<ncoord[0];i++){
              mlags[i]=(double *) Calloc(ncoord[0],double);
              if(mlags[i]==NULL) {*ismal=0; return;}}  // mlags[i][j] matrix of spatial distances
           }
      Space_Dist(coordx,coordy,grid[0],ia,idx,ismal,ja,srange[1],type[0]);
      if(!ismal[0]) return;
      }
     else { //spatio temporal case
    maxtime=(double *) Calloc(1,double); //temporal threshold
    if(maxtime==NULL) {*ismal=0; return;}
    maximtime=(double *) Calloc(1,double);//maximum temporal distance
    if(maximtime==NULL) {*ismal=0; return;}
    minimtime=(double *) Calloc(1,double);//minimum temporal distance
    if(minimtime==NULL) {*ismal=0; return;}
    ntime=(int *) Calloc(1,int);//number of times
    if(ntime==NULL) {*ismal=0; return;}
    *ntime=*times;
    *maximtime=0;// set the initial maximum time
    *minimtime=-LOW;// set the initial minimum time
    if(*istap) {// tapering case
           npairs[0]=pow(ncoord[0]*ntime[0],2);
           lags=(double *) Calloc(*npairs,double);
           if(lags==NULL){*ismal=0; return;}
           lagt=(double *) Calloc(*npairs,double);
           if(lagt==NULL){*ismal=0; return;}
           }
      else {
             int i=0;
             npairs[0]=(ncoord[0]*ntime[0])*(ncoord[0]*ntime[0]-1)*0.5;
          // allocates the matrix of temporal distances:
          mlagt= (double **) Calloc(ntime[0],double *);
          if(mlagt==NULL) {*ismal=0; return;}
          for(i=0;i<*ntime;i++){
              mlagt[i]=(double *) Calloc(ntime[0],double);
              if(mlagt[i]==NULL) {*ismal=0; return;}}
             // allocates the matrix of spatial distances:
          mlags= (double **) Calloc(ncoord[0],double *);
          if(mlags==NULL) {*ismal=0; return;}
          for(i=0;i<ncoord[0];i++){
              //mlags[i]=malloc(ncoord[0]*sizeof(double));
              mlags[i]=(double *) Calloc(ncoord[0],double);
              if(mlags[i]==NULL) {*ismal=0; return;}}
    }
       SpaceTime_Dist(coordx,coordy,coordt,grid,ia,idx,ismal,ja,tapmodel,srange[1],trange[1],type[0]);
     //Set the range of the temporal intervals:
     trange[0]=*minimtime;
     *maxtime=*maximtime;
     if(trange[1]!=0) {if(trange[1]>*minimtime && trange[1]<*maximtime)
	*maxtime=trange[1];}
     else trange[1]=*maxtime;
    }


  // Set the range of the spatial distances:
  srange[0]=*minimdista;
  *maxdist=*maximdista;
  if(srange[1]!=0) {if(srange[1]>*minimdista && srange[1]<*maximdista)
		  *maxdist=srange[1];}
  else srange[1]=*maxdist;
  npair[0]=npairs[0];
  return;
  }

// Delete all the global variables
void DeleteGlobalVar()
{
    // Delete all the global variables:
    int i=0;
    if(!*isst) {// spatial case
         if(*istap) { Free(lags);}         // tapering case
         else {
              for(i=0;i<ncoord[0];i++)  {Free(mlags[i]);}
              Free(mlags);}}
    else {// spatialtemporal  case
          if(*istap) {Free(lags);Free(lagt);}  // tapering case
          else  {
            for(i=0;i<ncoord[0];i++) { Free(mlags[i]);}
            Free(mlags);
            for(i=0;i<ntime[0];i++) { Free(mlagt[i]);}
            Free(mlagt);}
    Free(maxtime);
    Free(maximtime);
    Free(minimtime);
    Free(ntime);
    }

    Free(tapsep);
    Free(npairs);
    Free(nrep);

    Free(maxdist);
    Free(minimdista);
    Free(maximdista);
    Free(ncoord);
    Free(ncoordx);
    Free(ncoordy);

    Free(istap);
    Free(isst);
    return;
}
