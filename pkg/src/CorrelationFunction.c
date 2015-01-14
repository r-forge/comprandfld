#include "header.h"

// check the validity of the parameters' range:
double CheckCor(int *cormod, double *par)
{
  double power=0.0, power1=0.0, power2=0.0, power_s=0.0, power_t=0.0;
  double rho=0.0, sep=0, scale=0.0, smooth=0.0, scale_s=0.0, scale_t=0;

  switch(*cormod) // Correlation functions are in alphabetical order
    {
    case 1:// Cauchy correlation function
      power2=par[0];
      scale=par[1];
      if(scale<=0 || power2<=0) rho=-2;
      break;
    case 2:// Exponential correlation function
    case 3:// Gaussian correlation function
    case 5:// Sferical correlation function
      scale=par[0];
      if(scale<=0) rho=-2;
      break;
    case 4: // Generalised Cuachy correlation function
      power1=par[0];
      power2=par[1];
      scale=par[2];
      if(scale<=0 || power1<=0 || power1>2 ||power2<=0) rho=-2;
      break;
    case 6:// Stable correlation function
      power=par[0];
      scale=par[1];
      if(scale<=0 || power<0 || power>2) rho=-2;
      break;
    case 7://  Whittle-Matern correlation function
      scale=par[0];
      smooth=par[1];
      if(scale<=0 || smooth<=0) rho=-2;
      break;
      // START non-separable correlation functions:
    case 21: //Gneiting correlation model as in (14) Gneitint (2002) with tau=1
    case 23://Porcu model as in (4) of Porcu Bevilacqua et al (2010), with beta_1=beta_2=1
    case 25:
    case 26: //Gneiting correlation on the sphere
    case 30:
      power_s=par[0];
      power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      sep=par[4];
      if(scale_s<=0 || scale_t<=0 || power_s<0 || power_s>2 || power_t<0 || power_t>2 || sep<0 || sep>1) rho=-2;
      break;
  case 22:// Iaco-Cesare model as in (14) of Gneitint (2006): note that power parameters are in [0,1]
      power2=par[0];
      power_s=par[1];
      power_t=par[2];
      scale_s=par[3];
      scale_t=par[4];
      if(power2<=0||scale_s<=0 || scale_t<=0 || power_s<0 || power_s>2 || power_t<0 || power_t>2) rho=-2;
      break;
    case 24:// Stein model as in (16) of JASA (2005) with epsilon=0*/
      power_t=par[0];
      scale_s=par[1];
      scale_t=par[2];
      smooth=par[3];
      if(scale_s<=0 || scale_t<=0 || power_t<0 || power_t>2 || smooth<=0) rho=-2;
      break;
      // END non-separable correlation functions
      // START separable correlation functions:
    case 41:// Exp-Cauchy:
      power2=par[0];
      scale_s=par[1];
      scale_t=par[2];
      if(scale_s<=0 || scale_t<=0 || power2<=0) rho=-2;
      break;
    case 42:// Double exp:
    case 43:// Exp-Gauss:
      scale_s=par[0];
      scale_t=par[1];
      if(scale_s<=0 || scale_t<=0) rho=-2;
      break;
    case 45:// Matern-Cauchy:
      power2=par[0];
      scale_s=par[1];
      scale_t=par[2];
      smooth=par[3];
      if(scale_s<=0 || scale_t<=0 || power2<=0 || smooth<=0) rho=-2;
      break;
    case 46:// Matern-exp:
      scale_s=par[0];
      scale_t=par[1];
      smooth=par[2];
      if(scale_s<=0 || scale_t<=0 || smooth<=0) rho=-2;
      break;
    case 47:// Stable-stable:
      power_s=par[0];
      power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      if(scale_s<=0 || scale_t<=0 || power_s<0 || power_s>2 || power_t<0 || power_s>2) rho=-2;
      break;
      // END separable correlation functions:
    }
  return rho;
}
// list of spatial and spatial-temporal correlation functions:
double CorFct(int *cormod, double h, double u, double *par)
{
  double arg=0.0, power=0.0, power1=0.0, power2=0.0, power_s=0.0, power_t=0.0;
  double rho=0.0, sep=0, scale=0.0, smooth=0.0, scale_s=0.0, scale_t=0, x=0;

  switch(*cormod) // Correlation functions are in alphabetical order
    {
    case 1:// Cauchy correlation function
      power1=2;
      power2=par[0];
      scale=par[1];
      rho=CorFunCauchy(h, power2, scale);
      break;
    case 2:// Exponential correlation function
      power=1;
      scale=par[0];
      rho=CorFunStable(h, power, scale);
      break;
    case 3:// Gaussian correlation function
      power=2;
      scale=par[0];
      rho=CorFunStable(h, power, scale);
      break;
    case 4: // Generalised Cuachy correlation function
      power1=par[0];
      power2=par[1];
      scale=par[2];
      rho=CorFunGenCauchy(h, power1, power2, scale);
      break;
    case 5:// Sferical correlation function
      scale=par[0];
      rho=CorFunSferical(h, scale);
      break;
    case 6:// Stable correlation function
      power=par[0];
      scale=par[1];
      rho=CorFunStable(h, power, scale);
      break;
    case 7://  Whittle-Matern correlation function
      scale=par[0];
      smooth=par[1];
      rho=CorFunWitMat(h, scale, smooth);
      break;
 /***************** spatial tapers****************************/
   case 15:// Wendland1 for tap
      rho=CorFunWend1(h,maxdist[0]);
      break;
    case 16:// Wendland1 for tap
      rho=CorFunWend2(h,maxdist[0]);
      break;
    case 17:// Wendland1 for tap
      rho=CorFunWend3(h,maxdist[0]);
 /***************** end spatial tapers****************************/
      // START non-separable correlation functions:
    case 21:   //Gneiting correlation model as in (14) Gneitint (2002) with tau=1
      power_s=par[0];
      power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      sep=par[4];
      arg=1+pow(u/scale_t, power_t);
      rho=exp(-(pow(h/scale_s, power_s))*pow(arg, -0.5*sep*power_s))/pow(arg,1);
      break;
    case 22:// Iaco-Cesare model as in (14) of Gneitint (2006): note that power parameters are in [0,1]
      power2=par[0];
      power_s=par[1];
      power_t=par[2];
      scale_s=par[3];
      scale_t=par[4];
      rho=pow(1+pow(h/scale_s, power_s)+pow(u/scale_t, power_t),-power2);
      break;
    case 23://Porcu model as in (4) of Porcu Bevilacqua et al (2010), with beta_1=beta_2=1
      power_s=par[0];
      power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      sep=par[4];
      if(sep) rho=pow(0.5*pow(1+pow(h/scale_s, power_s),sep)+0.5*pow(1+pow(u/scale_t, power_t),sep),-1/sep);
      else rho=pow((1+pow(h/scale_s, power_s))*(1+pow(u/scale_t,power_t)),-1);
      break;
    case 24:// Stein model as in (16) of JASA (2005) with epsilon=0*/
      power_t=par[0];
      scale_s=par[1];
      scale_t=par[2];
      smooth=par[3];
      arg=smooth+R_pow(u, power_t)/scale_t;
      if(h==0) rho=1/(R_pow(2, arg)*gammafn(arg+1));
      else rho=R_pow(h/scale_s, arg)*bessel_k(h/scale_s, arg, 1)/(R_pow(2, arg)*gammafn(arg+1));
      break;
    case 25:   //Gneiting correlation with prac ranges "giusto"
      power_s=par[0];
      power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      sep=par[4];
      arg=1+pow(u, power_t)/scale_t;
      rho=exp(-(pow(h, power_s)/scale_s)*pow(arg, 0.5*sep*power_s))/pow(arg,1.5);
      break;
     case 26:// Gneiting correlation model valid on the sphere
      power_s=par[0];
      power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      sep=par[4];
      arg=1+pow(h/scale_s, power_s);
      rho=exp(-pow(u/scale_t, power_t)/(pow(arg, 0.5*sep*power_t)))/arg;
      break;
    case 30:// non sep compact suppo
      power_s=par[0];
      power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      sep=par[4];
      arg=pow((1+(pow(u/scale_t,power_t))),0.5);
      if(h<(scale_s/pow(arg,sep))) rho=pow(1-h*pow(arg,sep)/scale_s,1.5+power_s)/arg;
      else                   rho=0;
      break;
      // END non-separable correlation functions
      // START separable correlation functions:
    case 41:// Exp-Cauchy:
      power2=par[0];
      scale_s=par[1];
      scale_t=par[2];
      rho=CorFunStable(h,1,scale_s)*pow((1+pow(u/scale_t, 2)), -power2);
      break;
    case 42:// Double exp:
      scale_s=par[0];
      scale_t=par[1];
      rho=CorFunStable(h,1,scale_s)*CorFunStable(u,1,scale_t);
      break;
    case 43:// Exp-Gauss:
      scale_s=par[0];
      scale_t=par[1];
      rho=CorFunStable(h,1,scale_s)*CorFunStable(u,2,scale_t);
      break;
    case 45:// Matern-Cauchy:
      power2=par[0];
      scale_s=par[1];
      scale_t=par[2];
      smooth=par[3];
      if(h==0) arg=1;
      else arg=pow(2,1-smooth)/gammafn(smooth)*pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);
      rho=arg*pow((1+pow(u/scale_t, 2)),-power2);
      break;
    case 46:// Matern-exp:
      scale_s=par[0];
      scale_t=par[1];
      smooth=par[2];
      if(h==0) arg=1;
      else arg=pow(2,1-smooth)/gammafn(smooth)*pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);
      rho=arg*exp(-u/scale_t);
      break;
    case 47:// Stable-stab:
      power_s=par[0];
      power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      rho=CorFunStable(h,power_s,scale_s)*CorFunStable(u,power_t,scale_t);
      break;
    /******************************space time taper***************************************/
    case 100:
      scale_s=par[0];
      scale_t=par[1];
      rho=CorFunWend1(h,maxdist[0])*CorFunWend2(u,maxtime[0]);
      break;
        case 101:
      scale_s=par[0];
      scale_t=par[1];
      rho=CorFunWend1(h,maxdist[0])*CorFunWend2(u,maxtime[0]);
      break;
        case 102:
      scale_s=par[0];
      scale_t=par[1];
      rho=CorFunWend1(h,maxdist[0])*CorFunWend3(u,maxtime[0]);
      break;

      case 103:
      scale_s=par[0];
      scale_t=par[1];
      rho=CorFunWend2(h,maxdist[0])*CorFunWend2(u,maxtime[0]);
      break;
        case 104:
      scale_s=par[0];
      scale_t=par[1];
      rho=CorFunWend2(h,maxdist[0])*CorFunWend2(u,maxtime[0]);
      break;
        case 105:
      scale_s=par[0];
      scale_t=par[1];
      rho=CorFunWend2(h,maxdist[0])*CorFunWend3(u,maxtime[0]);
      break;
          case 106:
      scale_s=par[0];
      scale_t=par[1];
      rho=CorFunWend3(h,maxdist[0])*CorFunWend2(u,maxtime[0]);
      break;
        case 107:
      scale_s=par[0];
      scale_t=par[1];
      rho=CorFunWend3(h,maxdist[0])*CorFunWend2(u,maxtime[0]);
      break;
        case 108:
      scale_s=par[0];
      scale_t=par[1];
      rho=CorFunWend3(h,maxdist[0])*CorFunWend3(u,maxtime[0]);
     break;
     case 109:  /* non separable temporalquasi taper */
      scale_s=*maxdist;
      arg=pow(1+pow(h/scale_s, 1),*tapsep);
      x=u*arg/(*maxtime);
      if(u<=*maxtime/arg) rho=pow(arg,-6)*(1+7*x)*pow(1-x,7);
      else  rho=0;
     break;
     case 110:  /* non separable spatial quasi taper */
     scale_t=*maxtime;
      arg=pow(1+pow(u/scale_t, 1),*tapsep);
      x=h*arg/(*maxdist);
      if(h<= *maxdist/arg) rho=pow(arg,-6)*(1+7*x)*pow(1-x,7);
      else  rho=0;

     break;
  /******************************end space time taper***************************************/
      // END separable correlation functions:
    }
  return rho;
}
// Cauhcy class of correlation models:
double CorFunCauchy(double lag, double power2, double scale)
{
  double rho=0.0;
  // Computes the correlation:
  rho=pow((1+pow(lag/scale,2)),-power2);
  return rho;
}
// Generalised Cauhcy class of correlation models:
double CorFunGenCauchy(double lag, double power1, double power2, double scale)
{
  double rho=0.0;
  // Computes the correlation:
  rho=pow((1+pow(lag/scale,power1)), -power2/power1);
  return rho;
}
// Stable class of correlation models:
double CorFunStable(double lag, double power, double scale)
{
  double rho=0.0;
  // Computes the correlation:
  rho=exp(-pow(lag/scale,power));
  return rho;
}
// Double Stable class of correlation models:
/*double CorFunDobStable(double lag, double power_s, double power_t, double scale_s, double scale_t, double tsep)
{
  double rho=0.0;
  // Computes the correlation:
  rho=exp(-pow(lag/scale_s,power_s)-pow(tsep/scale_t,power_t));
  return rho;
  }*/
// Sferical class of correlation models:
double CorFunSferical(double lag, double scale)
{
  double rho=0.0;
  if(lag<=scale) rho=1-1.5*lag/scale+0.5*pow(lag/scale, 3);
  else rho=0;
  return rho;
}
// Whittle=matern class of correlation models:
double CorFunWitMat(double lag, double scale, double smooth)
{
  double rho=0.0;
  // Computes the correlation:
  if(lag) rho=pow(2,1-smooth)/gammafn(smooth)*pow(lag/scale,smooth)*
	bessel_k(lag/scale,smooth,1);
  else rho=1;
  return rho;
}
double CorFunWend1(double lag,double scale)
{
  double rho=0.0,x=0;
  x=lag*pow(scale,-1);
  if(x<=1) rho=pow(1-x,2)*(1+0.5*x);
  else rho=0;
  return rho;
}
double CorFunWend2(double lag,double scale)
{
  double rho=0.0,x=0;
   x=lag*pow(scale,-1);
  if(x<=1) rho=pow(1-x,4)*(1+4*x);
  else rho=0;
  return rho;
}
double CorFunWend3(double lag,double scale)
{
  double rho=0.0,x=0;
   x=lag*pow(scale,-1);
  if(x<=1) rho=pow(1-x,6)*(1+6*x+35*pow(x,2)/3);
  else rho=0;
  return rho;
}
// Computation of the upper (lower) triangular spatial correlation matrix:
void CorrelationMat(double *rho, int *cormod, double *nuis, double *par)
{
    
  int i=0,j=0,h=0;// check the paramaters range:
  if(nuis[1]<0 || nuis[2]<=0 || CheckCor(cormod,par)==-2){
    rho[0]=-2;
    return;}// compute the correlations:
     for(i=0;i<(*ncoord-1);i++){
	    for(j=(i+1);j<*ncoord;j++){
    rho[h]=CorFct(cormod,mlags[i][j],0,par);
           
    h++;
    }}
   
  return;
}

// Computation of the correlations for spatial tapering:
void CorrelationMat_tap(double *rho, int *cormod, double *nuis, double *par)
{
  int i=0;// check the paramaters range:
  if(nuis[1]<0 || nuis[2]<=0 || CheckCor(cormod,par)==-2){
    rho[0]=-2;
    return;}// compute the correlations:
  for(i=0;i<*npairs;i++) rho[i]=CorFct(cormod,lags[i],0,par);
  return;
}
// Computation of the upper (lower) triangular spatial-temporal correlation matrix:
void CorrelationMat_st(double *rho, int *cormod, double *nuis, double *par)
{
  int i=0,j=0,t=0,v=0,h=0;
  if(nuis[1]<0 || nuis[2]<=0 || CheckCor(cormod,par)==-2){
    rho[0]=-2;
    return;}
  for(i=0;i<*ncoord;i++){
    for(t=0;t<*ntime;t++){
      for(j=i;j<*ncoord;j++){
	if(i==j){
	  for(v=t+1;v<*ntime;v++){
    rho[h]=CorFct(cormod,0,mlagt[t][v],par);
    h++;}}
    else {
         for(v=0;v<*ntime;v++){
              rho[h]=CorFct(cormod,mlags[i][j],mlagt[t][v],par);
              h++;}}
    }}}
  return;
}

void CorrelationMat_st_tap(double *rho, int *cormod, double *nuis, double *par)
{
  int i=0;
  if(nuis[1]<0 || nuis[2]<=0 || CheckCor(cormod,par)==-2){
    rho[0]=-2;
    return;}
  // first spatial loop:
  for(i=0;i<*npairs;i++) rho[i]=CorFct(cormod,lags[i],lagt[i],par);

  return;
}

void Corr_c(double *cc,double *coordx, double *coordy, double *coordt, int *cormod, int *grid, double *locx,  double *locy,int *ncoord, int *nloc,int*tloc,
                 int *ntime, double *par, int *spt, double *time,int *type)
{
int i,j,h=0;
double dis=0.0;
if(!spt[0])	{   //spatial case
switch(*type){
case 0:// Euclidean distances:
	if(*grid){// in case of a equispaced grid of coordinates:
	         }
	else{// in case of an irregular grid of coordinates:
     for(j=0;j<(*nloc);j++){
	  for(i=0;i<(*ncoord);i++){
	      dis=hypot(coordx[i]-locx[j],coordy[i]-locy[j]);
	      cc[h]=CorFct(cormod,dis,0,par);
	      h++;
	      }}}
	      break;
case 1:// Chordal distances:
	if(*grid){// in case of a equispaced grid of coordinates:
	         }
	else{// in case of an irregular grid of coordinates:
	      for(j=0;j<(*nloc);j++){
	  for(i=0;i<(*ncoord);i++){
	      dis=Dist_chordal(coordx[i],coordy[i],locx[j],locy[j]);
	      cc[h]=CorFct(cormod,dis,0,par);h++;
	        }}}
	      break;
case 2:// Geodesic distances:
	 	if(*grid){// in case of a equispaced grid of coordinates:
	          }
	else{// in case of an irregular grid of coordinates:
	      for(j=0;j<(*nloc);j++){
	  for(i=0;i<(*ncoord);i++){
	      dis=Dist_geodesic(coordx[i],coordy[i],locx[j],locx[j]);
	      cc[h]=CorFct(cormod,dis,0,par);h++;
	     }}}
	      break;
	}}
else{    //spatio temporal  case

int t,v;double dit=0.0;
switch(*type)
    {
           case 0:// Euclidean distances:
	     if(*grid){  // in case of a equispaced grid of coordinates:
	      }
	      else{ // in case of an irregular grid of coordinates:
	            for(v=0;v<(*tloc);v++){
	           for(j=0;j<(*nloc);j++){
	      for(i=0;i<*ncoord;i++){
		  	dis=hypot(coordx[i]-locx[j],coordy[i]-locy[j]);
	      for(t=0;t<*ntime;t++){
		    dit=fabs(coordt[t]-time[v]);
		    cc[h]=CorFct(cormod,dis,dit,par);
		    h++;}}}}
		  }
		    break;
		  case 1: // Chordal distances:
		  	if(*grid){ // in case of a equispaced grid of coordinates:
	         }
	         else{ // in case of an irregular grid of coordinates:
	               for(v=0;v<(*tloc);v++){
	              for(j=0;j<(*nloc);j++){
	       	for(i=0;i<*ncoord;i++){
		    	dis=Dist_chordal(coordx[i],coordy[i],locx[j],locy[j]);
	        for(t=0;t<*ntime;t++){
		      dit=fabs(coordt[t]-time[v]);
		      cc[h]=CorFct(cormod,dis,dit,par);
		      h++;}}}}
		      }
	       break;
		  case 2:// Geodesic distances:
		  	if(*grid){ // in case of a equispaced grid of coordinates:
	         }
	          else{ // in case of an irregular grid of coordinates:
	                 for(v=0;v<(*tloc);v++){
                      for(j=0;j<(*nloc);j++){
	         	for(i=0;i<*ncoord;i++){
                  dis=Dist_geodesic(coordx[i],coordy[i],locx[j],locy[j]);
	          for(t=0;t<*ntime;t++){
		        dit=fabs(coordt[t]-time[v]);
		        cc[h]=CorFct(cormod,dis,dit,par);
		        h++;}} }}
	         }
	         break;
	}}}


void Corr_c_tap(double *cc,double *cc_tap,double *coordx, double *coordy, double *coordt, int *cormod, int *cormodtap, int *grid, double *locx,  double *locy,int *ncoord, int *nloc,int*tloc,
                 int *ntime, double *par, int *spt, double *time,int *type)
{
int i,j,h=0;
double dis=0.0;
if(!spt[0])	{   //spatial case

switch(*type){
case 0:// Euclidean distances:
	if(*grid){// in case of a equispaced grid of coordinates:
	         }
	else{// in case of an irregular grid of coordinates:
     for(j=0;j<(*nloc);j++){
	  for(i=0;i<(*ncoord);i++){
	      dis=hypot(coordx[i]-locx[j],coordy[i]-locy[j]);
	      cc[h]=CorFct(cormod,dis,0,par);
	      cc_tap[h]=CorFct(cormod,dis,0,par)*CorFct(cormodtap,dis,0,par);
	      h++;
	      }}}
	      break;
case 1:// Chordal distances:
	if(*grid){// in case of a equispaced grid of coordinates:
	         }
	else{// in case of an irregular grid of coordinates:
	      for(j=0;j<(*nloc);j++){
	  for(i=0;i<(*ncoord);i++){
	      dis=Dist_chordal(coordx[i],coordy[i],locx[j],locy[j]);
	      cc[h]=CorFct(cormod,dis,0,par);
	      cc_tap[h]=CorFct(cormod,dis,0,par)*CorFct(cormodtap,dis,0,par);
	      h++;
	        }}}
	      break;
case 2:// Geodesic distances:
	 	if(*grid){// in case of a equispaced grid of coordinates:
	          }
	else{// in case of an irregular grid of coordinates:
	      for(j=0;j<(*nloc);j++){
	  for(i=0;i<(*ncoord);i++){
	      dis=Dist_geodesic(coordx[i],coordy[i],locx[j],locx[j]);
	      cc[h]=CorFct(cormod,dis,0,par);
	      cc_tap[h]=CorFct(cormod,dis,0,par)*CorFct(cormodtap,dis,0,par);
	      h++;
	     }}}
	      break;
	}}
else{    //spatio temporal  case

int t,v;double dit=0.0;
switch(*type)
    {
           case 0:// Euclidean distances:
	     if(*grid){  // in case of a equispaced grid of coordinates:
	      }
	      else{ // in case of an irregular grid of coordinates:
	            for(v=0;v<(*tloc);v++){
	           for(j=0;j<(*nloc);j++){
	      for(i=0;i<*ncoord;i++){
		  	dis=hypot(coordx[i]-locx[j],coordy[i]-locy[j]);
	      for(t=0;t<*ntime;t++){
		    dit=fabs(coordt[t]-time[v]);
		    cc[h]=CorFct(cormod,dis,dit,par);
	         cc_tap[h]=cc[h]*CorFct(cormodtap,dis,dit,par);
		    h++;}}}}
		  }
		    break;
		  case 1: // Chordal distances:
		  	if(*grid){ // in case of a equispaced grid of coordinates:
	         }
	         else{ // in case of an irregular grid of coordinates:
	               for(v=0;v<(*tloc);v++){
	              for(j=0;j<(*nloc);j++){
	       	for(i=0;i<*ncoord;i++){
		    	dis=Dist_chordal(coordx[i],coordy[i],locx[j],locy[j]);
	        for(t=0;t<*ntime;t++){
		      dit=fabs(coordt[t]-time[v]);
		      cc[h]=CorFct(cormod,dis,dit,par);
	           cc_tap[h]=cc[h]*CorFct(cormodtap,dis,dit,par);
		      h++;}}}}
		      }
	       break;
		  case 2:// Geodesic distances:
		  	if(*grid){ // in case of a equispaced grid of coordinates:
	         }
	          else{ // in case of an irregular grid of coordinates:
	                 for(v=0;v<(*tloc);v++){
                      for(j=0;j<(*nloc);j++){
	         	for(i=0;i<*ncoord;i++){
                  dis=Dist_geodesic(coordx[i],coordy[i],locx[j],locy[j]);
	          for(t=0;t<*ntime;t++){
		        dit=fabs(coordt[t]-time[v]);
		        cc[h]=CorFct(cormod,dis,dit,par);
	             cc_tap[h]=cc[h]*CorFct(cormodtap,dis,dit,par);
		        h++;}} }}
	         }
	         break;
	}}}


// Derivatives with respect to power2 of the Cauchy correlation model:
double DCauchyPow(double power2, double scale, double rho)
{
  return -rho*log(pow(rho,-1/power2));
}
// Derivatives with respect to scale of the Cauchy correlation model:
double DCauchySc(double lag, double power2, double scale, double rho)
{
 return 2*power2*rho*pow(rho, 1/power2)*pow(lag,2)/pow(scale,3);
}
// Derivatives with respect to scale of the Exponential correlation model:
double DExpoSc(double lag, double scale, double rho)
{
 return rho*lag/pow(scale,2);
}
// Derivatives with respect to scale of the Gaussian correlation model:
double DGaussSc(double lag, double scale, double rho)
{
  return 2*rho*pow(lag,2)/pow(scale,3);
}
// Derivatives with respect to power1 of the generalised Cauchy correlation model:
double DGenCauP1(double lag, double power1, double power2, double scale, double rho)
{
  if(lag)return power2*rho/power1*(log(1+pow(lag/scale,power1))/power1-
				   pow(lag/scale,power1)*log(lag/scale)/
				   (1+pow(lag/scale,power1)));
  else return 0.0;
}
// Derivatives with respect to power2 of the generalised Cauchy correlation model:
double DGenCauP2(double lag, double power1, double scale, double rho)
{
  return -rho*log(1+pow(lag/scale,power1))/power1;
}
// Derivatives with respect to scale of the generalised Cauchy correlation model:
double DGenCauSc(double lag, double power1, double power2, double scale, double rho)
{
  if(lag) return rho/(1+pow(lag/scale,2))*power2*pow(lag,power1)/pow(scale,power1+1);
  else return 0.0;
}
// Derivatives with respect to scale of the sferical correlation model:
double DSferiSc(double lag, double scale)
{
  if(lag<=scale) return 1.5*lag*(1-pow(lag/scale, 2))/pow(scale, 2);
  else return 0.0;
}
// Derivatives with respect to power of the Stable correlation model:
double DStabPow(double lag, double power, double scale, double rho)
{
  if(lag) return -rho*pow(lag/scale,power)*log(lag/scale);
  else return 0.0;
}
// Derivatives with respect to scale of the Stable correlation model:
double DStabSc(double lag, double power, double scale, double rho)
{
  if(lag) return rho*power*pow(lag/scale,power)/scale;
  else return 0.0;
}
// Derivatives with respect to scale of the Whittle-Matern correlation model:
double DWhMatSc(double *eps, double lag, double scale, double smooth)
{
  if (lag){
    double pscale=0.0;
    pscale=(bessel_k(lag/(scale+*eps),smooth,1)-
	    bessel_k(lag/scale,smooth,1))/ *eps;
    return pow(2,1-smooth)/gammafn(smooth)*pow(lag/scale,smooth)*
      (pscale-smooth*bessel_k(lag/scale,smooth,1)/scale);}
  else return 0;
}

// Derivatives with respect to smooth of the Whittle-Matern correlation model:
double DWhMatSm(double *eps, double lag, double scale, double smooth)
{
  if (lag){
    double psmooth=0.0;
    psmooth=(bessel_k(lag/scale,smooth+ *eps,1)-
	     bessel_k(lag/scale,smooth,1))/ *eps;
    return pow(2,1-smooth)*pow(lag/scale,smooth)/gammafn(smooth)*
      ((log(lag/scale)-log(2)-digamma(smooth))*bessel_k(lag/scale,smooth,1)+psmooth);}
  else return 0;
}
// Derivatives with respect to the temporal scale parameter of the Gneiting correlation model:
double DGneiting_sc_t(double h,double u, double power_s,double power_t,
		      double scale_s,double scale_t,double sep)
{
  double arg=0.0,rho=0.0,a=0.0;
  arg=1+pow(u/scale_t, power_t);
  rho=exp(-pow(h/scale_s, power_s)/(pow(arg, 0.5*sep*power_s)))/arg;
  a=( power_t*pow(u/scale_t, power_t)*rho - 0.5*rho*pow(h/scale_s, power_s)*power_t*power_s*sep*pow(arg,-0.5*sep*power_s))/(scale_t*arg);
  return(a) ;
}
// Derivatives with respect to the spatial scale parameter of the Gneiting correlation model:
double DGneiting_sc_s(double h,double u,double power_s,double power_t,
		      double scale_s,double scale_t,double sep)
{
  double arg=0.0,rho=0.0;
  arg=1+pow(u/scale_t, power_t);
  rho=exp(-pow(h/scale_s, power_s)/(pow(arg, 0.5*sep*power_s)))/arg;
  return ((pow(h/scale_s, power_s)*power_s*rho*pow(arg,-0.5*power_s*sep))/scale_s);
}
// Derivatives with respect to the separable parameter of the Gneiting correlation model:
double DGneiting_sep(double h,double u, double power_s,double power_t,
		     double scale_s,double scale_t,double sep)
{
  double a=0,arg=0,rho=0;
  arg=1+pow(u/scale_t, power_t);
  rho=exp(-pow(h/scale_s, power_s)/(pow(arg, 0.5*sep*power_s)))/arg;
  if(arg) a=0.5*pow(h/scale_s, power_s)*pow(arg,-0.5*sep*power_s)*power_s*rho*log(arg);
  return(a);
}
// Derivatives with respect to spatial power of the Gneiting correlation model:
double DGneiting_pw_s(double h,double u, double power_s,double power_t,
		      double scale_s,double scale_t,double sep)
{
  double arg=0.0,rho=0.0,a=0.0;
  arg=1+pow(u/scale_t, power_t);
  rho=exp(-pow(h/scale_s, power_s)/(pow(arg, 0.5*sep*power_s)))/arg;
  if(h && arg){a=(-pow(h/scale_s, power_s)*pow(arg,-0.5*sep*power_s)*log(h/scale_s) +
                  0.5*pow(h/scale_s, power_s)*pow(arg,-0.5*sep*power_s)*sep*log(arg) )*rho;}
  return(a);
}

// Derivatives with respect to temporal power of the Gneiting correlation model:
double DGneiting_pw_t(double h,double u, double power_s,double power_t,
		      double scale_s,double scale_t,double sep)
{
  double arg=0.0,rho=0.0,a=0.0;
  arg=1+pow(u/scale_t, power_t);
  rho=exp(-pow(h/scale_s, power_s)/(pow(arg, 0.5*sep*power_s)))/arg;
  if(u){a=( -pow(u/scale_t, power_t)*log(u/scale_t)*rho+
             0.5*rho*pow(h/scale_s, power_s)*power_s*sep*log(u/scale_t)*pow(arg,-0.5*sep*power_s))/arg;}
  return(a) ;
}
// Derivatives with respect to the separable parameter of the Gneiting correlation model:
double DIaco_sc_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double power2)
{
 double rho=0,arg=0;
 arg=1+pow(h/scale_s, power_s)+pow(u/scale_t, power_t);
 rho=pow(arg,-power2);
 return (rho*power_s*power2*pow(h/scale_s, power_s))/(arg*scale_s);
 }

// Derivatives with respect to temporal scale of the Gneiting GC  correlation model:
double DGneiting_GC_sc_t(double h,double u, double power_s,double power_t,
		      double scale_s,double scale_t,double sep)
{
  double arg=0.0,rho=0.0;
      arg=1+pow(h/scale_s, power_s);
      rho=exp(-pow(u/scale_t, power_t)/(pow(arg, 0.5*sep*power_t)))/arg;
      return(pow(u,power_t)*rho*pow(arg,-0.5*sep*power_t-1)/(scale_t*scale_t));
}

// Derivatives with respect to the separable parameter of the Gneiting GC  correlation model:
double DGneiting_GC_sep(double h,double u, double power_s,double power_t,
		     double scale_s,double scale_t,double sep)
{
    double arg=0.0,rho=0.0;
    arg=1+pow(h/scale_s, power_s);
    rho=exp(-pow(u/scale_t, power_t)/(pow(arg, 0.5*sep*power_t)))/arg;
    return(0.5*log(arg)*pow(u,power_t)*power_t*rho*pow(arg,-0.5*sep*power_t-1)/(scale_t));
  }

// Derivatives with respect to temporal power of the Gneiting GC correlation model:
double DGneiting_GC_pw_t(double h,double u, double power_s,double power_t,
		      double scale_s,double scale_t,double sep)
{
double arg=0.0,rho=0.0;
    arg=1+pow(h/scale_s, power_s);
    rho=exp(-pow(u/scale_t, power_t)/(pow(arg, 0.5*sep*power_t)))/arg;
    if(u&&arg) return(-0.5*pow(u,power_t)*rho*pow(arg,-0.5*sep*power_t)*(2*log(u)-sep*log(arg))/(arg*scale_t));
    else return 0;
}

// Derivatives with respect to spatial scale of the Gneiting  GC correlation model:
double DGneiting_GC_sc_s(double h,double u,double power_s,double power_t,
		      double scale_s,double scale_t,double sep)
{
  double arg=0.0,rho=0.0,a=0.0,b=0.0;
  arg=1+pow(h/scale_s, power_s);
  rho=exp(-pow(u/scale_t, power_t)/(pow(arg, 0.5*sep*power_t)))/arg;
  a=-0.5*(pow(u,power_t)*sep*power_t*power_s*pow(h/scale_s, power_s)*rho)/(scale_s*scale_t*pow(arg,2)*pow(arg,0.5*sep*power_t));
  b=(power_s*pow(h/scale_s, power_s)*rho)/(scale_s*pow(arg,2));
  return(a+b);
}
// Derivatives with respect to spatial power of the Gneiting correlation model:
double DGneiting_GC_pw_s(double h,double u, double power_s,double power_t,
		      double scale_s,double scale_t,double sep)
{
 {
  double arg=0.0,rho=0.0,a=0.0,b=0.0;
  arg=1+pow(h/scale_s, power_s);
  rho=exp(-pow(u/scale_t, power_t)/(pow(arg, 0.5*sep*power_t)))/arg;
  if(h) a=0.5*(log(h/scale_s)*pow(u,power_t)*sep*power_t*pow(h/scale_s, power_s)*rho)/(scale_t*pow(arg,2)*pow(arg,0.5*sep*power_t));
  else  a=0;
  if(h) b=pow(h/scale_s, power_s)*(log(h/scale_s)*rho)/(pow(arg,2));
  else b=0;
  return(a-b);
}

}

 double DIaco_sc_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double power2)
{

  double rho=0,arg=0;
  arg=1+pow(h/scale_s, power_s)+pow(u/scale_t, power_t);
  rho=pow(arg,-power2);
  return (rho*power_t*power2*pow(u/scale_t, power_t))/(arg*scale_t);

 }
 double DIaco_pw_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double power2)
{
  double a=0,rho=0,arg=0;
  arg=1+pow(h/scale_s, power_s)+pow(u/scale_t, power_t);
  rho=pow(arg,-power2);
  if(h) a=-power2*log(h/scale_s)*pow(h/scale_s, power_s)*rho/arg;
  return a;

 }
 double DIaco_pw_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double power2)
{
  double a=0,rho=0,arg=0;
 arg=1+pow(h/scale_s, power_s)+pow(u/scale_t, power_t);
 rho=pow(arg,-power2);
 if(h) a=-power2*log(u/scale_t)*pow(u/scale_t, power_t)*rho/arg;
 return a;
 }

  double DIaco_pw2(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double power2)
{
 double a=0,rho=0,arg=0;
 arg=1+pow(h/scale_s, power_s)+pow(u/scale_t, power_t);
 rho=pow(arg,-power2);
 if(arg) a=-log(arg)*rho;
 return(a);
 }

double DMat_Cauchy_sc_t(double h,double u,double power2,double scale_s,double scale_t,double smooth)
{
  double arg=0,arg3=0;
  arg3=pow((1+pow(u/scale_t, 2)),-power2);
  if(h) arg=(pow(2, 1-smooth)/gammafn(smooth))*pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);
  else arg=1;
  return 2*pow(u,2)*power2*arg*arg3/(pow(scale_t,3)*(1+pow(u/scale_t, 2)));
}

double DMat_Cauchy_pw2(double h,double u,double power2,double scale_s,double scale_t,double smooth)
{
  double arg=0.0,arg2=0.0,arg3=0.0;
  arg3=pow((1+pow(u/scale_t, 2)),-power2);
  if(h) arg=(pow(2, 1-smooth)/gammafn(smooth))*pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);
  else  arg=1;
  if(1+pow(u/scale_t, 2)) arg2=arg*arg3*log(1+pow(u/scale_t, 2));
 return arg2;
}


double DMat_Cauchy_sc_s(double h,double u,double power2,double scale_s,double scale_t,double smooth)
{
  double arg1,arg2=0,arg3;
  arg3=pow((1+pow(u/scale_t, 2)),-power2);
  if(h) {  arg2=2*smooth*scale_s*bessel_k(h/scale_s, smooth,1)-h*bessel_k(h/scale_s, smooth+1,1);}
  arg1=pow(2, 1-smooth)*pow(h/scale_s, smooth)*arg3;
  return -arg1*arg2/(gammafn(smooth)*pow(scale_s,2));
}

double DMat_Cauchy_sm(double h,double u,double *eps, double power2, double scale_s,double scale_t,double smooth)
{
  double arg=0.0,arg2=0.0,arg3=0.0,psmooth=0.0;
  arg3=pow((1+pow(u/scale_t, 2)),-power2);
  psmooth=(bessel_k(h/scale_s,smooth+ *eps,1)-bessel_k(h/scale_s,smooth,1))/ *eps;
  if(h) arg=(pow(2, 1-smooth)/gammafn(smooth))*pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);
  else arg=1;
  if(h) arg2=-arg3*arg*(log(2)+digamma(smooth)-log(h/scale_s)-psmooth/bessel_k(h/scale_s, smooth, 1));
  else arg2=0;
  return arg2;
}

double DMat_Exp_sc_t(double h,double u,double scale_s,double scale_t,double smooth)
{
  double arg=0;
  if(h) arg=(pow(2, 1-smooth)/gammafn(smooth))*pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);
 else arg=1;
  return arg*u*exp(-u/scale_t)/R_pow(scale_t,2);
}

double DMat_Exp_sc_s(double h,double u,double scale_s,double scale_t,double smooth)
{
double arg1=0,arg2=0;

 if(h){arg1=(pow(2, 1-smooth)/gammafn(smooth))*pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);
   arg2=h*(pow(2, 1-smooth)/gammafn(smooth))*pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth+1, 1);}
 else {arg1=1;
   arg2=0;}

 return -arg1*2*smooth*exp(-u/scale_t)/scale_s + arg2*exp(-u/scale_t)/pow(scale_s,2);
}

double DMat_Exp_sm(double h,double u,double *eps,double scale_s,double scale_t,double smooth)
{
  double arg=0.0,psmooth=0.0,arg2=0.0;
  psmooth=(bessel_k(h/scale_s,smooth+ *eps,1)-bessel_k(h/scale_s,smooth,1))/ *eps;
  if(h){arg=(pow(2, 1-smooth)/gammafn(smooth))*pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);
    arg2=-exp(-u/scale_t)*arg*(log(2)+digamma(smooth)-log(h/scale_s)-psmooth/bessel_k(h/scale_s,smooth,1));}
  else arg2=0;
  return arg2;
}

double DPorcu_sc_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep)
{
  double arg=0,arg1=0,arg2=0;
  arg1=1+pow(h/scale_s, power_s);
  arg2=1+pow(u/scale_t, power_t);
  arg=0.5*pow(arg1, sep)+0.5*pow(arg2,sep);
  if(!sep) return((pow(h/scale_s,power_s)*power_s)/(arg1*arg1*arg2*scale_s));
  else    return((pow(arg,-(sep+1)/sep)*pow(arg1,sep-1)*pow(h/scale_s, power_s)*power_s)/(2*scale_s));
}
double DPorcu_sc_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep)
{
double arg=0,arg1=0,arg2=0;
  arg1=1+pow(h/scale_s, power_s);
  arg2=1+pow(u/scale_t, power_t);
  arg=0.5*pow(arg1, sep)+0.5*pow(arg2,sep);
  if(!sep) return((pow(u/scale_t, power_t)*power_t)/(arg1*arg2*arg2*scale_t));
  else    return((pow(arg,-(sep+1)/sep)*pow(arg1,sep-1)*pow(u/scale_t, power_t)*power_t)/(2*scale_t));
}
double DPorcu_pw_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep)
{
  double a=0,arg=0,arg1=0,arg2=0,rho=0;
  arg1=1+pow(h/scale_s, power_s);
  arg2=1+pow(u/scale_t, power_t);
  arg=pow(0.5*pow(arg1, sep)+0.5*pow(arg2,sep),-1/sep);
  if(sep) rho=pow(arg,-1/sep);
  else rho=pow(arg1*arg2,-1);
  if(h) a=-(0.5*pow(h, power_s)*rho*arg1*log(h))/(arg*arg1*scale_s);
  return a;
}
double DPorcu_pw_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep)
{
  double a=0,arg=0,arg1=0,arg2=0,rho=0;
  arg1=1+pow(h/scale_s, power_s);
  arg2=1+pow(u/scale_t, power_t);
  arg=pow(0.5*pow(arg1, sep)+0.5*pow(arg2,sep),-1/sep);
  if(sep) rho=pow(arg,-1/sep);
  else rho=pow(arg1*arg2,-1);
  if(h) a=-(0.5*pow(u, power_t)*rho*arg2*log(u))/(arg*arg2*scale_t);
  return a;
}

double DPorcu_sep(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep)
{
  double arg=0,arg1=0,arg2=0,rho=0,a=0;
   arg1=1+pow(h/scale_s, power_s);
  arg2=1+pow(u/scale_t, power_t);
  arg=pow(0.5*pow(arg1, sep)+0.5*pow(arg2,sep),-1/sep);
  if(sep) rho=pow(arg,-1/sep);
  else rho=pow(arg1*arg2,-1);
  if(arg1&&arg2)
    a=rho*(log(arg)/pow(sep,2)-((0.5*log(arg1)*pow(arg1,sep)+0.5*log(arg2)*pow(arg2,sep))))/(sep*arg);
  return a;
}
/*
double DExp_Exp_sc_s(double h,double u,double scale_s,double scale_t)
{
  double rho=0;
  rho=exp(-h/scale_s-u/scale_t);
  return h*rho/pow(scale_s,2);
}

double DExp_Exp_sc_t(double h,double u,double scale_s,double scale_t)
{
  double rho=0;
  rho=exp(-h/scale_s-u/scale_t);
  return u*rho/pow(scale_t,2);
}
*/
double DExp_Cauchy_sc_s(double h,double u,double power2,double scale_s,double scale_t)
{
  double rho=0;
  rho=exp(-h/scale_s)*pow(1+pow(u/scale_t, 2), -power2);
  return h*rho/pow(scale_s,2);
 }
 double DExp_Cauchy_sc_t(double h,double u,double power2,double scale_s,double scale_t)
{
  double rho=0;
  rho=exp(-h/scale_s)*pow(1+pow(u/scale_t, 2), -power2);
  return 2*rho*pow(u,2)*power2/(pow(scale_t,3)*(1+pow(u/scale_t,2)));
 }

  double DExp_Cauchy_pw2(double h,double u,double power2,double scale_s,double scale_t)
{
  double a=0,rho=0;
  rho=exp(-h/scale_s)*pow(1+pow(u/scale_t, 2), -power2);
  if(1+pow(u/scale_t,2))   a=-log(1+pow(u/scale_t,2))*rho;
  return a;
 }
/*
double DExp_Gauss_sc_s(double h,double u,double scale_s,double scale_t)
{
  double rho=0;
  rho=exp(-h/scale_s-pow(u/scale_t,2));
  return h*rho/pow(scale_s,2);
}
double DExp_Gauss_sc_t(double h,double u,double scale_s,double scale_t)
{
  double rho=0;
  rho=exp(-h/scale_s-pow(u/scale_t,2));
  return 2*pow(u,2)*rho/pow(scale_t,3);
  }*/
// compute the gradient matrix (numcoord...) for the spatial field:
void DCorrelationMat(int *cormod,double *drho,double *eps,int *flagcor,int *nparcor,double *parcor,double *rho)
{
  int i=0,j=0,m=0,k=0,h=0,npa=0;
  double *gradcor,*derho;
  //Initializes gradients:
   npa=0.5* *ncoord * (*ncoord-1);
  gradcor=(double *) R_alloc(*nparcor, sizeof(double));
  derho=(double *) R_alloc(npa * *nparcor, sizeof(double));

  k=0;
for(i=0; i<(*ncoord-1);i++){
    for(j=(i+1); j<*ncoord;j++){
    GradCorrFct(rho[h],cormod,eps,flagcor,gradcor,mlags[i][j],0,parcor);
    for(m=0;m<*nparcor;m++){
      derho[k]=gradcor[m];k++;}h++;}}
  k=0;
  for(i=0;i<*nparcor;i++){
    for(m=0;m<npa;m++){
      drho[k]=derho[i+m* *nparcor];k++;}}
  return;
}
// compute the gradient matrix (numcoord...) for the spatial field:
void DCorrelationMat_tap(int *cormod,double *drho,double *eps,int *flagcor,int *nparcor,double *parcor,double *rho)
{
  int i=0,m=0,k=0;
  double *gradcor,*derho;
  //Initializes gradients:
  gradcor=(double *) R_alloc(*nparcor, sizeof(double));
  derho=(double *) R_alloc(*npairs * *nparcor, sizeof(double));
  k=0;
    for(i=0;i<*npairs;i++){
    GradCorrFct(rho[i],cormod,eps,flagcor,gradcor,lags[i],0,parcor);
    for(m=0;m<*nparcor;m++){derho[k]=gradcor[m];k++;}}
  k=0;
  for(i=0;i<*nparcor;i++){
    for(m=0;m<*npairs;m++){
      drho[k]=derho[i+m* *nparcor];k++;}}
  return;
}

// compute the gradient matrix (numcoord*numtime...) for the spatial-temporal field:
void DCorrelationMat_st_tap(int *cormod,double *drho,double *eps,int *flagcor,
			int *nparcor,double *parcor,double *rho)
{
 int i=0,j=0,k=0,s=0;
 double *gradcor,*derho;
 gradcor=(double *) R_alloc(*nparcor, sizeof(double));
 derho=(double *) R_alloc(*npairs * nparcor[0], sizeof(double));

  for(i=0;i<*npairs;i++) {
     GradCorrFct(rho[i],cormod,eps,flagcor,gradcor,lags[i],lagt[i],parcor);
   for(s=0;s<*nparcor;s++){
     derho[k]=gradcor[s];
     k++;}
     }
 k=0;
 for(i=0;i<*nparcor;i++)
   for(j=0;j<*npairs;j++){
     drho[k]=derho[i+j* *nparcor];
     k++;}
 return;
}


// compute the gradient matrix (numcoord*numtime...) for the spatial-temporal field:
void DCorrelationMat_st(int *cormod,double *drho,double *eps,int *flagcor,
			int *nparcor,double *parcor,double *rho)
{
 int i=0,j=0,h=0,k=0,s=0,t=0,v=0,st=0,npa=0;
 double *gradcor,*derho;

 st=*ncoord * *ntime;
 npa=0.5*st*(st-1);
 gradcor=(double *) R_alloc(*nparcor, sizeof(double));
 derho=(double *) R_alloc(npa * nparcor[0], sizeof(double));

for(i=0;i<*ncoord;i++){
    for(t=0;t<*ntime;t++){
      for(j=i;j<*ncoord;j++){
	if(i==j){
	  for(v=t+1;v<*ntime;v++){
	        GradCorrFct(rho[h],cormod,eps,flagcor,gradcor,0,mlagt[t][v],parcor);
	        h++;
   for(s=0;s<*nparcor;s++){
     derho[k]=gradcor[s];
     k++;}}}
     else {
          for(v=0;v<*ntime;v++){
               GradCorrFct(rho[h],cormod,eps,flagcor,gradcor,mlags[i][j],mlagt[t][v],parcor);
               h++;
   for(s=0;s<*nparcor;s++){
     derho[k]=gradcor[s];
     k++;}}}}}}

 k=0;
 for(i=0;i<*nparcor;i++)
   for(j=0;j<npa;j++){
     drho[k]=derho[i+j* *nparcor];
     k++;}
 return;
}

// list of the gradients. Derivatives with respect ot the correlations parameters:
void GradCorrFct(double rho, int *cormod, double *eps, int *flag,
		 double *grad, double h, double u, double *par)
{
  int i=0;
  double power=0.0, power1=0.0, power2=0.0, power_s=0, power_t=0;
  double scale=0.0, scale_s=0, scale_t=0, smooth=0.0, sep=0;

  switch(*cormod)// Correlation functions are in alphabetical order
    {//spatial gradients of correlations:
    case 1:// Cauchy correlation function
      power2=par[0];
      scale=par[1];
      if(flag[0]==1){//power parameter
	  grad[i]=DCauchyPow(power2,scale,rho);i++;}
      if(flag[1]==1)//scale parameter
	grad[i]=DCauchySc(h, power2, scale, rho);
      break;
    case 2:// Exponential correlation function
      scale=par[0];//scale parameter
      if(flag[0] == 1) grad[i]=DExpoSc(h, scale, rho);
      break;
    case 3:// Gaussian correlation function
      scale=par[0];//scale parameter
      if(flag[0]==1) grad[i]=DGaussSc(h,scale,rho);
      break;
    case 4:// Generalised Cuachy correlation function
      power1=par[0];
      power2=par[1];
      scale=par[2];
      if(flag[0]==1){//power1 parameter
	grad[i]=DGenCauP1(h,power1,power2,scale,rho);i++;}
      if(flag[1]==1){//power2 parameter
	grad[i]=DGenCauP2(h,power1,scale,rho);i++;}
      if(flag[2]==1)//scale parameter
	grad[i]=DGenCauSc(h,power1,power2,scale,rho);
      break;
    case 5:// Sferical correlation function
      scale=par[0];//scale parameter
      if(flag[0]==1) grad[i]=DSferiSc(h,scale);
      break;
    case 6:// Stable correlation function
      power=par[0];
      scale=par[1];
      if(flag[0]==1){//power parameter
	grad[i]=DStabPow(h,power,scale,rho);i++;}
      //scale parameter
      if(flag[1]==1) grad[i]=DStabSc(h,power,scale,rho);
      break;
    case 7:// Whittle-Matern correlation function
      scale=par[0];
      smooth=par[1];
      if(flag[0]==1){//scale parameter
	grad[i]=DWhMatSc(eps,h,scale,smooth);i++;}
      //smooth parameter
      if(flag[1]==1) grad[i]=DWhMatSm(eps,h,scale,smooth);
      break;
      //spatio-temproal gradients of correlations:
      //separable correlations
    case 21://Gneiting spatio-temporal correlation
    case 26://Gneiting GC spatio-temporal correlation
      power_s=par[0];
      power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      sep=par[4];
      if(flag[0]==1){//spatial-power parameter
	grad[i]=DGneiting_pw_s(h,u,power_s,power_t,scale_s,scale_t,sep);i++;}
      if(flag[1]==1){//temporal-power parameter
	grad[i]=DGneiting_pw_t(h,u,power_s,power_t,scale_s,scale_t,sep);i++;}
      if(flag[2]==1){//spatial-scale parameter
	grad[i]=DGneiting_sc_s(h,u,power_s,power_t,scale_s,scale_t,sep);i++;}
      if(flag[3]==1){//temporal-scale parameter
	grad[i]=DGneiting_sc_t(h,u,power_s,power_t,scale_s,scale_t,sep);i++;}
      //separable parameter
      if(flag[4]==1) grad[i]=DGneiting_sep(h,u,power_s,power_t,scale_s,scale_t,sep);
      break;
    case 22://Decesare-Iaco spatio-temporal correlation
      power2=par[0];
      power_s=par[1];
      power_t=par[2];
      scale_s=par[3];
      scale_t=par[4];
      if(flag[0]==1){//power2 parameter
	grad[i]=DIaco_pw2(h,u,power_s,power_t,scale_s,scale_t,power2);i++;}
      if(flag[1]==1){//spatial-power parameter
	grad[i]=DIaco_pw_s(h,u,power_s,power_t,scale_s,scale_t,power2);i++;}
      if(flag[2]==1){//temporal-power parameter
	grad[i]=DIaco_pw_t(h,u,power_s,power_t,scale_s,scale_t,power2);i++;}
      if(flag[3]==1){//spatial-scale parameter
	grad[i]=DIaco_sc_s(h,u,power_s,power_t,scale_s,scale_t,power2);i++;}
      //temporal-scale parameter
      if(flag[4]==1) grad[i]=DIaco_sc_t(h,u,power_s,power_t,scale_s,scale_t,power2);
      break;
    case 23://Porcu spatio-temporal correlation
      power_s=par[0];
      power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      sep=par[4];
      if(flag[0]==1){//spatial-power parameter
	grad[i]=DPorcu_pw_s(h,u,power_s,power_t,scale_s,scale_t,sep);i++;}
      if(flag[1]==1){//temporal-power parameter
	grad[i]=DPorcu_pw_t(h,u,power_s,power_t,scale_s,scale_t,sep);i++;}
      if(flag[2]==1){//spatial-scale parameter
	grad[i]=DPorcu_sc_s(h,u,power_s,power_t,scale_s,scale_t,sep);i++;}
      if(flag[3]==1){//temporal-power parameter
	grad[i]=DPorcu_sc_t(h,u,power_s,power_t,scale_s,scale_t,sep);i++;}
      //separable parameter
      if(flag[4]==1) grad[i]=DPorcu_sep(h,u,power_s,power_t,scale_s,scale_t,sep);
      break;
    case 24://Stein spatio-temporal correlation???
      power_t=par[0];
      scale_s=par[1];
      scale_t=par[2];
      smooth=par[3];
      /*      if(flag[0]==1)
        {
	  grad[i]=DStein_pw_t(h,u,power_t,scale_s,scale_t,smooth);
	  i++;
        }
      if(flag[1]==1)
        {
	  grad[i]=DStein_sc_s(h,u,power_t,scale_s,scale_t,smooth);
	  i++;
        }
      if(flag[2]==1)
        {
	  grad[i]=DStein_sc_t(h,u,power_t,scale_s,scale_t,smooth);
	  i++;
        }
      if(flag[3]==1)
      grad[i]=DStein_sm(h,u,power_t,scale_s,scale_t,smooth);*/

      break;
      //Separable correlations
    case 41://Exponential-Cauchy
      power2=par[0];
      scale_s=par[1];
      scale_t=par[2];
      if(flag[0]==1){//power2 parameter
	grad[i]=DExp_Cauchy_pw2(h,u,power2,scale_s,scale_t);i++;}
      if(flag[1]==1){//spatial-scale parameter
	grad[i]=DExp_Cauchy_sc_s(h,u,power2,scale_s,scale_t);i++;}
      //temporal-scale parameter
      if(flag[2]==1) grad[i]=DExp_Cauchy_sc_t(h,u,power2,scale_s,scale_t);
      break;
    case 42://Double Exponential
      scale_s=par[0];
      scale_t=par[1];
      if(flag[0]==1){//spatial-scale parameter
	grad[i]=DStabSc(h,1,scale_s,rho);i++;}
      //temporal-scale parameter
      if(flag[1]==1) grad[i]=DStabSc(u,1,scale_t,rho);
      break;
    case 43://Exponential-Gaussian
      scale_s=par[0];
      scale_t=par[1];
      if(flag[0]==1){//spatial-scale parameter
	grad[i]=DStabSc(h,1,scale_s,rho);i++;}
      //temporal-scale parameter
      if(flag[1]==1) grad[i]=DStabSc(u,2,scale_t,rho);
      break;
    case 45://Matern-Cauchy
      power2=par[0];
      scale_s=par[1];
      scale_t=par[2];
      smooth=par[3];
      if(flag[0]==1){//power2 parameter
	grad[i]=DMat_Cauchy_pw2(h, u, power2, scale_s, scale_t, smooth);i++;}
      if(flag[1]==1){//spatial-scale parameter
	grad[i]=DMat_Cauchy_sc_s(h, u, power2, scale_s, scale_t, smooth);i++;}
      if(flag[2]==1){//temporal-scale parameter
	grad[i]=DMat_Cauchy_sc_t(h, u, power2, scale_s, scale_t, smooth);i++;}
      //smoothing parameter
      if(flag[3]==1) grad[i]=DMat_Cauchy_sm(h, u, eps, power2, scale_s, scale_t, smooth);
      break;
    case 46://Matern-Exponential
      scale_s=par[0];
      scale_t=par[1];
      smooth=par[2];
      if(flag[0]==1){//spatial-scale parameter
	grad[i]=DMat_Exp_sc_s(h, u, scale_s, scale_t, smooth);i++;}
      if(flag[1]==1){//temporal-scale parameter
	grad[i]=DMat_Exp_sc_t(h, u, scale_s, scale_t, smooth);i++;}
      //smoothing parameter
      if(flag[2]==1) grad[i]=DMat_Exp_sm(h, u, eps, scale_s, scale_t, smooth);
      break;
    case 47://Stable-Stable
      power_s=par[0];
      power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      if(flag[0]==1){//spatial-power parameter
	grad[i]=DStabPow(h,power_s,scale_s,rho);i++;}
      if(flag[1]==1){//temporal-power parameter
	grad[i]=DStabPow(u,power_t,scale_t,rho);i++;}
      if(flag[2]==1){//spatial-scale parameter
	grad[i]=DStabSc(h,power_s,scale_s,rho);i++;}
      if(flag[3]==1){//temporal-scale parameter
	grad[i]=DStabSc(u,power_t,scale_t,rho);}
      return;
    }
}


// list of spatial and spatial-temporal correlation functions:
void Comp_supp(double *c_supp,int *cormod, double h,double u, double *par)
{
  switch(*cormod) // Correlation functions are in alphabetical order
    {

    case 15:// Wendland1
    case 16:// Wendland2
    case 17:// Wendland3
      c_supp[0]=par[0];
      c_supp[1]=-LOW;
      break;
    case 100: // separable spacetime wendland
    case 101:
    case 102:
    case 103:
    case 104:
    case 105:
    case 106:
    case 107:
    case 108:
      c_supp[0]=par[0];
      c_supp[1]=par[1];
       break;
    case  109:   // quasi taper in time
      c_supp[0]=-LOW;
      c_supp[1]=par[1]/pow(1+pow(h/par[0],1),*tapsep);
      break;
    case  110: // quasi taper in space
      c_supp[1]=-LOW;
      c_supp[0]=par[0]/pow(1+pow(u/par[1],1),*tapsep);
      break;


/****************************************************************************************/
    }
}


// Computes the spatio-temporal variogram:
double Variogram(int *cormod, double h, double u, double *nuis, double *par)
{
  double vario=0.0;
  //Computes the variogram
  vario=nuis[1]+nuis[2]*(1-CorFct(cormod,h,u,par));
  return vario;
}
//***Define variograms for the Brown-Resnick process***
double VarioFct(int *cormod, double h, double *par, double u)
{
  double power_s=0.0, power_t=0.0, power1=0, power2=0;
  double scale_s=0.0, scale_t, vario=0.0;

  switch(*cormod) // Variogram functions are in alphabetical order
    {
    case 1:// Cauchy correlation function
      break;
    case 2:// Exponential correlation function
      scale_s=par[0];
      vario=VarioDobStable(h,1,1,scale_s,1,u);
      break;
    case 3:// Gaussian correlation function
      scale_s=par[0];
      vario=VarioDobStable(h,2,1,scale_s,1,u);
      break;
    case 4: // Generalised Cuachy correlation function
      power1=par[0];
      power2=par[1];
      scale_s=par[2];
      vario=VarioGCauchy(h,power1,power2,scale_s);
      break;
    case 6:// Stable correlation function
      power_s=par[0];
      scale_s=par[1];
      vario=VarioDobStable(h,power_s,1,scale_s,1,u);
      break;
    case 7://  Whittle-Matern correlation function
      break;
    case 47:
      power_s=par[0];
      power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      vario=VarioDobStable(h,power_s,power_t,scale_s,scale_t,u);
      break;
    }
  return vario;
}
// Variogram associated to the Stable class of correlation models:
double VarioDobStable(double lag, double power_s, double power_t, double scale_s, double scale_t, double tsep)
{
  double vario=0.0;
  // Computes the correlation:
  vario=pow(lag/scale_s,power_s)+pow(tsep/scale_t,power_t);
  return 2*vario;
}
// Variogram associated to the Stable class of correlation models:
double VarioGCauchy(double lag, double power1, double power2, double scale)
{
  double vario=0.0;
  // Computes the correlation:
  vario=pow(sqrt(power2/power1)*lag/scale,power1);
  return 2*vario;
}
// list of the gradients. Derivatives with respect ot the variograms parameters:
void GradVarioFct(double vario, int *cormod, double *eps, int *flag,
		  double *grad, double lag, double *par, double tsep)
{
  int i=0;
  double power=0.0, power1=0, power2=0, power_s=0, power_t=0;
  double scale=0.0, scale_s=0, scale_t=0;

  switch(*cormod)// Variogram functions are in alphabetical order
    {
    case 1:// Cauchy variogram function
      break;
    case 2:// Exponential correlation function
      scale=par[0];
      if(flag[0] == 1)
	grad[i]=-vario/scale;
      break;
    case 3:// Gaussian variogram function
      scale=par[0];
      if(flag[0]==1)
	grad[i]=-2*vario/scale;
      break;
    case 4:// Generalised Cuachy variogram function
      power1=par[0];
      power2=par[1];
      scale=par[2];
      if(flag[0]==1)
	grad[i]=-2*sqrt(power2/power1)*vario/scale;
       break;
    case 6:// Stable variogram function
      power=par[0];
      scale=par[1];
      if(flag[0]==1){//power parameter
	  grad[i]=vario*log(lag/scale);i++;}
      //scale parameter
      if(flag[1]==1) grad[i]=-vario*power/scale;
      break;
    case 7:// Whittle-Matern variogram function
      break;
    case 47:// Stable-Stable variogram function
      power_s=par[0];
      power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      if(flag[0]==1){//spatial-power parameter
	grad[i]=pow(lag/scale_s,power_s)*log(lag/scale_s);i++;}
      if(flag[1]==1){//temporal-power parameter
	grad[i]=pow(tsep/scale_t,power_t)*log(lag/scale_t);i++;}
      if(flag[2]==1){//spatial-scale parameter
	grad[i]=-pow(lag/scale_s,power_s)*power_s/scale_s;i++;}
      if(flag[3]==1){//temporal-scale parameter
	grad[i]=-pow(tsep/scale_t,power_t)*power_t/scale_t;i++;}
      break;
    }
  return;
}
// Vector of spatio-temporal correlations:
void VectCorrelation(double *rho, int *cormod, double *h, int *nlags, int *nlagt, double *par, double *u)
{
  int i,j,t=0;
  for(j=0;j<*nlagt;j++)
    for(i=0;i<*nlags;i++){
      rho[t]=CorFct(cormod, h[i], u[j], par);
      t++;}
  return;
}
// Vector of spatial extremal coefficients:
void ExtCoeff(int *cormod, double *extc, double *lag, int *model,
	      int *nlags, double *nuis, double *par)
{
  int i;
  double df1=nuis[0]+1, rho=0;
  switch(*model){// Call to the extremal coefficients:
  case 3:// Brown-Resnick process
    for(i=0;i<*nlags;i++)
      extc[i]=2*pnorm(0.5*sqrt(VarioFct(cormod,lag[i],par,0)),0,1,1,0);//compute the extremal coeff.
    break;
  case 4:// Extremal Gaussian process
    for(i=0;i<*nlags;i++)
      extc[i]=1+sqrt(0.5*(1-nuis[0]*CorFct(cormod,lag[i],0,par)));//compute the extremal-g coeff.
    break;
  case 5:// Extremal T process:
    for(i=0;i<*nlags;i++){
      rho=nuis[1]*CorFct(cormod,lag[i],0,par);//compute the correlation
      extc[i]=2*pt(sqrt(df1*(1-rho)/(1+rho)),df1,1,0);}//compute  the extremal-t coeff.
    break;}
  return;
}
