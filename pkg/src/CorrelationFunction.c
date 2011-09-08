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
      if(scale_s<=0 || scale_t<=0 || power_s<0 || power_s>2 || power_t<0 || power_t>2) rho=-2;
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
      // END separable correlation functions:
    }
  return rho;
}
// list of spatial and spatial-temporal correlation functions:
double CorFct(int *cormod, double h, double u, double *par)
{
  double arg=0.0, power=0.0, power1=0.0, power2=0.0, power_s=0.0, power_t=0.0;
  double rho=0.0, sep=0, scale=0.0, smooth=0.0, scale_s=0.0, scale_t=0;

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
    case 15:// Wendland1
      rho=CorFunWend1(h);
      break;
    case 16:// Wendland2
      rho=CorFunWend2(h);
      break;
    case 17:// Wendland3
      rho=CorFunWend3(h);
      break;
      // START non-separable correlation functions:
    case 21:   //Gneiting correlation model as in (14) Gneitint (2002) with tau=1
      power_s=par[0];
      power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      sep=par[4];
      arg=1+pow(u, power_t)/scale_t;
      rho=exp(-pow(h, power_s)/(scale_s*pow(arg, 0.5*sep*power_s)))/arg;
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
      if(sep) rho=pow(0.5*pow(1+pow(h, power_s)/scale_s,sep)+0.5*pow(1+pow(u, power_t)/scale_t,sep),-1/sep);
      else rho=pow((1+pow(h, power_s)/scale_s)*(1+pow(u,power_t)/scale_t),-1);
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
      // END non-separable correlation functions
      // START separable correlation functions:
    case 41:// Exp-Cauchy:
      power2=par[0];
      scale_s=par[1];
      scale_t=par[2];
      rho=exp(-h/scale_s)*pow((1+pow(u/scale_t, 2)), -power2);
      break;
    case 42:// Double exp:
      scale_s=par[0];
      scale_t=par[1];
      rho=exp(-h/scale_s-u/scale_t);
      break;
    case 43:// Exp-Gauss:
      scale_s=par[0];
      scale_t=par[1];
      rho=exp(-h/scale_s-pow(u/scale_t,2));
      break;
    case 45:// Matern-Cauchy:
      power2=par[0];
      scale_s=par[1];
      scale_t=par[2];
      smooth=par[3];
      if(h==0) arg=1;
      else arg=pow(2, 1-smooth)/gammafn(smooth)*pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);
      rho=arg*pow((1+pow(u/scale_t, 2)),-power2);
      break;
    case 46:// Matern-exp:
      scale_s=par[0];
      scale_t=par[1];
      smooth=par[2];
      if(h==0) arg=1;
      else arg=pow(2, 1-smooth)/gammafn(smooth)*pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);
      rho=arg*exp(-u/scale_t);
      break;
      // END separable correlation functions:
    }
  return rho;
}
// Cauhcy class of correlation models:
double CorFunCauchy(double lag, double power2, double scale)
{
  double rho=0.0;
  // Computes the correlation:
  rho=pow((1+pow(lag/scale,2)),- power2);
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
  rho=pow(2,1-smooth)/gammafn(smooth)*pow(lag/scale,smooth) *
	bessel_k(lag/scale,smooth,1);
  return rho;
}
double CorFunWend1(double lag)
{
  double rho=0.0,x=0;
  x=lag*pow(*maxdist,-1);
  if(x<=1) rho=pow(1-x,2)*(1+0.5*x);
  else rho=0;
  return rho;
}
double CorFunWend2(double lag)
{
  double rho=0.0,x=0;
   x=lag*pow(*maxdist,-1);
  if(x<=1) rho=pow(1-x,4)*(1+4*x);
  else rho=0;
  return rho;
}
double CorFunWend3(double lag)
{
  double rho=0.0,x=0;
   x=lag*pow(*maxdist,-1);
  if(x<=1) rho=pow(1-x,6)*(1+6*x+35*pow(x,2)/3);
  else rho=0;
  return rho;
}
// Computation of the upper (lower) triangular spatial correlation matrix:
void CorrelationMat(double *rho, int *cormod, double *nuis, double *par)
{
  int i;// check the paramaters range:
  if(nuis[1]<0 || nuis[2]<=0 || CheckCor(cormod,par)==-2){
    rho[0]=-2;
    return;}// compute the correlations:
  for(i=0;i<*npairs;i++)
    rho[i]=CorFct(cormod,lags[i],0,par);
  return;
}
// Computation of the upper (lower) triangular spatial-temporal correlation matrix:
void CorrelationMat_st(double *rho, int *cormod, double *nuis, double *par)
{
  int i=0, j=0, t=0, v=0, p=0;
  if(nuis[1]<0 || nuis[2]<=0 || CheckCor(cormod,par)==-2){
    rho[0]=-2;
    return;}
  // first spatial loop:
  for(i=0;i<*ncoord;i++)
    {// first temporal loop:
      for(t=0;t<*ntime;t++)
	{// second spatial loop:
	  for(j=i;j<*ncoord;j++)
	    {
	      if(i==j)
		for(v=(t+1);v<*ntime;v++)// second temporal loop:
		  {//marignal temporal correlations
		  rho[p]=CorFct(cormod, 0, mlagt[t][v], par);p++;}
	      else
		for(v=0;v<*ntime;v++){// second temporal loop:
		  rho[p]=CorFct(cormod, mlags[i][j], mlagt[t][v], par);//spatial-temporal correlations
		  p++;}
	    }
	}
    }
  return;
}
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
  return power2*rho/power1*(log(1+pow(lag/scale,power1))/power1-
			    pow(lag/scale,power1)*log(lag/scale)/
			    (1+pow(lag/scale,power1)));
}
// Derivatives with respect to power2 of the generalised Cauchy correlation model:
double DGenCauP2(double lag, double power1, double scale, double rho)
{
  return -rho*log(1+pow(lag/scale,power1))/power1;
}
// Derivatives with respect to scale of the generalised Cauchy correlation model:
double DGenCauSc(double lag, double power1, double power2, double scale, double rho)
{
  return rho/(1+pow(lag/scale,2))*power2*pow(lag,power1)/pow(scale,power1+1);
}
// Derivatives with respect to scale of the sferical correlation model:
double DSferiSc(double lag, double scale)
{
  if(lag<=scale)
    return 1.5*lag*(1-pow(lag/scale, 2))/pow(scale, 2);
  else
    return 0.0;
}
// Derivatives with respect to power of the Stable correlation model:
double DStabPow(double lag, double power, double scale, double rho)
{
  return -rho*pow(lag/scale,power)*log(lag/scale);
}
// Derivatives with respect to scale of the Stable correlation model:
double DStabSc(double lag, double power, double scale, double rho)
{
  return rho*pow(lag/scale,power-1)*power*lag/pow(scale,2);
}
// Derivatives with respect to scale of the Whittle-Matern correlation model:
double DWhMatSc(double *eps, double lag, double scale, double smooth)
{
  double pscale=0.0;
  pscale=(bessel_k(lag/(scale+*eps),smooth,1)-
	  bessel_k(lag/scale,smooth,1))/ *eps;
  return pow(2,1-smooth)/gammafn(smooth)*pow(lag/scale,smooth)*
    (pscale-smooth*bessel_k(lag/scale,smooth,1)/scale);
}
// Derivatives with respect to smooth of the Whittle-Matern correlation model:
double DWhMatSm(double *eps, double lag, double scale, double smooth)
{
  double psmooth=0.0;
  psmooth=(bessel_k(lag/scale,smooth+ *eps,1)-
	   bessel_k(lag/scale,smooth,1))/ *eps;
  return pow(2,1-smooth)*pow(lag/scale,smooth)/
    gammafn(smooth)*(log(lag/scale)-log(2)-
		     digamma(smooth)*bessel_k(lag/scale,smooth,1)+psmooth);
}
// Derivatives with respect to spatial scale of the Gneiting correlation model:
double DGneiting_sc_s(double h,double u,double power_s,double power_t,
		      double scale_s,double scale_t,double sep)
{
  double arg=0.0,rho=0.0;
  arg=1+pow(u, power_t)/scale_t;
  rho=exp(-pow(h, power_s)/(scale_s*pow(arg, 0.5*sep*power_s)));
  return pow(h, power_s)*rho/(pow(arg, 0.5*sep*power_s+1)*pow(scale_s,2));
}
// Derivatives with respect to temporal scale of the Gneiting correlation model:
double DGneiting_sc_t(double h,double u, double power_s,double power_t,
		      double scale_s,double scale_t,double sep)
{
  double arg=0.0,rho=0.0,a=0.0,b=0.0;
  arg=1+pow(u, power_t)/scale_t;
  rho=exp(-pow(h, power_s)/(scale_s*pow(arg, 0.5*sep*power_s)));
  a=(pow(u, power_t)*rho)/pow(arg*scale_t,2);
  b=-(0.5*pow(h, power_s)*sep*power_s*pow(u, power_t)*rho)/(pow(arg,2+0.5*sep*power_s)*scale_s*pow(scale_t,2));
  return b+a;
}
// Derivatives with respect to spatial power of the Gneiting correlation model:
double DGneiting_pw_s(double h,double u, double power_s,double power_t,
		      double scale_s,double scale_t,double sep)
{
  double arg=0.0,rho=0.0,a=0.0,b=0.0;
  arg=1+pow(u, power_t)/scale_t;
  rho=exp(-pow(h, power_s)/(scale_s*pow(arg, 0.5*sep*power_s)));
  if(h && arg){
    a=-(pow(h, power_s)*log(h))/(scale_s*pow(arg,0.5*sep*power_s));
    b=(0.5*pow(h, power_s)*sep*log(arg))/(scale_s*pow(arg, 0.5*sep*power_s));}
  return (a+b)*(rho/arg);
}
// Derivatives with respect to temporal power of the Gneiting correlation model:
double DGneiting_pw_t(double h,double u, double power_s,double power_t,
		      double scale_s,double scale_t,double sep)
{
  double arg=0.0,rho=0.0,a=0.0,b=0.0;
  arg=1+pow(u, power_t)/scale_t;
  rho=exp(-pow(h, power_s)/(scale_s*pow(arg, 0.5*sep*power_s)));
  if(u){
    a=(0.5*pow(h, power_s)*sep*power_s*pow(u, power_t)*log(u)*rho)/(pow(arg,2+0.5*sep*power_s)*scale_s*scale_t);
    b=-(rho*pow(u, power_t)*log(u))/(pow(arg,2)*scale_t);}
  return a+b;
}
// Derivatives with respect to the separable parameter of the Gneiting correlation model:
double DGneiting_sep(double h,double u, double power_s,double power_t,
		     double scale_s,double scale_t,double sep)
{
  double arg=0,rho=0,a=0;
  arg=1+pow(u, power_t)/scale_t;
  rho=exp(-pow(h, power_s)/(scale_s*pow(arg, 0.5*sep*power_s)));
  if(arg)   a=(0.5*pow(h, power_s)*power_s*log(arg)*rho)/(pow(arg, 0.5*sep*power_s+1)*scale_s);
  return a;
}
// Derivatives with respect to the separable parameter of the Gneiting correlation model:
double DIaco_sc_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double power2)
{
 double rho=0,arg=0;
 arg=1+pow(h/scale_s, power_s)+pow(u/scale_t, power_t);
 rho=pow(arg,-power2);
 return (rho*power_s*power2*pow(h/scale_s, power_s))/(arg*scale_s);
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
  else  arg=1;
  if(h)  arg2=-arg3*arg*(log(2)+digamma(smooth)
               -log(h/scale_s)-psmooth/bessel_k(h/scale_s, smooth, 1));
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

 if(h) {arg1=(pow(2, 1-smooth)/gammafn(smooth))*pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);
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
  arg2=-exp(-u/scale_t)*arg*(log(2)+digamma(smooth)
               -log(h/scale_s)-psmooth/bessel_k(h/scale_s, smooth, 1));}
 else arg2=0;
 return arg2;
}

double DPorcu_sc_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep)
{
  double arg=0,arg1=0,arg2=0,rho=0;
  arg1=1+pow(h, power_s)/scale_s;
  arg2=1+pow(u, power_t)/scale_t;
  arg=pow(0.5*pow(arg1, sep)+0.5*pow(arg2,sep),-1/sep);
  if(sep) rho=pow(arg,-1/sep);
  else rho=pow(arg1*arg2,-1);
  return (0.5*pow(h, power_s)*rho*arg1)/(arg*arg1*pow(scale_s,2));
}
double DPorcu_sc_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep)
{
  double arg=0,arg1=0,arg2=0,rho=0;
  arg1=1+pow(h, power_s)/scale_s;
  arg2=1+pow(u, power_t)/scale_t;
  arg=pow(0.5*pow(arg1, sep)+0.5*pow(arg2,sep),-1/sep);
  if(sep) rho=pow(arg,-1/sep);
  else rho=pow(arg1*arg2,-1);
  return (0.5*pow(u, power_t)*rho*arg2)/(arg*arg2*pow(scale_t,2));
}
double DPorcu_pw_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep)
{
  double a=0,arg=0,arg1=0,arg2=0,rho=0;
  arg1=1+pow(h, power_s)/scale_s;
  arg2=1+pow(u, power_t)/scale_t;
  arg=pow(0.5*pow(arg1, sep)+0.5*pow(arg2,sep),-1/sep);
  if(sep) rho=pow(arg,-1/sep);
  else rho=pow(arg1*arg2,-1);
  if(h) a=-(0.5*pow(h, power_s)*rho*arg1*log(h))/(arg*arg1*scale_s);
  return a;
}
double DPorcu_pw_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep)
{
  double a=0,arg=0,arg1=0,arg2=0,rho=0;
  arg1=1+pow(h, power_s)/scale_s;
  arg2=1+pow(u, power_t)/scale_t;
  arg=pow(0.5*pow(arg1, sep)+0.5*pow(arg2,sep),-1/sep);
  if(sep) rho=pow(arg,-1/sep);
  else rho=pow(arg1*arg2,-1);
  if(h) a=-(0.5*pow(u, power_t)*rho*arg2*log(u))/(arg*arg2*scale_t);
  return a;
}


double DPorcu_sep(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep)
{
  double arg=0,arg1=0,arg2=0,rho=0,a=0;
  arg1=1+pow(h, power_s)/scale_s;
  arg2=1+pow(u, power_t)/scale_t;
  arg=pow(0.5*pow(arg1, sep)+0.5*pow(arg2,sep),-1/sep);
  if(sep) rho=pow(arg,-1/sep);
  else rho=pow(arg1*arg2,-1);
  if(arg1&&arg2)
    a=rho*(log(arg)/pow(sep,2)-((0.5*log(arg1)*pow(arg1,sep)+0.5*log(arg2)*pow(arg2,sep))))/(sep*arg);
  return a;
}

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
}

// list of the gradients. Derivatives with respect ot the correlations parameters:
void GradCorrFct(double rho, int *cormod, double *eps, int *flag,
		 double *grad, double h, double u, double *par)
{
  int i=0;
  double power=0.0, power1=0.0, power2=0.0, power_s=0, power_t=0;
  double scale=0.0, scale_s=0, scale_t=0, smooth=0.0, sep=0;

  switch(*cormod)// Correlation functions are in alphabetical order
    {
    case 1:// Cauchy correlation function
      power2=par[0];
      scale=par[1];
      if(flag[0]==1)
	{
	  grad[i]=DCauchyPow(power2, scale, rho);
	  i++;
	}
      if(flag[1]==1)
	grad[i]=DCauchySc(h, power2, scale, rho);
      break;
    case 2:// Exponential correlation function
      scale=par[0];
      if(flag[0] == 1)
	grad[i]=DExpoSc(h, scale, rho);
      break;
    case 3:// Gaussian correlation function
      scale=par[0];
      if(flag[0]==1)
	grad[i]=DGaussSc(h,scale,rho);
      break;
    case 4:// Generalised Cuachy correlation function
      power1=par[0];
      power2=par[1];
      scale=par[2];
      if(flag[0]==1)
	{
          grad[i]=DGenCauP1(h,power1,power2,scale,rho);
	  i++;
	}
      if(flag[1]==1)
	{
	  grad[i]=DGenCauP2(h,power1,scale,rho);
	  i++;
	}
      if(flag[2]==1)
	grad[i]=DGenCauSc(h,power1,power2,scale,rho);
      break;
    case 5:// Sferical correlation function
      scale=par[0];
      if(flag[0]==1)
	grad[i]=DSferiSc(h,scale);
      break;
    case 6:// Stable correlation function
      power=par[0];
      scale=par[1];
      if(flag[0]==1)
	{
	  grad[i]=DStabPow(h,power,scale,rho);
	  i++;
	}
      if(flag[1]==1)
	grad[i]=DStabSc(h,power,scale,rho);
      break;
    case 7:// Whittle-Matern correlation function
      scale=par[0];
      smooth=par[1];
      if(flag[0]==1)
	{
	  grad[i]=DWhMatSc(eps,h,scale,smooth);
	  i++;
	}
      if(flag[1]==1)
	grad[i]=DWhMatSm(eps,h,scale,smooth);
      break;
      //  space-time gradient correlation function:
    case 21:
      power_s=par[0];
      power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      sep=par[4];
      if(flag[0]==1)
        {
	  grad[i]=DGneiting_pw_s(h,u,power_s,power_t,scale_s,scale_t,sep);
	  i++;
        }
      if(flag[1]==1)
        {
	  grad[i]=DGneiting_pw_t(h,u,power_s,power_t,scale_s,scale_t,sep);
	  i++;
        }
      if(flag[2]==1)
        {
	  grad[i]=DGneiting_sc_s(h,u,power_s,power_t,scale_s,scale_t,sep);
	  i++;
        }
      if(flag[3]==1)
        {
	  grad[i]=DGneiting_sc_t(h,u,power_s,power_t,scale_s,scale_t,sep);
	  i++;
        }
      if(flag[4]==1)
	grad[i]=DGneiting_sep(h,u,power_s,power_t,scale_s,scale_t,sep);
      break;
    case 22:
      power2=par[0];
      power_s=par[1];
      power_t=par[2];
      scale_s=par[3];
      scale_t=par[4];
      if(flag[0]==1)
        {
	  grad[i]=DIaco_pw2(h,u,power_s,power_t,scale_s,scale_t,power2);
	  i++;
        }
      if(flag[1]==1)
        {
	  grad[i]=DIaco_pw_s(h,u,power_s,power_t,scale_s,scale_t,power2);
	  i++;
        }
      if(flag[2]==1)
        {
	  grad[i]=DIaco_pw_t(h,u,power_s,power_t,scale_s,scale_t,power2);
	  i++;
        }
      if(flag[3]==1)
        {
	  grad[i]=DIaco_sc_s(h,u,power_s,power_t,scale_s,scale_t,power2);
	  i++;
        }
      if(flag[4]==1)
        grad[i]=DIaco_sc_t(h,u,power_s,power_t,scale_s,scale_t,power2);
      break;

    case 23:
      power_s=par[0];
      power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      sep=par[4];
      if(flag[0]==1)
        {
	  grad[i]=DPorcu_pw_s(h,u,power_s,power_t,scale_s,scale_t,sep);
	  i++;
        }
      if(flag[1]==1)
        {
	  grad[i]=DPorcu_pw_t(h,u,power_s,power_t,scale_s,scale_t,sep);
	  i++;
        }
      if(flag[2]==1)
        {
	  grad[i]=DPorcu_sc_s(h,u,power_s,power_t,scale_s,scale_t,sep);
	  i++;
        }
      if(flag[3]==1)
        {
	  grad[i]=DPorcu_sc_t(h,u,power_s,power_t,scale_s,scale_t,sep);
	  i++;
        }
      if(flag[4]==1)
        grad[i]=DPorcu_sep(h,u,power_s,power_t,scale_s,scale_t,sep);
      break;

    case 24:
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

    case 41:
      power2=par[0];
      scale_s=par[1];
      scale_t=par[2];
      if(flag[0]==1)
        {
	  grad[i]=DExp_Cauchy_pw2(h,u,power2,scale_s,scale_t);
	  i++;
        }
      if(flag[1]==1)
        {
	  grad[i]=DExp_Cauchy_sc_s(h,u,power2,scale_s,scale_t);
	  i++;
        }
      if(flag[2]==1)
        grad[i]=DExp_Cauchy_sc_t(h,u,power2,scale_s,scale_t);
      break;

    case 42:
      scale_s=par[0];
      scale_t=par[1];

      if(flag[0]==1)
        {
	  grad[i]=DExp_Exp_sc_s(h, u,  scale_s, scale_t);
	  i++;
        }
      if(flag[1]==1)
        grad[i]=DExp_Exp_sc_t(h, u,  scale_s, scale_t);

      break;

    case 43:
      scale_s=par[0];
      scale_t=par[1];
      if(flag[0]==1)
        {
	  grad[i]=DExp_Gauss_sc_s(h, u, scale_s,scale_t);
	  i++;
        }
      if(flag[1]==1)
        grad[i]=DExp_Gauss_sc_t(h, u, scale_s, scale_t);
      break;

    case 45:
      power2=par[0];
      scale_s=par[1];
      scale_t=par[2];
      smooth=par[3];
      if(flag[0]==1)
        {
	  grad[i]=DMat_Cauchy_pw2(h, u, power2, scale_s, scale_t, smooth);
	  i++;
        }
      if(flag[1]==1)
        {
	  grad[i]=DMat_Cauchy_sc_s(h, u, power2, scale_s, scale_t, smooth);
	  i++;
        }
      if(flag[2]==1)
        {
	  grad[i]=DMat_Cauchy_sc_t(h, u, power2, scale_s, scale_t, smooth);
	  i++;
        }
      if(flag[3]==1)
        grad[i]=DMat_Cauchy_sm(h, u, eps, power2, scale_s, scale_t, smooth);

      break;

    case 46:
      scale_s=par[0];
      scale_t=par[1];
      smooth=par[2];
      if(flag[0]==1)
        {
	  grad[i]=DMat_Exp_sc_s(h, u, scale_s, scale_t, smooth);
	  i++;
        }
      if(flag[1]==1)
        {
	  grad[i]=DMat_Exp_sc_t(h, u, scale_s, scale_t, smooth);
	  i++;
        }
      if(flag[2]==1)
	grad[i]=DMat_Exp_sm(h, u, eps, scale_s, scale_t, smooth);

      break;
    }
  return;
}


// list of the gradients. Derivatives with respect ot the correlations parameters:
/*void GradCorrFct(double rho, int *cormod, double *eps, int *flag,
		 double *grad, double lag, double *par)
{
  int i=0;
  double power=0.0, power1=0.0, power2=0.0, scale=0.0, smooth=0.0;
  double parscale=0.0, parsmooth=0.0;

  switch(*cormod)// Correlation functions are in alphabetical order
    {
    case 1:// Cauchy correlation function
      power2=par[0];
      scale=par[1];
      if(flag[0]==1)
	{
	  grad[i]=DCauchyPow(power2, scale, rho);
	  i++;
	}
      if(flag[1]==1)
	grad[i]=DCauchySc(lag, power2, scale, rho);
      break;
    case 2:// Exponential correlation function
      scale=par[0];
      if(flag[0] == 1)
	grad[i]=DExpoSc(lag, scale, rho);
      break;
    case 3:// Gaussian correlation function
      scale=par[0];
      if(flag[0]==1)
	grad[i]=DGaussSc(lag,scale,rho);
      break;
    case 4:// Generalised Cuachy correlation function
      power1=par[0];
      power2=par[1];
      scale=par[2];
      if(flag[0]==1)
	{
          grad[i]=DGenCauP1(lag,power1,power2,scale,rho);
	  i++;
	}
      if(flag[1]==1)
	{
	  grad[i]=DGenCauP2(lag,power1,scale,rho);
	  i++;
	}
      if(flag[2]==1)
	grad[i]=DGenCauSc(lag,power1,power2,scale,rho);
      break;
    case 5:// Sferical correlation function
      scale=par[0];
      if(flag[0]==1)
	grad[i]=DSferiSc(lag,scale);
      break;
    case 6:// Stable correlation function
      power=par[0];
      scale=par[1];
      if(flag[0]==1)
	{
	  grad[i]=DStabPow(lag,power,scale,rho);
	  i++;
	}
      if(flag[1]==1)
	grad[i]=DStabSc(lag,power,scale,rho);
      break;
    case 7:// Whittle-Matern correlation function
      scale=par[0];
      smooth=par[1];
      if(flag[0]==1)
	{
	  grad[i]=DWhMatSc(eps,lag,scale,smooth);
	  i++;
	}
      if(flag[1]==1)
	grad[i]=DWhMatSm(eps,lag,scale,smooth);
      break;
    }
  return;
  }*/
// Computes the spatial-temporal variogram:
double Variogram(int *cormod, double h, double u, double *nuis, double *par)
{
  double vario=0.0;
  //Computes the variogram
  vario=nuis[1]+nuis[2]*(1-CorFct(cormod,h,u,par));
  return vario;
}
double VarioFct(int *cormod, double lag, double *par)
{
  double power=0.0, scale=0.0, vario=0.0;

  switch(*cormod) // Variogram functions are in alphabetical order
    {
    case 1:// Cauchy correlation function
      break;
    case 2:// Exponential correlation function
      power=1;
      scale=par[0];
      vario=VarioStable(lag, power, scale);
      break;
    case 3:// Gaussian correlation function
      power=2;
      scale=par[0];
      vario=VarioStable(lag, power, scale);
      break;
    case 4: // Generalised Cuachy correlation function
      break;
    case 6:// Stable correlation function
      power=par[0];
      scale=par[1];
      vario=VarioStable(lag, power, scale);
      break;
    case 7://  Whittle-Matern correlation function
      break;
    }
  return vario;
}
// Variogram associated to the Stable class of correlation models:
double VarioStable(double lag, double power, double scale)
{
  double vario=0.0;
  // Computes the correlation:
  vario=pow(lag/scale,power);
  return 2*vario;
}

// list of the gradients. Derivatives with respect ot the variograms parameters:
void GradVarioFct(double vario, int *cormod, double *eps, int *flag,
		  double *grad, double lag, double *par)
{
  int i=0;
  double power=0.0, scale=0.0;

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
       break;
    case 6:// Stable variogram function
      power=par[0];
      scale=par[1];
      if(flag[0]==1)
	{
	  grad[i]=vario*log(lag/scale);
	  i++;
	}
      if(flag[1]==1)
	grad[i]=-vario*power/scale;
      break;
    case 7:// Whittle-Matern variogram function
      break;
    }
  return;
}

void VectCorrelation(double *rho, int *cormod, double *lag, int *nlags, double *par)
{
  int i;

  for(i=0;i<*nlags;i++)
    rho[i]=CorFct(cormod, lag[i], 0, par);
  return;
}

void ExtCoeff(int *cormod, double *extc, double *lag, int *model,
	      int *nlags, double *nuis, double *par)
{
  int i;
  double df1=nuis[0]+1, rho=0;
  switch(*model){// Call to the extremal coefficients:
  case 3:// Brown-Resnick process
    for(i=0;i<*nlags;i++)
      extc[i]=2*pnorm(0.5*sqrt(VarioFct(cormod,lag[i],par)),0,1,1,0);//compute the extremal coeff.
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
