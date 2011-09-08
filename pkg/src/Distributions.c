#include "header.h"

// Procedures are in alphabetical order.
// Student-t probability density function:
double dts(double x, double df)
{
  double res=0.0, df1=df+1;
  res=gammafn(0.5*df1)/sqrt(df*M_PI)/gammafn(0.5*df)*
    pow(1+pow(x,2)/df, - 0.5*df1);
  return res;
}
// First order derivative of the Student-t density with respect to the argument:
double d1x_dt(double x, double df)
{
  double res=0.0;
  res=-dt(x,df,0)*(df+1)*x/
    (1+pow(x,2)/df)/df;
  return res;
}
// Second order derivative of the Student-t density with respect to the argument:
double d2x_dt(double x, double df)
  {
    double res=0.0, y=0.0;
    y=1+pow(x,2)/df;
    res=-d1x_dt(x,df)*(x*(df+3)/y/df-1/x);
    return res;
  }
// First order derivative of the Student-t distribution
// with respect to the degree of freedom:
double ddf_pt(double x, double df)
  {
    double epsabs=0.0, epsrel=0.0, integ=0.0, tinteg=0.0;
    double abserr=0.0, origin=0.0, q=0.0, res=0.0, *work;
    int inf=0.0, neval=0.0, ier=0.0, limit=0.0, lenw=0.0;
    int last=0.0, *iwork;
    // Initialize the parameters and variable for the integrate fun:
    inf=-1;
    epsabs=1e-5;
    epsrel=1e-5;
    limit=100;
    lenw=4*limit;
    iwork=(int *) R_alloc(limit,sizeof(int));
    work=(double *) R_alloc(lenw,sizeof(double));
    // Checks the sign of the argument:
    if(x <= 0)
      Rdqagi(integr_pt,(void*)&df,&x,&inf,&epsabs,&epsrel,
	     &integ,&abserr,&neval,&ier,&limit,&lenw,&last,iwork,work);
    else
      {
	q=-x;
	Rdqagi(integr_pt,(void*)&df,&origin,&inf,&epsabs,&epsrel,
	       &tinteg,&abserr,&neval,&ier,&limit,&lenw,&last,iwork,work);
	Rdqagi(integr_pt,(void*)&df,&q,&inf,&epsabs,&epsrel,&integ,&abserr,
	       &neval,&ier,&limit,&lenw,&last,iwork,work);
	integ=2*tinteg-integ;
	}
    res=0.5*pt(x,df,1,0)*(digamma(0.5*(df+1))-digamma(0.5*df)-1/df)+integ;
    return res;
  }
// Integrand function (derivatives of the Student-t cdf):
double int_pt(double x, double df)
  {
    double res=0.0, x2=0.0, y=0.0;
    x2=pow(x,2);
    y=1+x2/df;
    res=0.5*dt(x,df,0)*((df+1)*x2/pow(df,2)/y-log(y));
    return res;
  }
// Vectorised integrand function:
void integr_pt(double *x, int n, void *ex)
{
  int i=0;
  double d=0.0;
  
  //*df=*((double *)ex);
  //d=*df;
  d=*((double *)ex);
  for(i=0;i<n;i++)
    x[i]=int_pt(x[i],d);
  return;
}
// First order derivative of the Student-t probability density function
// with respect to the degree of freedom with argument g(x):
double ddf_t_dt(double x, double clc, double df, double somc2)
  {
    double df1=df+1, res=0.0, nu=df-1, y=0.0;
    y=1+pow(x,2)/df;
    res=0.5*dt(x,df,0)*(digamma(0.5*df1)-digamma(0.5*df)-
			log(y)-1/df+2*df1*x*clc/nu/sqrt(df)/somc2/y);
    return res;
  }
// First order derivative of the first order derivative of the 
// Student-t probability density function (with respect to the arg)
// with respect to the degree of freedom with argument g(x):
double ddf_t_d1x_dt(double x, double clc, double df, double somc2)
  {
    double df1=df+1, nu=df-1, res=0.0, sdf=sqrt(df), y=0.0;
    y=1+pow(x,2)/df;    
    res=0.5*d1x_dt(x,df)*(digamma(0.5*df1)-digamma(0.5*df)-
			  log(y)-2/df+2*x*clc*(df1+2)/nu/sdf/somc2/y+
			  2/df1-2*clc*sdf/nu/x/somc2);
    return res;
  }
// Switch from a GEV distribution to another GEV distribution:
void Dist2Dist(double *data, double *eloc, double *escale, double *eshape, 
	       int *ndata, int *nsite, double *ploc, double *pscale, 
	       double *pshape, int *type, double *res)
{
  int i=0, k=0;

  switch(*type)
    {
    case 0: // Transform from GEV to UNIFORM
      for(k = 0; k < *nsite; k++)
	for(i = 0; i < *ndata; i++)
	  res[i + k * *ndata] = pgev(data[i + k * *ndata], eloc[k], escale[k], eshape[k]);
      break;
    case 1: // Transform from GEV to unit FRECHET
      for(k = 0; k < *nsite; k++)
	for(i = 0; i < *ndata; i++)
	  res[i + k * *ndata] = qgev(pgev(data[i + k * *ndata], eloc[k], escale[k], eshape[k]), 1, 1, 1);
      break;
    case 2: // Transform from GEV to unit GUMBEL
      for(k = 0; k < *nsite; k++)
	for(i = 0; i < *ndata; i++)
	  res[i + k * *ndata] = qgev(pgev(data[i + k * *ndata], eloc[k], escale[k], eshape[k]), 0, 1, 0);
      break;
    case 3: // Transform from GEV to unit WEIBULL
      for(k = 0; k < *nsite; k++)
	for(i = 0; i < *ndata; i++)
	  res[i + k * *ndata] = qgev(pgev(data[i + k * *ndata], eloc[k], escale[k], eshape[k]), 1, 1, -1);
      break;
    case 4: // Transform from GEV to GEV
      for(k = 0; k < *nsite; k++)
	for(i = 0; i < *ndata; i++)
	  res[i + k * *ndata] = qgev(pgev(data[i + k * *ndata], eloc[k], escale[k], eshape[k]), ploc[k], pscale[k], pshape[k]);
      break;
    case 5: // Transform from unit FRECHET to GEV
      for(k = 0; k < *nsite; k++)
	for(i = 0; i < *ndata; i++)
	  res[i + k * *ndata] = qgev(pgev(data[i + k * *ndata], 1, 1, 1), ploc[k], pscale[k], pshape[k]);
      break;
    case 6: // Transform from unit GUMBEL to GEV
      for(k = 0; k < *nsite; k++)
	for(i = 0; i < *ndata; i++)
	  res[i + k * *ndata] = qgev(pgev(data[i + k * *ndata], 0, 1, 0), ploc[k], pscale[k], pshape[k]);
      break;
    case 7: // Transform from unit WEIBULL to GEV
      for(k = 0; k < *nsite; k++)
	for(i = 0; i < *ndata; i++)
	  res[i + k * *ndata] = qgev(pgev(data[i + k * *ndata], 1, 1, -1), ploc[k], pscale[k], pshape[k]);
    }
   
  return;
}

double dgev(double x, double loc, double scale, double shape)
{
  double y=0.0, res=0.0;

  y = (x - loc) / scale;

  if(shape==0)
    res = exp(-exp(-y) - y) / scale;
  else
    res = exp(-pow(fmax(1 + shape * y, 0), - 1 / shape)) * 
      pow(fmax(1 + shape * y, 0), - 1 / shape - 1) / scale;

  return res;
}

void GevLogLik(double *data, int *ndata, double *par, double *res)
{
  int n=0;
  
  if(par[1] <= 0)
    {
      *res = LOW;
      return;
    }

  for(n = 0; n < *ndata; n++)
    *res += log(dgev(data[n], par[0], par[1], par[2]));

  if(!R_FINITE(*res))
    *res = LOW;

  return;
}

double pgev(double x, double loc, double scale, double shape)
{
  double y=0.0, result=0.0;

  y = (x - loc) / scale;

  if(shape==0)
    result = exp(-exp(-y));
  else
    result = exp(-pow(fmax(1 + shape * y, 0), - 1 / shape));

  return result;
}

double qgev(double x, double loc, double scale, double shape)
{
  double res=0.0;

  if(shape==0)
    res = loc - scale * log(-log(x));
  else
    res = loc + scale * (pow(-log(x), -shape) - 1) / shape;

  return res;
}

double d2norm(double x, double y, double rho)
{
  double res=0.0, omr=1-pow(rho,2);

  res=0.5*exp(-0.5*(pow(x,2)-2*rho*x*y+pow(y,2))/omr)/sqrt(omr)/M_PI;

  return res;
}
