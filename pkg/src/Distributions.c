#include "header.h"

// Procedures are in alphabetical order.

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
  
  //  if((par[1] <= 0) || (par[2] < -1))
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

/*
double qgumbel(double x, double loc, double scale)
{
  double result=0.0;

  result = loc - scale * log(-log(x));

  return result;
}


double frechet2gev(double x, double alpha, double beta, double gamma,
		   double loc, double scale, double shape)
{
  double u=0.0, result=0.0;

  u = pgev(x, alpha, beta, gamma);

  result = qgev(u, loc, scale, shape);

  return result;

}
*/

/*
double gumbel2gev(double x, double loc, double scale, double shape)
{
  double y=0.0, result=0.0;

  y = exp(x);

  result = unitfrechet2gev(y, loc, scale, shape);

  return result;

}

double gev2gumbel(double x, double loc, double scale, double shape)
{
  double u=0.0, result=0.0;

  u = pgev(x, loc, scale, shape);
  // Transform to unit Gumbel
  result = qgumbel(u, 0, 1);

  return result;
}

double gev2unitfrechet(double x, double loc, double scale, double shape)
{
  double u=0.0, result=0.0;

  u = pgev(x, loc, scale, shape);
  // Transform to unit Frechet
  result = qgev(u, 1, 1, 1);

  return result;

}

double unitfrechet2gev(double x, double loc, double scale, double shape)
{
  double result=0.0;

  result = scale * (pow(x, 1 / shape) - 1) /  shape + loc;

  return result;

}
*/
