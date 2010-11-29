#include "header.h"

// binned variogram:
void Binned_Variogram(double *bins, double *data, int *lenbins,  double *maxdist, 
		      double *moments, int *nbins, int *ndata, int *nsite)
{
  int h=0, i=0, j=0, n=0, p=0;
  double step=0.0;

  step = (*maxdist - *minimdista) / (*nbins - 1);
  bins[0] = *minimdista;

  //define bins:
  for(h = 1; h < *nbins; h++)
    bins[h] = bins[h - 1] + step;

  for(i = 0; i < (*nsite - 1); i++)
    for(j = (i + 1); j < *nsite; j++)
      {
	if(lags[p] <= *maxdist)
	  {
	    for(h = 0; h < (*nbins - 1); h++)
	      if((bins[h] <= lags[p]) && (lags[p] < bins[h + 1]))
		{
		  for(n = 0; n < *ndata; n++)
		    {
		      moments[h] += pow(data[n + i * *ndata] - data[n + j * *ndata], 2) / 2;
		      lenbins[h] += 1;
		    }
		}
	  }
	p++;
      }

  return;
}

// binned madogram:
void Binned_Madogram(double *bins, double *data, int *lenbins,  double *maxdist, 
		      double *moments, int *nbins, int *ndata, int *nsite)
{
  int h=0, i=0, j=0, n=0, p=0;
  double step=0.0;

  step = (*maxdist - *minimdista) / (*nbins - 1);
  bins[0] = *minimdista;

  //define bins:
  for(h = 1; h < *nbins; h++)
    bins[h] = bins[h - 1] + step;

  for(i = 0; i < (*nsite - 1); i++)
    for(j = (i + 1); j < *nsite; j++)
      {
	if(lags[p] <= *maxdist)
	  {
	    for(h = 0; h < (*nbins - 1); h++)
	      if((bins[h] <= lags[p]) && (lags[p] < bins[h + 1]))
		{
		  for(n = 0; n < *ndata; n++)
		    {
		      moments[h] += fabs(data[n + i * *ndata] - data[n + j * *ndata]) / 2;
		      lenbins[h] += 1;
		    }
		}
	  }
	p++;
      }

  return;
}

// variogram cloud:
void Cloud_Variogram(double *bins, double *data, int *lenbins, 
		      double *moments, int *ndata, int *nsite)
{
  int h=0, i=0, j=0, n=0;

  for(i = 0; i < (*nsite - 1); i++)
    for(j = (i + 1); j < *nsite; j++)
      {
	bins[h] = lags[h];
	for(n = 0; n < *ndata; n++)
	  moments[h] += pow(data[n + i * *ndata] - data[n + j * *ndata], 2) / 2;
	lenbins[h] = *ndata;
	h++;
      }

  return;
}

// madogram cloud:
void Cloud_Madogram(double *bins, double *data, int *lenbins, 
		      double *moments, int *ndata, int *nsite)
{
  int h=0, i=0, j=0, n=0;

  for(i = 0; i < (*nsite - 1); i++)
    for(j = (i + 1); j < *nsite; j++)
      {
	bins[h] = lags[h];
	for(n = 0; n < *ndata; n++)
	  moments[h] += fabs(data[n + i * *ndata] - data[n + j * *ndata]) / 2;
	lenbins[h] = *ndata;
	h++;
      }

  return;
}

// empirical variogram (main structure):
void Empiric_Variogram(double *bins, int *cloud, double *data, 
		       int *lenbins, double *maxdist, double *moments, 
		       int *ndata, int *nsite, int *nbins, int *type)
{
  int h=0, i=0, j=0, n=0, p=0;
  double step=0.0;

  if(*maxdist == 0)
    *maxdist = *maximdista;

  switch(*type)
    {
    case 1:
      switch(*cloud)
	{
	case 0:// BINNED VARIOGRAM
	  Binned_Variogram(bins, data, lenbins,  maxdist, moments, nbins, ndata, nsite);
	  break;
	case 1:// VARIOGRAM CLOUD
	  Cloud_Variogram(bins, data, lenbins, moments, ndata, nsite);
	  break;
	}
      break;
    case 2:
      switch(*cloud)
	{
	case 0:// BINNED VARIOGRAM
	  Binned_Madogram(bins, data, lenbins,  maxdist, moments, nbins, ndata, nsite);
	  break;
	case 1:// VARIOGRAM CLOUD
	  Cloud_Madogram(bins, data, lenbins, moments, ndata, nsite);
	  break;
	}
      break;
    case 3:
      switch(*cloud)
	{
	case 0:// BINNED VARIOGRAM
	  Binned_Madogram(bins, data, lenbins,  maxdist, moments, nbins, ndata, nsite);
	  break;
	case 1:// VARIOGRAM CLOUD
	  Cloud_Madogram(bins, data, lenbins, moments, ndata, nsite);
	  break;
	}
      break;
    case 4:
      switch(*cloud)
	{
	case 0:// BINNED VARIOGRAM
	  Binned_Madogram(bins, data, lenbins,  maxdist, moments, nbins, ndata, nsite);
	  break;
	case 1:// VARIOGRAM CLOUD
	  Cloud_Madogram(bins, data, lenbins, moments, ndata, nsite);
	  break;
	}
      break;
    }

  return;
}

// Least square method for Gaussian model:
void LeastSquare_G(double *bins, int *corrmod, double *par, int *nbins, 
		   double *moments, double *lenbins, double *nuisance, 
		   double *res)
{
  int h=0;
  double variogram=0.0, vario=0.0;

  for(h = 0; h < (*nbins - 1); h++)
    if(lenbins[h])
      {
	variogram = moments[h] / lenbins[h];
	vario = nuisance[0] + nuisance[1] * //nugget+sill*(1-corr) 
	  (1 - CorrelationFct(corrmod, (bins[h] + bins[h + 1]) / 2, par));
	if(vario)
	  *res = *res - pow(variogram - vario, 2) * lenbins[h];
      }

  return;
}

// Weighted least square method for Gaussian model:
void WLeastSquare_G(double *bins, int *corrmod, double *par, int *nbins, 
		    double *moments, double *lenbins, double *nuisance, 
		    double *res)
{
  int h=0;
  double variogram=0.0, vario=0.0;

  for(h = 0; h < (*nbins - 1); h++)
    if(lenbins[h])
      {
	variogram = moments[h] / lenbins[h];
	vario = nuisance[0] + nuisance[1] * //nugget+sill*(1-corr) 
	  (1 - CorrelationFct(corrmod, (bins[h] + bins[h + 1]) / 2, par));
	if(vario)
	  *res = *res - pow(variogram / vario - 1, 2) * lenbins[h];
      }

  return;
}
/*
// Least square method for max-stable models:
void LeastSquare_M(double *bins, int *corrmod, double *par, int *nbins, 
		   double *moments, double *lenbins, double *nuisance, 
		   double *res)
{
  int h=0;
  double variogram=0.0, extcoeff=0.0, rextcoeff=0.0;

  for(h = 0; h < (*nbins - 1); h++)
    if(lenbins[h])
      {
	variogram = moments[h] / lenbins[h];
	rextcoeff = (1 + 2 * variogram) / (1 - 2 * variogram);
	extcoeff = 1;
	vario = nuisance[0] + nuisance[1] * //nugget+sill*(1-corr) 
	  (1 - CorrelationFct(corrmod, (bins[h] + bins[h + 1]) / 2, par));
	*res = *res - pow(rextcoeff - extcoeff, 2) * lenbins[h];
      }

  return;
}

// Weighted least square method for max-stable models:
void WLeastSquare_M(double *bins, int *corrmod, double *par, int *nbins, 
		    double *moments, double *lenbins, double *nuisance, 
		    double *res, double *weights)
{
  int h=0;
  double variogram=0.0, extcoeff=0.0, rextcoeff=0.0;

  for(h = 0; h < (*nbins - 1); h++)
    if(lenbins[h])
      {
	variogram = moments[h] / lenbins[h];
	rextcoeff = (1 + 2 * variogram) / (1 - 2 * variogram);
	extcoeff = 1;
	vario = nuisance[0] + nuisance[1] * //nugget+sill*(1-corr) 
	  (1 - CorrelationFct(corrmod, (bins[h] + bins[h + 1]) / 2, par));
	*res = *res - pow(variogram / vario - 1, 2) * lenbins[h];
      }

  return;
}
*/

