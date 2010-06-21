#include "header.h"

// empirical variogram

void Empiric_Variogram(double *bins, double *coordx, double *coordy, double *data, double *lags, 
		       double *lenbins, double *maxdist, double *moments, int *npairs, int *nsite, 
		       int *nbins)
{
  int h=0, i=0, j=0, p=0;
  double lower=0.0, step=0.0;

  lower = Minima(lags, npairs);

  if(*maxdist == 0)
    *maxdist = Maxima(lags, npairs) / 2;

  step = (*maxdist - lower) / (*nbins - 1);
  bins[0] = lower;

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
		  moments[h] = moments[h] + pow(data[i] - data[j], 2) / 2;
		  lenbins[h] = lenbins[h] + 1;
		}
	  }
	p++;
      }

  return;
}


// squared differences

void Wls(double *bins, int *corrmod, double *par, int *nbins, double *moments, 
	 double *lenbins, double *nuisance, int *weighted, double *res)
{
  int h=0;
  double mean=0.0, vario=0.0;

  for(h = 0; h < (*nbins - 1); h++)
    if(lenbins[h])
      {
	mean = moments[h] / lenbins[h];
	vario = nuisance[0] + nuisance[1] * //nugget+sill*(1-corr) 
	  (1 - CorrelationFct(corrmod, (bins[h] + bins[h + 1]) / 2, par));
	if(vario)
	  {
	    if(*weighted == 0) *res = *res - pow(mean - vario, 2) * lenbins[h];
	    if(*weighted == 1) *res = *res - pow(mean / vario - 1, 2) * lenbins[h];
	  }
      }

  return;
}


