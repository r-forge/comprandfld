#include "header.h"

void Distances(double *coordx, double *coordy, int *nsite, double *lags)
{
  int i=0, j=0, h=0;

  for(i = 0; i < (*nsite - 1); i++)
    for(j = (i + 1); j < *nsite; j++)
      {
	lags[h] = pythag(coordx[i] - coordx[j], coordy[i] - coordy[j]);
	h++;
      }

  return;
}

double Maxima(double *x, int size)
{
  double res=0.0;
  int i=0;

  res = x[0];

  for(i = 1; i < size; i++)
    res = fmax(res, x[i]);

  return res;
}

double Minima(double *x, int size)
{
  double res=0.0;
  int i=0;

  res = x[0];

  for(i = 1; i < size; i++)
    res = fmin(res, x[i]);

  return res;
}
