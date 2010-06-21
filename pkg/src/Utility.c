#include "header.h"
#define REARTH 6378.388


void Distances(double *coordx, double *coordy, double *lags, int *nsite, int *type)
{
  int i=0, j=0, h=0;


  if(*type)
    {
      for(i = 0; i < (*nsite - 1); i++)
	for(j = (i + 1); j < *nsite; j++)
	  {
	    lags[h] = Dist_geodesic(coordx[i],coordy[i], coordx[j],coordy[j]);
	    h++;
	  }
    }
  else
    {
      for(i = 0; i < (*nsite - 1); i++)
	for(j = (i + 1); j < *nsite; j++)
	  {
	    lags[h] = pythag(coordx[i] - coordx[j], coordy[i] - coordy[j]);
	    h++;
	  }
    }
  return;
}

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
