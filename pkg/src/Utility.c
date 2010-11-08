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
	    lags[h] = Dist_geodesic(coordx[i], coordy[i], coordx[j], coordy[j]);
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

void SetSampling(double *coordx, double *coordy, double *data, int n,
		 int *ndata, int *npts, double *scoordx, double *scoordy, 
		 double *sdata, int *size, double xmax, double xmin, 
		 double ymax, double ymin)
{
  int i=0, j=0;

  for(i = 0; i < *size; i++)
    if((xmin <= coordx[i]) && (coordx[i] <= xmax) && 
       (ymin <= coordy[i]) && (coordy[i] <= ymax))
      {
	scoordx[j] = coordx[i];
	scoordy[j] = coordy[i];
	sdata[j] = data[n * *ndata + i];
	j++;
      }
  *npts = j;

  return;
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
