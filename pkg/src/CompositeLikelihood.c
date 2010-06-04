#include "header.h"

// Composite log-likelihood for extremal models:

void CompLikelihood(double *coordx, double *coordy, int *corrmod, double *data, 
		    int *model, double *nuisance, int *ndata, int *nsite, 
		    double *par, double *res, int *type)
{
  int i=0, j=0, n=0; 
  double corr=0.0, lag=0.0, s1=0.0, s12=0.0, s1s=0.0;

  s1 = nuisance[1] + nuisance[2];//set nugget + sill
  s1s = pow(s1, 2);
 
  for(i = 0; i < (*nsite - 1); i++)
    for(j = (i + 1); j < *nsite; j++)
      {
	//pairwise Euclidean distance
	lag = pythag(coordx[i] - coordx[j], coordy[i] - coordy[j]);
	corr = CorrelationFct(corrmod, lag, par); // pairwise correlation
	for(n = 0; n < *ndata; n++)
	  *res += PairLikelihood(corr, nuisance, s1, s1s, data[(n + i * *ndata)], 
				 data[(n + j * *ndata)], type);
      }    

  if(!R_FINITE(*res))
    *res = LOW;

  return;
}

// Pairwise log-likelihood for extremal models on frechet scale:

double PairLikelihood(double corr, double *nuisance, double s1, double s1s, 
		      double u, double v, int *type)
{
  double det=0.0, res=0.0, s12=0.0, vario=0.0;

  switch(*type)
    {
    case 1:// Difference likelihood for given pair
      vario = nuisance[1] + nuisance[2] * (1 - corr); //nugget+sill*(1-corr)
      res = -.5 * (log(vario) + pow(u - v ,2) / (2 * vario));
      break;
    case 2: // Pairwise likelihood for given pair
      s12 = nuisance[2] * corr; //sill*corr
      det = s1s - pow(s12, 2);
      u = u - nuisance[0]; //u-mean
      v = v - nuisance[0];
      res = -.5 * log(det) - .5 * (s1 * (pow(u, 2) + pow(v, 2)) - 
				   2 * s12 * u * v) / det;
      break;
    }
  
  return res;
}


