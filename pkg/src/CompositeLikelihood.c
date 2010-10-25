#include "header.h"

// Composite log-likelihood for Gaussian models:

void CompLikelihood(double *coordx, double *coordy, int *corrmod, double *data, 
		    double *dista, double *lags, int *like, int *model, double *nuisance, 
		    int *ndata, int *nsite, double *par, double *res, int *type)
{
  int i=0, h=0, j=0, n=0; 
  double corr=0.0, s1=0.0, s12=0.0, s1s=0.0;

  s1 = nuisance[1] + nuisance[2];//set nugget + sill
  s1s = pow(s1, 2);
 
  for(i = 0; i < (*nsite - 1); i++)
    for(j = (i + 1); j < *nsite; j++)
      {
	//pairwise Euclidean distance
	if(lags[h] <= *dista)
	  {
	    corr = CorrelationFct(corrmod, lags[h], par); // pairwise correlation
	    for(n = 0; n < *ndata; n++)
	      *res += PairLikelihood(corr, like, nuisance, s1, s1s, 
				     data[(n + i * *ndata)], data[(n + j * *ndata)], type);
	  }
	h++;
      }    

  if(!R_FINITE(*res))
    *res = LOW;

  return;
}

// Pairwise log-likelihood for Gaussian model:

double PairLikelihood(double corr, int *like, double *nuisance, 
		      double s1, double s1s, double u, double v, int *type)
{
  double det=0.0, res=0.0, s12=0.0, vario=0.0;
  double u2=0.0, v2=0.0, lij=0.0, li=0.0, lj=0.0;
 
  switch(*like)
    {
    case 1:// Conditional likelihood:
      switch(*type)
      	{
      	case 2: // Conditional Pairwise likelihood for a given pair
	  s12 = nuisance[2] * corr; //sill*corr
	  det = s1s - R_pow(s12, 2);
	  u = u - nuisance[0]; //u-mean
	  v = v - nuisance[0]; //v-mean
	  u2 = R_pow(u, 2);
	  v2 = R_pow(v, 2);
	  res =  - log(2 * M_PI) - log(det) + log(s1) + 
	    (u2 + v2) * (0.5 / s1 - s1 / det) + 2 * s12 * u * v / det;
	  break;
	}
      break;
    case 3: // Marginal likelihood:
      switch(*type)
	{
	case 1:// Difference likelihood for given pair
	  vario = nuisance[1] + nuisance[2] * (1 - corr); //nugget+sill*(1-corr)
	  res = -.5 * (log(2 * M_PI) + log(vario) + pow(u - v ,2) / (2 * vario));
	  break;
	case 2: // Pairwise likelihood for given pair
	  s12 = nuisance[2] * corr; //sill*corr
	  det = s1s - pow(s12, 2);
	  u = u - nuisance[0]; //u-mean
	  v = v - nuisance[0]; //v-mean
	  res = -.5 * (2 * log(2 * M_PI) + log(det) + 
		       (s1 * (pow(u, 2) + pow(v, 2)) - 2 * s12 * u * v) / det);
	  break;
	}
      break;
    }  
  return res;
}


