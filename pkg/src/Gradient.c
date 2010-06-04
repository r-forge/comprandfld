#include "header.h"

// Gradient of the composite-likelihood objects:

void CompScore(double *coordx, double *coordy, int *corrmod, double *data, 
	       double *eps, int *flag, int *flagcorr, int *model, int *ndata, 
	       int *ngrc, int *npar, int *nsite, double *par, double *parcorr, 
	       double *res, int *type, int *weight)
{
  int d=0, i=0, j=0, n=0;
  double corr=0.0, lag=0.0, *grc, *score;

  grc = (double *) R_alloc(*ngrc, sizeof(double));
  score = (double *) R_alloc(*npar, sizeof(double));

  for(n = 0; n < *ndata; n++)
    {
      for(i = 0; i < (*nsite - 1); i++)
        for(j = (i + 1); j < *nsite; j++)
	  {
            lag = pythag(coordx[i] - coordx[j], coordy[i] - coordy[j]);
	    corr = CorrelationFct(corrmod, lag, parcorr);
	    GradientCorrFct(corr, corrmod, eps, flagcorr, grc, lag, parcorr);
	    switch(*model)
	      {
	      case 1:// Gaussian models
		switch(*type)
		  {
		  case 1: // Difference score for given pair
		    Score_Gauss_Diff(corr, flag, grc, score, npar, par, 
				     data[(n + i * *ndata)], data[(n + j * *ndata)]);
		    break;
		  case 2:// Pairwise score for given pair
		    break;
		  }
		break;
	      }
	    // Summation of the pairwise gradients:
	    for(d = 0; d < *npar; d++)
	      res[d] = res[d] + score[d];

	  }
    }

  return;
}

// Estimation of the Senstive (H) and Variability (J) components of
// the Godambe matrix
void GodambeMat(double *coordx, double *coordy, int *corrmod, double *data, 
		double *eps, int *flagcorr, int *flagnuis, int *model, int *ndata, 
		int *npar, int *nparc, int *nsite, double *parcorr, double *nuisance, 
		double *godambe, int *type)
{
  int d=0, i=0, j=0, k=0, n=0, nmat=0;
  double corr, *gradcorr, *gradient, lag, *score;

  nmat = pow(*npar, 2);// Set the dimension of the Godambe matrix

  gradcorr = (double *) R_alloc(*nparc, sizeof(double));
  gradient = (double *) R_alloc(*npar, sizeof(double));
  score = (double *) R_alloc(*npar, sizeof(double));

  for(n = 0; n < *ndata; n++)
    {
      for(i = 0; i < *npar; i++)// Initialize the gradient vector
        score[i] = 0;

      for(i = 0; i < (*nsite - 1); i++)
        for(j = (i + 1); j < *nsite; j++)
	  {// Set the lag for given pair
	    lag = pythag(coordx[i] - coordx[j], coordy[i] - coordy[j]);
	    corr = CorrelationFct(corrmod, lag, parcorr);// Compute the correlation function
	    // Compute the gradient of a given correlation model
	    GradientCorrFct(corr, corrmod, eps, flagcorr, gradcorr, lag, parcorr);
	    switch(*model)// Compute the gradient of a log likelihood object
	      {
	      case 1:// Gaussian model 
		switch(*type)
		  {
		  case 1: // Gradient of the log difference likelihood
		    Score_Gauss_Diff(corr, flagnuis, gradcorr, gradient, npar, nuisance,
				     data[(n + i * *ndata)], data[(n + j * *ndata)]);
		    break;
		  case 2:// Gradient of the log pairwise likelihood
		    Score_Gauss_Pair(corr, flagnuis, gradcorr, gradient, npar, nuisance, 
				     data[(n + i * *ndata)], data[(n + j * *ndata)]);
		    break;
		  }
		break;
	      }
	    // Set the sensitivity matrix:
	    for(d = 0; d < *npar; d++)
	      {
		score[d] = score[d] + gradient[d];

		for(k = 0; k < *npar; k++)
		  godambe[d * *npar + k] = godambe[d * *npar + k] +
		    gradient[d] * gradient[k];
	      }
	  }
      // Set the variability matrix:
      for(i = 0; i < *npar; i++)
	for(j = 0; j < *npar; j++)
	  godambe[(i * *npar + j) + nmat] = godambe[(i * *npar + j) + nmat] + 
	    score[i] * score[j];
    }

  return;
}

// Score of the Gaussian model via pairwise:
void Score_Gauss_Pair(double corr, int *flag, double *gradcorr, double *gradient, int *npar,
		      double *par, double u, double v)
{
  // Initialization variables:
  double mean=par[0], nugget=par[1], sill=par[2];
  double a=0, b=0, d=0;
  int h=0, i=0, j=0;

  d = pow(nugget + sill, 2) - pow(sill * corr, 2);
  a = (nugget + sill) / d;
  b = sill * corr / d;

  // Derivative of the pairwise respect with the mean
  if(flag[0] == 1)
    {
      gradient[i] = (u + v - 2 * mean) / ((1 + corr) * sill + nugget);
      i++;
    }

  if(flag[1] == 1)
    {
      gradient[i] = (pow(u, 2) + pow(v, 2)) * (pow(a, 2) - .5 / d) -
	a * (2 * b * u * v + 1);
      i++;
    }

  if(flag[2] == 1)
    {
      gradient[i] = -a + (pow(u, 2) + pow(v, 2)) * (a * (a - pow(b, 2) * d / sill) - .5 / d) +
	b * (u * v * (1 / sill + 2 * (pow(b, 2) * d / sill - a)) + b * d / sill);
      i++;
    }

  for(j = i; j < *npar; j++)
    {
      gradient[j] = b * (b * d * (1 - a * (pow(u, 2) + pow(v, 2))) + u * v * (1 + 2 * pow(b, 2) * d))  * 
	gradcorr[h] / corr;
      h++;
    }

  return;
}

// Score of the Gaussian model via difference:
void Score_Gauss_Diff(double corr, int *flag, double *gradcorr, double *gradient, int *npar,
		      double *par, double u, double v)
{
  // Initialization variables:
  double mean=par[0], nugget=par[1], sill=par[2];
  double vario=0.0, sh=0.0;
  int h=0, i=0, j=0;

  vario = nugget + sill * (1 - corr);
  sh = 0.5 * (0.5 * pow(u - v ,2) / vario - 1) / vario;

  if(flag[1] == 1)
    {
      gradient[i] = sh;
      i++; 
    }

  if(flag[2] == 1)
    {
      gradient[i] = (1 - corr) * sh;
      i++;
    }

  for(j = i; j < *npar; j++)
    {
      gradient[j] = - sill * gradcorr[h] * sh;
      h++;
    }

  return;
}

