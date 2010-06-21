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
		    Score_Gauss_Diff(corr, 1, flag, grc, score, npar, par, 
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
void GodambeMat_emp(double *coordx, double *coordy, int *corrmod, double *data, 
		    double *eps, int *flagcorr, int *flagnuis, double *lags, int *model, 
		    int *ndata, int *npar, int *nparc, int *nsite, double *parcorr, 
		    double *nuisance, double *godambe, int *type)
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
		    Score_Gauss_Diff(corr, 1, flagnuis, gradcorr, gradient, npar, nuisance,
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

void GodambeMat_teo(double *coordx, double *coordy, int *corrmod, double *dista,
                    double *eps, int *flagcorr, int *flagnuis, double *lags, 
		    int *model, int *npar, int *nparc, int *nsite, double *parcorr, 
		    double *nuisance, double *sens, double *vari, int *type)
{
  double mean=nuisance[0], nugget=nuisance[1],  sill=nuisance[2];
  double *gradcorr_ij, *gradient_ij, *gradcorr_lk, *gradient_lk;
  double corr_ij=0.0, corr_lk=0.0, crosscorr=0.0, *vario;
  int s=0, l=0, k=0 ,h=0, i=0, j=0, m=0, n=0;
  int ij=0, lk=0;

  switch(*model)// Compute the Sensitivity matrix
    {
    case 1:// Gaussian model 
      switch(*type)
	{
	case 1: // variability of the log-score difference likelihood
	  gradcorr_ij = (double *) R_alloc(*nparc, sizeof(double));
	  gradcorr_lk = (double *) R_alloc(*nparc, sizeof(double));
	  gradient_ij = (double *) R_alloc(*npar, sizeof(double)); 
	  gradient_lk = (double *) R_alloc(*npar, sizeof(double));
	  vario = (double *) R_alloc(6, sizeof(double));

	  for(i = 0; i < (*nsite - 1); i++) 
	    {
	      for(j = (i + 1); j < *nsite; j++)
		{
		  lk=0;
		  if(lags[ij] <= *dista)
		    {// Compute the correlation function
		      corr_ij = CorrelationFct(corrmod, lags[ij], parcorr);
		      // Compute the gradient of a given correlation model
		      GradientCorrFct(corr_ij, corrmod, eps, flagcorr, gradcorr_ij, lags[ij], parcorr);
		      Score_Gauss_Diff(corr_ij, 0, flagnuis, gradcorr_ij, gradient_ij, npar, nuisance, 0, 0);

		      vario[0] = Variogram(corrmod, lags[ij], nuisance, parcorr);

		      h=0;

		      for(m = 0; m < *npar; m++)
			for(n = m; n < *npar; n++)
			  {
			    sens[h] = sens[h] + 0.5 * gradient_ij[m] * gradient_ij[n];
			    h++;
			  }

		      for(l = 0; l < (*nsite - 1); l++)
			{
			  vario[2] = Variogram(corrmod, pythag(coordx[i] - coordx[l], coordy[i] - coordy[l]),
					       nuisance, parcorr);
			  vario[3] = Variogram(corrmod, pythag(coordx[j] - coordx[l], coordy[j] - coordy[l]),
					       nuisance, parcorr);

			  for(k = (l + 1); k < *nsite; k++)
			    {
			      if(lags[lk] <= *dista)
				{// Compute the correlation function
				  corr_lk = CorrelationFct(corrmod, lags[lk], parcorr);
				  // Compute the gradient of a given correlation model

				  GradientCorrFct(corr_lk, corrmod, eps, flagcorr, gradcorr_lk, lags[lk], parcorr);
				  Score_Gauss_Diff(corr_lk, 0, flagnuis, gradcorr_lk, gradient_lk, npar, nuisance, 0, 0);

				  vario[1] = Variogram(corrmod, lags[lk], nuisance, parcorr);
				  vario[4] = Variogram(corrmod, pythag(coordx[i] - coordx[k], coordy[i] - coordy[k]),
						       nuisance, parcorr);
				  vario[5] = Variogram(corrmod, pythag(coordx[j] - coordx[k], coordy[j] - coordy[k]),
						       nuisance, parcorr);

				  crosscorr=(R_pow(vario[2] - vario[3] - vario[4] + vario[5],2)) / (4 * vario[0] * vario[1]);
				  s = 0;

				  for(m = 0; m < *npar; m++)
				    for(n = m; n < *npar; n++)
				      {
					vari[s] = vari[s] +  0.5 * gradient_ij[m] * gradient_lk[n] * crosscorr;
				      s++;
				      }
				}
			      lk++;
			    }
			}
		    }
		  ij++;
		}
	    }
	  break;
	}
      break;
    }
  return;
}


/*
void GodambeMat_teo(double *coordx, double *coordy, int *corrmod, double *dista, double *eps, 
		    int *flagcorr, int *flagnuis, double *lags, int *model, int *npar,
                    int *nparc, int *nsite, double *parcorr, double *nuisance, double *sens, 
		    double *vari, int *type)
{
  double mean=nuisance[0], nugget=nuisance[1],  sill=nuisance[2];
  double *gradcorr_ij, *gradient_ij, *gradcorr_lk, *gradient_lk;
  double corr_ij=0.0, corr_lk=0.0 , crosscorr=0.0, *vario;
  int s=0, l=0, k=0 ,h=0, i=0, j=0, m=0, n=0;
  int ij=0, il=0, jl=0, lk=0, ik=0, jk=0;

  switch(*model)// Compute the Sensitivity matrix
    {
    case 1:// Gaussian model 
      switch(*type)
	{
	case 1: // variability of the log-score difference likelihood
	  gradcorr_ij = (double *) R_alloc(*nparc, sizeof(double));
	  gradcorr_lk = (double *) R_alloc(*nparc, sizeof(double));
	  gradient_ij = (double *) R_alloc(*npar, sizeof(double)); 
	  gradient_lk = (double *) R_alloc(*npar, sizeof(double));
	  vario = (double *) R_alloc(6, sizeof(double));

	  for(i = 0; i < (*nsite - 1); i++)
	    for(j = (i + 1); j < *nsite; j++)
	      {
		//lag_ij = pythag(coordx[i] - coordx[j], coordy[i] - coordy[j]);
                if(lags[ij] <= *dista)
		  {
		    corr_ij = CorrelationFct(corrmod, lags[ij], parcorr);// Compute the correlation function
		    // Compute the gradient of a given correlation model
		    GradientCorrFct(corr_ij, corrmod, eps, flagcorr, gradcorr_ij, lags[ij], parcorr);
		    Score_Gauss_Diff(corr_ij, 0, flagnuis, gradcorr_ij, gradient_ij, npar, nuisance, 0, 0);

		    vario[0] = Variogram(corrmod, lags[ij], nuisance, parcorr);

		    h=0;

		    for(m = 0; m < *npar; m++)
		      for(n = m; n < *npar; n++)
			{
			  sens[h] = sens[h] + 0.5 * gradient_ij[m] * gradient_ij[n];
			  h++;
			}

		    for(l = 0; l < (*nsite - 1); l++)
		      {
			vario[2] = Variogram(corrmod, lags[il], nuisance, parcorr);
			vario[3] = Variogram(corrmod, lags[jl], nuisance, parcorr);

			for(k = (l + 1); k < *nsite; k++)
			  {
			    //lag_lk = pythag(coordx[l] - coordx[k], coordy[l] - coordy[k]);
		            if(lags[lk] <= *dista)
			      {
				corr_lk = CorrelationFct(corrmod, lags[lk], parcorr);// Compute the correlation function
				// Compute the gradient of a given correlation model

				GradientCorrFct(corr_lk, corrmod, eps, flagcorr, gradcorr_lk, lags[lk], parcorr);
				Score_Gauss_Diff(corr_lk, 0, flagnuis, gradcorr_lk, gradient_lk, npar, nuisance, 0, 0);

				vario[1] = Variogram(corrmod, lags[lk], nuisance, parcorr);
				vario[4] = Variogram(corrmod, lags[ik], nuisance, parcorr);
				vario[5] = Variogram(corrmod, lags[jk], nuisance, parcorr);

				crosscorr=(R_pow(vario[2] - vario[3] - vario[4] + vario[5], 2)) / (4 * vario[0] * vario[1]);
				s = 0;

				for(m = 0; m < *npar; m++)
				  for(n = m; n < *npar; n++)
				    {
				      vari[s] = vari[s] +  0.5 * gradient_ij[m] * gradient_lk[n] * crosscorr;
				      s++;
				    }
				ik++;
				jk++;
			      }
			    lk++;
			  }
			il++;
			jl++;
		      }
		  }
		ij++;
	      }
	  break;
	}
      break;
    }
 return;
}
*/

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
void Score_Gauss_Diff(double corr, int expval, int *flag, double *gradcorr, 
		      double *gradient, int *npar, double *par, double u, double v)
{
  // Initialization variables:
  double mean=par[0], nugget=par[1], sill=par[2];
  double vario=0.0, sh=0.0;
  int h=0, i=0, j=0;

  vario = nugget + sill * (1 - corr);
  if(expval)
    sh = 0.5 * (0.5 * pow(u - v ,2) / vario - 1) / vario;
  else
    sh = 1 / vario;

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


void Sensitivity(double *coordx, double *coordy, int *corrmod, double *dista,
                 double *eps, int *flagcorr, int *flagnuis, int *model, int *npar, 
                 int *nparc, int *nsite, double *parcorr, double *nuisance, 
                 double *sens, int *type)
{
  // Initialization variables:
  double *gradcorr, *gradient;
  double corr=0.0, lag=0.0;
  int h=0, i=0, j=0, m=0, n=0;

  switch(*model)// Compute the Sensitivity matrix
    {
    case 1:// Gaussian model 
      switch(*type)
	{
	case 1: // Sensitivity of the log difference likelihood

	  gradcorr = (double *) R_alloc(*nparc, sizeof(double));
	  gradient = (double *) R_alloc(*npar, sizeof(double));

	  for(i = 0; i < (*nsite - 1); i++)
	    for(j = (i + 1); j < *nsite; j++)
	      {// Set the lag for given pair

		lag = pythag(coordx[i] - coordx[j], coordy[i] - coordy[j]);
		if(lag <= *dista)
		  {
		    corr = CorrelationFct(corrmod, lag, parcorr);// Compute the correlation function
		    // Compute the gradient of a given correlation model
		    GradientCorrFct(corr, corrmod, eps, flagcorr, gradcorr, lag, parcorr);
		    Score_Gauss_Diff(corr, 0, flagnuis, gradcorr, gradient, npar, nuisance, 0, 0);

		    h = 0;

		    for(m = 0; m < *npar; m++)
		      for(n = m; n < *npar; n++)
			{
			  sens[h] = sens[h] + 0.5 * gradient[m] * gradient[n];
			  h++;
			}
		  }
	      }
	  break;
	}
      break;
    }
  return;
}

