#include "header.h"

/*###################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Email: simone.padoan@epfl.ch.
### Institute: EPFL.
### File name: Godambe.c
### Description:
### This file contains a set of procedures
### for the computation of the Godambe matrix of
### random fields.
### Last change: 25/10/2010.
##################################################*/


// Compute the score vector of a random field:
void CompScore(double *coordx, double *coordy, int *corrmod, double *data,
	       double *eps, int *flag, int *flagcorr, int *model, int *like,
               int *ndata, int *ngrc, int *npar, int *nsite, double *par,
               double *parcorr, double *res, int *type, int *weight)
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
		    Grad_Diff_Gauss(corr, flag, grc, score, npar, par,
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

// Empirical estimation of the Senstive (H) and Variability (J) components of
// the Godambe matrix:
void GodambeMat_emp(double *coordx, double *coordy, int *corrmod, double *data,
		    double *eps, int *flagcorr, int *flagnuis, int *like, int *model,
                    int *ndata, int *npar, int *nparc, int *nsite, double *parcorr,
                    double *nuisance, double *godambe, double *score, int *type)
{
  int d=0, i=0, j=0, k=0, n=0, nmat=0;
  double corr, *gradcorr, *gradient, lag; //*score;

  nmat = pow(*npar, 2);// Set the dimension of the Godambe matrix

  gradcorr = (double *) R_alloc(*nparc, sizeof(double));
  gradient = (double *) R_alloc(*npar, sizeof(double));
  //  score = (double *) R_alloc(*npar, sizeof(double));

  for(n = 0; n < *ndata; n++)
    {
      //     for(i = 0; i < *npar; i++)// Initialize the gradient vector
	//   score[i] = 0;

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
		switch(*like)
		  {
		  case 1:// Conditional likelihood:
		    break;
		  case 3: // Marginal likelihood:
		    switch(*type)
		      {
		      case 1: // Gradient of the log difference likelihood
			Grad_Diff_Gauss(corr, flagnuis, gradcorr, gradient, npar, nuisance,
					data[(n + i * *ndata)], data[(n + j * *ndata)]);
			break;
		      case 2:// Gradient of the log pairwise likelihood
			Grad_Pair_Gauss(corr, flagnuis, gradcorr, gradient, npar, nuisance,
					data[(n + i * *ndata)], data[(n + j * *ndata)]);

			break;
		      }
		    break;
		  }
		break;
	      }
	    // Set the sensitivity matrix:
	    for(d = 0; d < *npar; d++)
	      {
		score[d] = score[d] + gradient[d];

		//		for(k = 0; k < *npar; k++)
		//	  godambe[d * *npar + k] = godambe[d * *npar + k] +
		//	    gradient[d] * gradient[k];
	      }
	  }
      // Set the variability matrix:
      //      for(i = 0; i < *npar; i++)
      //	for(j = 0; j < *npar; j++)
      //	  godambe[(i * *npar + j) + nmat] = godambe[(i * *npar + j) + nmat] +
      //	    score[i] * score[j];
    }

  return;
}

// The exact Senstive (H) and Variability (J) components of
// the Godambe matrix:
void GodambeMat(double *coordx, double *coordy, int *corrmod, double *data, double *dista, 
		double *eps, int *flagcorr, int *flagnuis, double *lags, int *like, 
		int *lonlat, int *model, int *npar, int *nparc, int *nsite, double *parcorr, 
		double *nuisance, double *sensmat, int *type, double *varimat, int *vartype, 
		double *winc)
{
  int *npair;

  npair = (int *) R_alloc(1, sizeof(int));
  *npair = 0;

  //---------- COMPUTATION OF THE GODAMBE MATRIX ---------//
  switch(*vartype)
    {
    case 1://------------ START THEORETICAL COMPUTATION ------------//
      switch(*like)
	{
	case 1://----------- CONDITIONAL LIKELIHOOD ------------------------//
	  //------- GODAMBE FOR THE CONDITIONAL CASE -------//
	  break;
	case 3://----------- MARGINAL LIKELIHOOD ---------------------------//
	  switch(*type)
	    {
	    case 1://------- GODAMBE FOR THE DIFFERENCE CASE -------//
	      GodambeMat_Diff(coordx, coordy, corrmod, dista, eps, flagcorr, flagnuis,
			      lags, model, npar, nparc, nsite, parcorr, nuisance,
			      sensmat, varimat);
	      break;
	    case 2://------- GODAMBE FOR THE PAIRWISE CASE -------//
	      break;
	    }
	  break;
	}
      break;//------------ END THEORETICAL COMPUTATION ------------//
    case 2://------------ START SUB-SAMPLE ESTIMATION ------------//
      Sensitivity(coordx, coordy, corrmod, dista, eps, flagcorr,
		  flagnuis, lags, like, model, npair, npar, nparc,
		  nsite, parcorr, nuisance, sensmat, type);
      Vari_SubSamp(coordx, coordy, corrmod, data, dista, eps, flagcorr, 
		   flagnuis, like, lonlat, npair, npar, nparc, nsite, 
		   nuisance, parcorr, type, varimat, winc);
      break;//------------ END SUB-SAMPLE ESTIMATION ------------//
    case 3://------------ START EMPIRICAL ESTIMATION ------------//
      break;//------------ END EMPIRICAL ESTIMATION ------------//
    }
  return;
}


void GodambeMat_Diff(double *coordx, double *coordy, int *corrmod, double *dista, double *eps,
		     int *flagcorr, int *flagnuis, double *lags, int *model, int *npar, int *nparc,
		     int *nsite, double *parcorr, double *nuisance, double *sensmat, double *varimat)
{

  double *gradcorr_ij, *gradient_ij, *gradcorr_lk, *gradient_lk;
  double corr_ij=0.0, corr_lk=0.0, crosscorr=0.0, *vario;
  int s=0, l=0, k=0 , i=0, j=0, m=0, n=0;
  int ij=0, lk=0;

  gradcorr_ij = (double *) R_alloc(*nparc, sizeof(double));
  gradcorr_lk = (double *) R_alloc(*nparc, sizeof(double));
  gradient_ij = (double *) R_alloc(*npar, sizeof(double));
  gradient_lk = (double *) R_alloc(*npar, sizeof(double));
  vario = (double *) R_alloc(6, sizeof(double));

  //--------------- GAUSSIAN MODEL --------------------------------//
  for(i = 0; i < (*nsite - 1); i++)
    {
      for(j = (i + 1); j < *nsite; j++)
	{
	  lk=0;
	  if(lags[ij] <= *dista)
	    {// Compute the correlation function
	      corr_ij = CorrelationFct(corrmod, lags[ij], parcorr);
	      // Compute the gradient for a given correlation value
	      GradientCorrFct(corr_ij, corrmod, eps, flagcorr, gradcorr_ij, lags[ij], parcorr);
	      // Compute the gradient of the variogram for the difference Gaussian likelihood
	      Grad_Diff_Vario(corr_ij, flagnuis, gradcorr_ij, gradient_ij, npar, nuisance);
	      // COMPUTATION OF THE SENSITIVITY MATRIX
	      Sens_Diff_Gauss_ij(gradient_ij, npar, sensmat);

	      // START COMPUTATION OF THE VARIABILITY MATRIX
	      vario[0] = Variogram(corrmod, lags[ij], nuisance, parcorr);

	      for(l = 0; l < (*nsite - 1); l++)
		{// Insert the case of great circle distance
		  vario[2] = Variogram(corrmod, pythag(coordx[i] - coordx[l], coordy[i] - coordy[l]), nuisance, parcorr);
		  vario[3] = Variogram(corrmod, pythag(coordx[j] - coordx[l], coordy[j] - coordy[l]), nuisance, parcorr);

		  for(k = (l + 1); k < *nsite; k++)
		    {
		      if(lags[lk] <= *dista)
			{// Compute the correlation function
			  corr_lk = CorrelationFct(corrmod, lags[lk], parcorr);
			  // Compute the gradient of a given correlation model

			  GradientCorrFct(corr_lk, corrmod, eps, flagcorr, gradcorr_lk, lags[lk], parcorr);
			  // Compute the gradient of the variogram for the LK difference Gaussian likelihood
			  Grad_Diff_Vario(corr_lk, flagnuis, gradcorr_lk, gradient_lk, npar, nuisance);

			  vario[1] = Variogram(corrmod, lags[lk], nuisance, parcorr);
			  vario[4] = Variogram(corrmod, pythag(coordx[i] - coordx[k], coordy[i] - coordy[k]), nuisance, parcorr);
			  vario[5] = Variogram(corrmod, pythag(coordx[j] - coordx[k], coordy[j] - coordy[k]), nuisance, parcorr);

			  crosscorr = (R_pow(vario[2] - vario[3] - vario[4] + vario[5],2)) / (4 * vario[0] * vario[1]);
			  s = 0;

			  for(m = 0; m < *npar; m++)
			    for(n = m; n < *npar; n++)
			      {
				varimat[s] = varimat[s] +  0.5 * gradient_ij[m] * gradient_lk[n] * crosscorr;
				s++;
			      }
			  // END COMPUTATION OF THE VARIABILITY MATRIX
			}
		      lk++;
		    }
		}
	    }
	  ij++;
	}
    }

  return;
}


// Compute the gradient vector of the conditional pairwise log likelihood for a Gaussian model:
void Grad_Cond_Gauss(double corr, int *flag, double *gradcorr, double *gradient, int *npar,
		     double *par, double u, double v)
{
  // Initialization variables:
  double mean=par[0], nugget=par[1], sill=par[2];
  double C=0.0, L=0.0, M=0.0, Q=0.0, R=0.0, S=0.0;
  double a=0.0, b=0.0, c=0.0, d=0.0, da=0.0, e=0.0;
  double m=0.0, pu=0.0, pv=0.0, q=0.0, su=0.0, sv=0.0;
  double suv=0.0, uv=0.0;
  int h=0, i=0, j=0;

  a = nugget + sill;
  b = sill * corr;
  c = pow(a, 2) - pow(b, 2);
  d = a / c;
  da = 2 * a;
  e = b / c;
  u = u - mean;
  v = v - mean;
  uv = u + v;
  pu = pow(u, 2);
  pv = pow(v, 2);
  R = pu + pv;
  S = .5 * R / c;
  L = u * v;
  M = L / c;
  Q = d * R - 2 * e * L - 1;
  su = (-1 + pu / a) / da; //first statistics: first component
  sv = (-1 + pv / a) / da; //second statistics: second component
  suv = su + sv;

  // Derivatives of the conditional respect with the mean
  if(flag[0] == 1)
    {
      gradient[i] = uv * (2 / (a + b) - 1 / a);
      i++;
    }
  // Derivative of the conditional respect with the nugget
  if(flag[1] == 1)
    {
      gradient[i] = 2 * (d * Q - S) - suv;
      i++;
    }
  // Derivative of the conditional respect with the sill
  if(flag[2] == 1)
    {
      gradient[i] = 2 * ((d - e) * Q - S - corr * M) - suv;
      i++;
    }
  // Derivative of the conditional respect with the correlation parameters
  h = 0;
  C = sill * (M - e * Q);
  for(j = i; j < *npar; j++)
    {
      gradient[j] = 2 * C  * gradcorr[h];
      h++;
    }

  return;
}


// Compute the gradient vector of the difference log likelihood for a Gaussian model :
void Grad_Diff_Gauss(double corr, int *flag, double *gradcorr,
		     double *gradient, int *npar, double *par, double u, double v)
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


// Compute the gradient vector of the variogram for a Gaussian model :
void Grad_Diff_Vario(double corr, int *flag, double *gradcorr,
		     double *gradient, int *npar, double *par)
{
  // Initialization variables:
  double mean=par[0], nugget=par[1], sill=par[2];
  double vario=0.0, sh=0.0;
  int h=0, i=0, j=0;

  vario = nugget + sill * (1 - corr);
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

// Compute the gradient vector of the pairwise log likelihood for a Gaussian model:
void Grad_Pair_Gauss(double corr, int *flag, double *gradcorr, double *gradient, int *npar,
		      double *par, double u, double v)
{
  // Initialization variables:
  double mean=par[0], nugget=par[1], sill=par[2];
  double C=0.0, L=0.0, M=0.0, Q=0.0, R=0.0, S=0.0;
  double a=0.0, b=0.0, c=0.0, d=0.0, e=0.0;
  int h=0, i=0, j=0;

  a = nugget + sill;
  b = sill * corr;
  c = pow(a, 2) - pow(b, 2);
  d = a / c;
  e = b / c;
  u = u - mean;
  v = v - mean;
  R = pow(u, 2) + pow(v, 2);
  S = .5 * R / c;
  L = u * v;
  M = L / c;
  Q = d * R - 2 * e * L - 1;

  // Derivative of the pairwise respect with the mean
  if(flag[0] == 1)
    {
      gradient[i] = (u + v) / (a + b);
      i++;
    }
  // Derivative of the pairwise respect with the nugget
  if(flag[1] == 1)
    {
      gradient[i] = d * Q - S;
      i++;
    }
  // Derivative of the pairwise respect with the sill
  if(flag[2] == 1)
    {
      gradient[i] = (d - e) * Q - S - corr * M;
      i++;
    }
  // Derivative of the pairwise respect with the correlation parameters
  C = sill * (M - e * Q);
  for(j = i; j < *npar; j++)
    {
      gradient[j] = C  * gradcorr[h];
      h++;
    }

  return;
}

// Compute the Sensitivity matrix of a random field:
void Sensitivity(double *coordx, double *coordy, int *corrmod, double *dista, 
		 double *eps, int *flagcorr, int *flagnuis, double *lags, 
		 int *like, int *model, int *npair, int *npar, int *nparc, 
		 int *nsite, double *parcorr, double *nuisance, double *sensmat, 
		 int *type)
{
  // Initialization variables:
  double *gradcorr, *gradient;
  double corr=0.0, lag=0.0;
  int h=0, i=0, j=0, m=0, n=0;

  switch(*model)// Compute the Sensitivity matrix
    {
    case 1:// Gaussian model
      switch(*like)
	{
	case 1:// Conditional likelihood:
	  Sens_Cond_Gauss(corrmod, dista, eps, flagcorr,
			  flagnuis, lags, nsite, nuisance,
			  npair, npar, nparc, parcorr, sensmat);
	  break;
	case 3: // Marginal likelihood:
	  switch(*type)
	    {
	    case 1: // Sensitivity for the difference likelihood case
	      Sens_Diff_Gauss(corrmod, dista, eps, flagcorr,
			      flagnuis, lags, nsite, nuisance,
			      npair, npar, nparc, parcorr, sensmat);
	      break;
	    case 2: // Sensitivity for the pairwise likelihood case
	      Sens_Pair_Gauss(corrmod, dista, eps, flagcorr,
			      flagnuis, lags, nsite, nuisance,
			      npair, npar, nparc, parcorr, sensmat);
	      break;
	    }
	  break;
	}
      break;
    }
  return;
}

// Compute the Sensitivity matrix for the composite Gaussian difference likelihood:
void Sens_Diff_Gauss(int *corrmod, double *dista, double *eps, int *flagcorr,
		     int *flagnuis, double *lags, int *nsite, double *nuisance,
		     int *npair, int *npar, int *nparc, double *parcorr, double *sensmat)
{
  // Initialization variables:
  int i=0, j=0, h=0;
  double corr=0.0, *gradcorr, *gradient;

  gradcorr = (double *) R_alloc(*nparc, sizeof(double));// Correlation gradient
  gradient = (double *) R_alloc(*npar, sizeof(double));// Overall gradient

  for(i = 0; i < (*nsite - 1); i++)
    for(j = (i + 1); j < *nsite; j++)
      {
	if(lags[h] <= *dista)
	  {
	    // Compute the correlation function for the elements i,j
	    corr = CorrelationFct(corrmod, lags[h], parcorr);
	    // Compute the gradient for the given correlation
	    GradientCorrFct(corr, corrmod, eps, flagcorr, gradcorr, lags[h], parcorr);
	    // Compute the gradient of the variogram for the given pair
	    Grad_Diff_Vario(corr, flagnuis, gradcorr, gradient, npar, nuisance);
	    // ADD TO THE SENSITIVITY MATRIX THE CONTRIBUTE OF THE GIVEN PAIR
	    Sens_Diff_Gauss_ij(gradient, npar, sensmat);
	    *npair = *npair + 1;
	  }
	h++;
      }

  return;
}


// Compute the Sensitivity matrix for a single Gaussian difference likelihood component:
void Sens_Diff_Gauss_ij(double *gradient, int *npar, double *sensmat)
{
  // Initialization variables:
  int h=0,i=0, j=0;

  for(i = 0; i < *npar; i++)
    for(j = i; j < *npar; j++)
      {
	sensmat[h] = sensmat[h] - 0.5 * gradient[i] * gradient[j];
	h++;
      }

  return;
}

// Compute the Sensitivity matrix for the pairwise composite Gaussian likelihood:
void Sens_Pair_Gauss(int *corrmod, double *dista, double *eps, int *flagcorr,
		     int *flagnuis, double *lags, int *nsite, double *nuisance,
		     int *npair, int *npar, int *nparc, double *parcorr, double *sensmat)
{
  // Initialization variables:
  int h=0, i=0, l=0, nsens=0, j=0;
  double corr=0.0, *gradcorr, *gradient, *sens;

  nsens = *npar * (*npar + 1) / 2;
  gradcorr = (double *) R_alloc(*nparc, sizeof(double));// Correlation gradient
  gradient = (double *) R_alloc(*npar, sizeof(double));// Overall gradient
  sens = (double *) R_alloc(nsens, sizeof(double));// One sensitive contribute

  for(i = 0; i < (*nsite - 1); i++)
    for(j = (i + 1); j < *nsite; j++)
      {
	if(lags[h] <= *dista)
	  {
	    // Compute the correlation function for the elements i,j
	    corr = CorrelationFct(corrmod, lags[h], parcorr);
	    // Compute the gradient for the given correlation
	    GradientCorrFct(corr, corrmod, eps, flagcorr, gradcorr, lags[h], parcorr);
	    // ADD TO THE SENSITIVITY MATRIX THE CONTRIBUTE OF THE GIVEN PAIR
	    Sens_Pair_Gauss_ij(corr, flagnuis, gradcorr, npar, nparc, nuisance, sens);
	    *npair = *npair + 1;
	    for(l = 0; l < nsens; l++)
	      sensmat[l] = sensmat[l] + sens[l];
	  }
	h++;
      }

  return;
}

// Compute the Sensitivity matrix for the Gaussian difference likelihood:
void Sens_Pair_Gauss_ij(double corr, int *flag, double *gradcorr, int *npar,
			int *nparc, double *par, double *sensmat)
{
  // Initialization variables:
  double mean=par[0], nugget=par[1], sill=par[2];
  double C=0.0, D=0.0, L=0.0, M=0.0, Q=0.0, R=0.0;
  double S=0.0, ER=0.0, EL=0.0, ES=0.0, EM=0.0, EQ=0.0;
  double a=0.0, b=0.0, c=0.0, d=0.0, e=0.0;
  double f=0.0, g=0.0, q=0.0, t=0.0;
  int h=0, i=0, l=0, j=0, p=0;

  a = nugget + sill;
  b = sill * corr;
  c = pow(a, 2) - pow(b, 2);
  d = a / c;
  e = b / c;
  f = a + b;
  q = d - corr * e;
  t = 1 - pow(corr, 2);

  // Expected values:
  ER = 2 * a;
  EL = b;
  ES = .5 * ER / c;
  EM = EL / c;
  EQ = d * ER - 2 * e * EL - 1;


  //---- START COMPUTATION THE SENSITIVITY MATRIX-----//

  // Derivatives of the PLLik respect with the mean
  if(flag[0] == 1)
    {
      sensmat[p] = - 2 / f; // Second derivative
       i++; p++; l++;
      // Mixed derivatives:
      if(flag[1] == 1) // nugget
	{
	  sensmat[p] = 0;
	  p++; l++;
	}
      if(flag[2] == 1) // sill
	{
	  sensmat[p] = 0;
	  p++; l++;
	}
      for(j = l; j < *npar; j++) // correlation parameters
	{
	  sensmat[p] = 0;
	  p++;
	}
    }
  // Derivative of the PLLik respect with the nugget
  if(flag[1] == 1)
    {
      l = i; i++;
      sensmat[p] = EQ * (1 / c - 4 * pow(d, 2)) + (2 * d * ER - 1) / c;
      p++; l++;
      if(flag[2] == 1) // sill
	{
	  sensmat[p] = 2 * (ES * (3 * d - corr * e) - EM * (d * corr + e) -
			 d * (d - corr * e) * (2 * EQ + 1)) - 1 / c;
	  p++; l++;
	}
      h = 0;
      for(j = l; j < *npar; j++) // correlation parameters
	{
	  sensmat[p] = 2 * sill * (d * e * (2 * EQ + 1) - e * ES - d * EM) * gradcorr[h];
	  h++; p++;
	}
    }
  // Derivative of the PLLik respect with the sill
  if(flag[2] == 1)
    {
      l = i; i++;
      sensmat[p] = - 2 * pow(q, 2) * (2 * EQ + 1) + 2 * (ES * (2 * q + d * t) -
						      EM * (2 * corr * q + e * t)) - t / c;
      p++; l++;
      h = 0;
      for(j = l; j < *npar; j++) // correlation parameters
	{
	  sensmat[p] = (2 * e * (sill * (q * (2 * EQ + 1) - ES) - EQ + e * b) +
		     EM * (1 - 2 * sill * q))  * gradcorr[h];
	  h++; p++;
	}
    }
  // Derivative of the PLLik respect with the correlation parameters*/
  D =  - pow(sill / c, 2) * (2 * pow(b,2) + c);

  for(j = 0; j < *nparc; j++)
    for(l = j; l < *nparc; l++)
      {
	sensmat[p] = D * gradcorr[j] * gradcorr[l];
	p++;

      }

  //---- END COMPUTATION OF THE SENSITIVITY MATRIX-----//

  return;
}

// Compute the Sensitivity matrix for the conditional composite Gaussian likelihood:
void Sens_Cond_Gauss(int *corrmod, double *dista, double *eps, int *flagcorr,
		     int *flagnuis, double *lags, int *nsite, double *nuisance,
		     int *npair, int *npar, int *nparc, double *parcorr, double *sensmat)
{
  // Initialization variables:
  int h=0, i=0, l=0, nsens=0, j=0;
  double corr=0.0, *gradcorr, *gradient, *sens;

  nsens = *npar * (*npar + 1) / 2;
  gradcorr = (double *) R_alloc(*nparc, sizeof(double));// Correlation gradient
  gradient = (double *) R_alloc(*npar, sizeof(double));// Overall gradient
  sens = (double *) R_alloc(nsens, sizeof(double));// One sensitive contribute

  for(i = 0; i < (*nsite - 1); i++)
    for(j = (i + 1); j < *nsite; j++)
      {
	if(lags[h] <= *dista)
	  {
	    // Compute the correlation function for the elements i,j
	    corr = CorrelationFct(corrmod, lags[h], parcorr);
	    // Compute the gradient for the given correlation
	    GradientCorrFct(corr, corrmod, eps, flagcorr, gradcorr, lags[h], parcorr);
	    // ADD TO THE SENSITIVITY MATRIX THE CONTRIBUTE OF THE GIVEN PAIR
	    Sens_Cond_Gauss_ij(corr, flagnuis, gradcorr, npar, nparc, nuisance, sens);
	    *npair = *npair + 1;
	    for(l = 0; l < nsens; l++)   
	      sensmat[l] = sensmat[l] + sens[l];
	  }
	h++;
      }

  return;
}


// Compute the Sensitivity matrix for a single Gaussian conditional likelihood:
void Sens_Cond_Gauss_ij(double corr, int *flag, double *gradcorr, int *npar,
			int *nparc, double *par, double *sensmat)
{
 // Initialization variables:
 double mean=par[0], nugget=par[1], sill=par[2];
 double C=0.0, D=0.0, L=0.0, M=0.0, Q=0.0, R=0.0, S=0.0;
 double ER=0.0, EL=0.0, ES=0.0, EM=0.0, EQ=0.0;
 double a=0.0, b=0.0, c=0.0, d=0.0;
 double e=0.0, f=0.0, g=0.0, q=0.0, t=0.0;
 int h=0, i=0, l=0, j=0, p=0;

 a = nugget + sill;
 b = sill * corr;
 c = pow(a, 2) - pow(b, 2);
 d = a / c;
 e = b / c;
 f = a + b;
 q = d - corr * e;
 t = 1 - pow(corr, 2);

 // Expected values:
 ER = 2 * a;
 EL = b;
 ES = .5 * ER / c;
 EM = EL / c;
 EQ = d * ER - 2 * e * EL - 1;


 //---- START COMPUTATION THE SENSITIVITY MATRIX-----//

 // Derivatives of the PLLik respect with the mean
 if(flag[0] == 1)
   {
     sensmat[p] = - 4 / f + 2/a; // Second derivative
      i++; p++; l++;
     // Mixed derivatives:
     if(flag[1] == 1) // nugget
	{
	  sensmat[p] = 0;
	  p++; l++;
	}
     if(flag[2] == 1) // sill
	{
	  sensmat[p] = 0;
	  p++; l++;
	}
     for(j = l; j < *npar; j++) // correlation parameters
	{
	  sensmat[p] = 0;
	  p++;
	}
   }
 // Derivative of the PLLik respect with the nugget
 if(flag[1] == 1)
   {
     l = i; i++;
     sensmat[p] = 2 * (EQ * (1 / c - 4 * pow(d, 2)) + (2 * d * ER - 1) / c) - 1 / pow(a,2);
     p++; l++;
     if(flag[2] == 1) // sill
	{
	  sensmat[p] = 2 * ( 2 * (ES * (3 * d - corr * e) - EM * (d * corr + e) -
			 d * (d - corr * e) * (2 * EQ + 1)) - 1 / c) - 1 / pow(a,2) ;
	  p++; l++;
	}
     h = 0;
     for(j = l; j < *npar; j++) // correlation parameters
	{
	  sensmat[p] = 4 * sill * (d * e * (2 * EQ + 1) - e * ES - d * EM) * gradcorr[h];
	  h++; p++;
	}
   }
 // Derivative of the PLLik respect with the sill
 if(flag[2] == 1)
   {
     l = i; i++;
     sensmat[p] = 2 * (- 2 * pow(q, 2) * (2 * EQ + 1) + 2 * (ES * (2 * q + d * t) -
						      EM * (2 * corr * q + e * t)) - t / c) + 1 / pow(a,2);
     p++; l++;
     h = 0;
     for(j = l; j < *npar; j++) // correlation parameters
	{
	  sensmat[p] = 2 * ((2 * e * (sill * (q * (2 * EQ + 1) - ES) - EQ + e * b) +
		     EM * (1 - 2 * sill * q))  * gradcorr[h]);
	  h++; p++;
	}
   }
 // Derivative of the PLLik respect with the correlation parameters
 D =  - pow(sill / c, 2) * (2 * pow(b,2) + c) ;


 for(j = 0; j < *nparc; j++)
   for(l = j; l < *nparc; l++)
     {
	sensmat[p] = 2 * D * gradcorr[j] * gradcorr[l];
	p++;
     }

 //---- END COMPUTATION OF THE SENSITIVITY MATRIX-----//
 return;
}

void Vari_SubSamp(double *coordx, double *coordy, int *corrmod, double *data, 
		  double *dista, double *eps, int *flagcorr, int *flagnuis, 
		  int *like, int *lonlat, int *npair, int *npar, int *nparc, 
		  int *nsite, double *nuisance, double *parcorr, int *type,
		  double *varimat, double *winc)
{
  double corr=0.0, lag=0.0, *gradcorr, *gradient, *rangex, *rangey;
  double *scoordx, *scoordy, *sdata, *sumgrad, *xgrid, *ygrid;
  double deltax=0.0, deltay=0.0, dimwinx=0.0, dimwiny=0.0;
  int *npts, nsubsam=0, numintx=0, numinty=0;
  int h=0, i=0, l=0, m=0, nwpair=0, p=0, q=0, j=0;

  gradcorr = (double *) R_alloc(*nparc, sizeof(double));
  gradient = (double *) R_alloc(*npar, sizeof(double));
  sumgrad = (double *) R_alloc(*npar, sizeof(double));

  npts = (int *) R_alloc(1, sizeof(int));

  rangex = (double *) R_alloc(2, sizeof(double));
  rangey = (double *) R_alloc(2, sizeof(double));

  scoordx = (double *) R_alloc(*nsite, sizeof(double));
  scoordy = (double *) R_alloc(*nsite, sizeof(double));
  sdata = (double *) R_alloc(*nsite, sizeof(double));

  Range(coordx, rangex, nsite);// range of the x-coordinate
  Range(coordy, rangey, nsite);// range of the y-coordinate

  // set the sub-sampling window based on prototype unit window (R_0)
  // and scaling factor (lambda_n)
  deltax = rangex[1] - rangex[0];// R_n = lambda_n * R_0
  deltay = rangey[1] - rangey[0];

  dimwinx = *winc * sqrt(deltax);// lambda*_n = constant * sqrt(lambda_n)
  dimwiny = *winc * sqrt(deltay);
  numintx = (int) deltax / dimwinx + 1;
  numinty = (int) deltay / dimwiny + 1;

  xgrid = (double *) R_alloc(numintx, sizeof(double));
  ygrid = (double *) R_alloc(numinty, sizeof(double));

  Seq(rangex, numintx, xgrid);
  Seq(rangey, numinty, ygrid);

  for(i = 0; i < numintx; i++)
    for(j = 0; j < numinty; j++)
      {
	*npts = 0;
	SetSampling(coordx, coordy, data, npts, scoordx, scoordy, sdata, nsite,
		    xgrid[i + 1], xgrid[i], ygrid[j + 1], ygrid[j]);
	if(*npts > 2)
	  {
	    for(h = 0; h < *npar; h++)
	      sumgrad[h] = 0;

	    nwpair = *npts * (*npts - 1) / 2;
	    for(l = 0; l < (*npts - 1); l++)
	      for(m = (l + 1); m < *npts; m++)
		{
		  if(*lonlat)
		    lag = Dist_geodesic(scoordx[l], scoordy[l], scoordx[m], scoordy[m]);
		  else 
		    lag = pythag(scoordx[l] - scoordx[m], scoordy[l] - scoordy[m]);
		  if(lag <= *dista)
		    {
		      corr = CorrelationFct(corrmod, lag, parcorr);
		      GradientCorrFct(corr, corrmod, eps, flagcorr, gradcorr, lag, parcorr);
		      switch(*like)
			{
			case 1:// Conditional likelihood:
			  switch(*type)
			    {
			    case 2:
			      Grad_Cond_Gauss(corr, flagnuis, gradcorr, gradient, 
					      npar, nuisance, sdata[l], sdata[m]);
			      break;
			    }
			  break;
			case 3: // Marginal likelihood:
			  switch(*type)
			    {
			    case 1: // Gradient of the log difference likelihood
			      Grad_Diff_Gauss(corr, flagnuis, gradcorr, gradient, 
					      npar, nuisance, sdata[l], sdata[m]);
			      break;
			    case 2: // Gradient of the log pairwise likelihood
			      Grad_Pair_Gauss(corr, flagnuis, gradcorr, gradient, 
					      npar, nuisance, sdata[l], sdata[m]);
			      break;
			    }
			  break;
			}
		    }
		  for(h = 0; h < *npar; h++)
		    sumgrad[h] = sumgrad[h] + gradient[h];
		}
	    h = 0;
	    for(p = 0; p < *npar; p++)
	      for(q = p; q < *npar; q++)
		{
		  varimat[h] = varimat[h] + sumgrad[p] * sumgrad[q] / nwpair;
		  h++;
		}
	    nsubsam++;
	  }
      }
  h = 0;
  for(p = 0; p < *npar; p++)
    for(q = p; q < *npar; q++)
      {
	varimat[h] = *npair * varimat[h] / nsubsam;
	h++;
      }

  return;
}

/*
void Vari_SubSamp(double *coordx, double *coordy, int *corrmod, double *data, 
		  int *dimwin, double *dista, double *eps, int *flagcorr, 
		  int *flagnuis, int *like, int *lonlat, int *npair, int *npar, 
		  int *nparc, int *nsite, double *nuisance, int *numwin, 
		  double *parcorr, double *varimat, int *type)
{

  double corr=0.0, lag=0.0, *gradcorr, *gradient, *rangex, *rangey;
  double *scoordx, *scoordy, *sdata, *sumgrad, *xgrid, *ygrid;
  int *npts, nsubsam=0, numint=*numwin+1, nsteps=*numwin-*dimwin;
  int h=0, i=0, l=0, m=0, nwpair=0, p=0, q=0, j=0;

  gradcorr = (double *) R_alloc(*nparc, sizeof(double));
  gradient = (double *) R_alloc(*npar, sizeof(double));
  sumgrad = (double *) R_alloc(*npar, sizeof(double));

  npts = (int *) R_alloc(1, sizeof(int));

  rangex = (double *) R_alloc(2, sizeof(double));
  rangey = (double *) R_alloc(2, sizeof(double));

  scoordx = (double *) R_alloc(*nsite, sizeof(double));
  scoordy = (double *) R_alloc(*nsite, sizeof(double));
  sdata = (double *) R_alloc(*nsite, sizeof(double));

  xgrid = (double *) R_alloc(numint, sizeof(double));
  ygrid = (double *) R_alloc(numint, sizeof(double));

  Range(coordx, rangex, nsite);
  Range(coordy, rangey, nsite);

  Seq(rangex, numint, xgrid);
  Seq(rangey, numint, ygrid);

  for(i = 0; i < nsteps; i++)
    for(j = 0; j < nsteps; j++)
      {
	*npts = 0;
	SetSampling(coordx, coordy, data, npts, scoordx, scoordy, sdata, nsite,
		    xgrid[i + *dimwin], xgrid[i], ygrid[i + *dimwin], ygrid[i]);
	if(*npts > 2)
	  {
	    for(h = 0; h < *npar; h++)
	      sumgrad[h] = 0;

	    nwpair = *npts * (*npts - 1) / 2;
	    for(l = 0; l < (*npts - 1); l++)
	      for(m = (l + 1); m < *npts; m++)
		{
		  if(*lonlat)
		    lag = Dist_geodesic(scoordx[l], scoordy[l], scoordx[m], scoordy[m]);
		  else 
		    lag = pythag(scoordx[l] - scoordx[m], scoordy[l] - scoordy[m]);
		  if(lag <= *dista)
		    {
		      corr = CorrelationFct(corrmod, lag, parcorr);
		      GradientCorrFct(corr, corrmod, eps, flagcorr, gradcorr, lag, parcorr);
		      switch(*like)
			{
			case 1:// Conditional likelihood:
			  switch(*type)
			    {
			    case 2:
			      Grad_Cond_Gauss(corr, flagnuis, gradcorr, gradient, 
					      npar, nuisance, sdata[l], sdata[m]);
			      break;
			    }
			  break;
			case 3: // Marginal likelihood:
			  switch(*type)
			    {
			    case 1: // Gradient of the log difference likelihood
			      Grad_Diff_Gauss(corr, flagnuis, gradcorr, gradient, 
					      npar, nuisance, sdata[l], sdata[m]);
			      break;
			    case 2: // Gradient of the log pairwise likelihood
			      Grad_Pair_Gauss(corr, flagnuis, gradcorr, gradient, 
					      npar, nuisance, sdata[l], sdata[m]);
			      break;
			    }
			  break;
			}
		    }
		  for(h = 0; h < *npar; h++)
		    sumgrad[h] = sumgrad[h] + gradient[h];
		}
	    h = 0;
	    for(p = 0; p < *npar; p++)
	      for(q = p; q < *npar; q++)
		{
		  varimat[h] = varimat[h] + sumgrad[p] * sumgrad[q] / nwpair;
		  h++;
		}
	    nsubsam++;
	  }
      }
  h = 0;
  for(p = 0; p < *npar; p++)
    for(q = p; q < *npar; q++)
      {
	varimat[h] = *npair * varimat[h] / nsubsam;
	h++;
      }

  return;
}
*/
