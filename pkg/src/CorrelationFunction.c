#include "header.h"

// list of correlation function:
double CorrelationFct(int *corrmod, double lag, double *par)
{
  double corr=0.0, power=0.0, power1=0.0;
  double power2=0.0, scale=0.0, smooth=0.0;

  switch(*corrmod) // Correlation functions are in alphabetical order
    {
    case 1:// Cauchy correlation function
      power2 = par[0];
      scale = par[1];
      corr = pow((1 + pow(lag / scale, 2)), - power2);
      break;
    case 2:// Exponential correlation function
      scale = par[0];
      corr = exp(- lag / scale);
      break;
    case 3:// Gaussian correlation function
      scale = par[0];
      corr = exp(-pow(lag / scale, 2));
      break;
    case 4: // Generalised Cuachy correlation function 
      power1 = par[0];
      power2 = par[1];
      scale = par[2];
      corr = pow((1 + pow(lag / scale, power1)), - power2 / power1);
      break;
    case 5:// Stable correlation function
      power = par[0];
      scale = par[1];
      corr = exp(-pow(lag / scale, power));
      break;
    case 6://  Whittle-Matern correlation function
      scale = par[0];
      smooth = par[1];
      corr = pow(2, 1 - smooth) / gammafn(smooth) * pow(lag / scale, smooth) * 
	bessel_k(lag / scale, smooth, 1);
      break;
    }

  return corr;
}

// Computation of the lower (upper) triangular correlation matrix
void CorrelationMat(double *corr, int *corrmod, double *lags, int *npairs, double *par)
{
  int i;

  for(i = 0; i < *npairs; i++)
    corr[i] = CorrelationFct(corrmod, lags[i], par);

  return;
}

// list of the derivatives (respect with the parameters) of the correlation functions:
void GradientCorrFct(double corr, int *corrmod, double *eps, int *flag, double *grad, 
		     double lag, double *par)
{
  int i=0;
  double power=0.0, power1=0.0, power2=0.0, scale=0.0, smooth=0.0;
  double parscale=0.0, parsmooth=0.0;

  switch(*corrmod)// Correlation functions are in alphabetical order
    {
    case 1:// Cauchy correlation function
      power2 = par[0];
      scale = par[1];
      if(flag[0] == 1)
	{
	  grad[i] = - corr * log(pow(corr, - 1 / power2));
	  i++;
	}
      if(flag[1] == 1)
	{
	  grad[i] = 2 * power2 * corr * pow(corr, 1 / power2) * 
	    pow(lag, 2) / pow(scale, 3);
	}
      break;
    case 2:// Exponential correlation function
      if(flag[0] == 1)
	{
	  scale = par[0];
	  grad[i] = corr * lag / pow(scale, 2);
	  i++;
	}
      break;
    case 3:// Gaussian correlation function
      if(flag[0] == 1)
	{
	  scale = par[0];
	  grad[i] = 2 * corr * pow(lag, 2) / pow(scale, 3);
	}
	  break;
    case 4:// Generalised Cuachy correlation function
      power1 = par[0];
      power2 = par[1];
      scale = par[2];
      if(flag[0] == 1)
	{
          grad[i] = power2 * corr / power1 * (log(1 + pow(lag / scale, power1)) / 
					      power1 - pow(lag / scale, power1) * 
					      log(lag / scale) / (1 + pow(lag/ scale, power1)));
	  i++;
	}
      if(flag[1] == 1)
	{
	  grad[i] = - corr * log(1 + pow(lag / scale, power1)) / power1;
	  i++;
	}
      if(flag[2] == 1)
	{
	  grad[i] = corr / (1 + pow(lag / scale, 2)) * power2 * 
	    pow(lag, power1) / pow(scale, power1 + 1);
	}
      break;
    case 5:// Stable correlation function
      power = par[0];
      scale = par[1];
      if(flag[0] == 1)
	{
	  grad[i] = - corr * pow(lag / scale, power) * log(lag / scale);
	  i++;
	}
      if(flag[1] == 1)
	{
	  grad[i] = corr * pow(lag / scale, power - 1) * 
	    power * lag / pow(scale, 2);
	}
      break;
    case 6:// Whittle-Matern correlation function
      scale = par[0];
      smooth = par[1];
      if(flag[0] == 1)
	{
	   parscale = (bessel_k(lag / (scale + *eps), smooth, 1) - 
		       bessel_k(lag / scale, smooth, 1)) / *eps;
	   grad[i] = pow(2, 1 - smooth) / gammafn(smooth) * 
	     pow(lag / scale, smooth) * (parscale - smooth * 
					 bessel_k(lag / scale, smooth, 1) / scale);

	   i++;
	}
      if(flag[1] == 1)
	{
	  parsmooth = (bessel_k(lag / scale, smooth + *eps, 1) - 
		       bessel_k(lag / scale, smooth, 1)) / *eps;
	      
	  grad[i] = pow(2, 1 - smooth) * pow(lag / scale, smooth) / gammafn(smooth) * 
	    (log(lag / scale) - log(2) - digamma(smooth) * 
	     bessel_k(lag / scale, smooth, 1) + parsmooth);
	}
      break;
    }
 
  return;
}

