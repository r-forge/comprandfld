#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>

#define LOW -1.0e15


// 1)
/*----------------------------------------------------------------
File name: CorrelationFunction.c
Description: procedures for computation of correlation functions
Start
 ---------------------------------------------------------------*/

double CorrelationFct(int *corrmod, double lag, double *par);

void CorrelationMat(double *corr, int *corrmod, double *lags, 
		    int *npairs, double *par);

void GradientCorrFct(double corr, int *corrmod, double *eps, int *flag, 
		     double *grad,  double lag, double *par);

/*----------------------------------------------------------------
File name: CorrelationFunction.c
Description: procedures for computation of correlation functions
End
 ---------------------------------------------------------------*/

// 2)
/*----------------------------------------------------------------
File name: CompositeLikelihood.c
Description: functions for composite log-likelihood evaluation
Start
 ---------------------------------------------------------------*/


void CompLikelihood(double *coordx, double *coordy, int *corrmod, double *data, 
		    int *model, double *nuisance, int *ndata, int *nsite, 
		    double *par, double *res, int *type);

double PairLikelihood(double corr, double *nuisance, double s1, double s1s, 
		      double u, double v, int *type);


/*----------------------------------------------------------------
File name: CompositeLikelihood.c
Description: functions for composite log-likelihood evaluation
End
 ---------------------------------------------------------------*/

// 3)
/*----------------------------------------------------------------
File name: Gradient.c
Description: procedures for computation of gradients quantities
Start
 ---------------------------------------------------------------*/

void CompScore(double *coordx, double *coordy, int *corrmod, double *data, 
	       double *eps, int *flag, int *flagcorr, int *model, int *ndata, 
	       int *ngrc, int *npar, int *nsite, double *par, double *parcorr, 
	       double *res, int *type, int *weight);

void GodambeMat(double *coordx, double *coordy, int *corrmod, double *data, 
		double *eps, int *flagcorr, int *flagnuis, int *model, int *ndata, 
		int *npar, int *nparc, int *nsite, double *parcorr, double *nuisance,
		double *godambe, int *type);

void Score_Gauss_Diff(double corr, int *flag, double *gradcorr, double *gradient, int *npar,
		      double *par, double u, double v);

void Score_Gauss_Pair(double corr, int *flag, double *gradcorr, double *gradient, int *npar,
		      double *par, double u, double v);


/*----------------------------------------------------------------
File name: Gradient.c
Description: procedures for computation of gradients quantities
End
 ---------------------------------------------------------------*/

// 4)
/*----------------------------------------------------------------
File name: Likelihood.c
Description: procedures for computation of the full likelihood
Start
 ---------------------------------------------------------------*/

// Insert here

/*----------------------------------------------------------------
File name: Likelihood.c
Description: procedures for computation of the full likelihood
End
 ---------------------------------------------------------------*/


// 5)
/*----------------------------------------------------------------
File name: WeightedLeastSquare.c
Description: procedures for the estimation of the model parameters
via the weighted least square method.
Start
 ---------------------------------------------------------------*/

void Empiric_Variogram(double *bins, double *coordx, double *coordy, double *data, 
		       double *lenbins, double *maxdist, double *moments, int *nsite, int *nbins);

void Wls(double *bins, int *corrmod, double *par, int *nbins, double *moments, 
	 double *lenbins, double *nuisance, int *weighted, double *res);

/*----------------------------------------------------------------
File name: WeightedLeastSquare.c
Description: procedures for the estimation of the model parameters
via the weighted least square method.
End
 ---------------------------------------------------------------*/

// 6)
/*----------------------------------------------------------------
File name: Utility.c
Description: procedures for the computation of useful quantities.
Start
 ---------------------------------------------------------------*/

void Distances(double *coordx, double *coordy, int *nsite, double *distances);

double Maxima(double *x, int size);

double Minima(double *x, int size);

/*----------------------------------------------------------------
File name: Utility.c
Description: procedures for the computation of useful quantities.
End
 ---------------------------------------------------------------*/
