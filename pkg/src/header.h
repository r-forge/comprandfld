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

double Variogram(int *corrmod, double lag, double *nuisance, double *par);

void VectCorrelation(int *corrmod, double *corrfun, double *lags, int *npair, double *par);

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
		    double *dista, double *lags, int *like, int *model, double *nuisance, 
		    int *ndata, int *nsite, double *par, double *res, int *type);

double PairLikelihood(double corr, int *like, double *nuisance, 
		      double s1, double s1s, double u, double v, int *type);


/*----------------------------------------------------------------
File name: CompositeLikelihood.c
Description: functions for composite log-likelihood evaluation
End
 ---------------------------------------------------------------*/

// 3)
/*----------------------------------------------------------------
File name: Godambe.c
Description: procedures for the computation of the Godambe matrix
Start
 ---------------------------------------------------------------*/

void CompScore(double *coordx, double *coordy, int *corrmod, double *data, 
	       double *eps, int *flag, int *flagcorr, int *model, int *like,
               int *ndata, int *ngrc, int *npar, int *nsite, double *par, 
               double *parcorr, double *res, int *type, int *weight);

void GodambeMat_emp(double *coordx, double *coordy, int *corrmod, double *data, 
		    double *eps, int *flagcorr, int *flagnuis, int *like, int *model, 
                    int *ndata, int *npar, int *nparc, int *nsite, double *parcorr, 
                    double *nuisance, double *godambe, double *score, int *type);

void GodambeMat(double *coordx, double *coordy, int *corrmod, double *data, double *dista, 
		double *eps, int *flagcorr, int *flagnuis, double *lags, int *like, 
		int *lonlat, int *model, int *npar, int *nparc, int *nsite, double *parcorr, 
		double *nuisance, double *sensmat, int *type, double *varimat, int *vartype, 
		double *winc);

void GodambeMat_Diff(double *coordx, double *coordy, int *corrmod, double *dista, double *eps, 
		     int *flagcorr, int *flagnuis, double *lags, int *model, int *npar, int *nparc, 
		     int *nsite, double *parcorr, double *nuisance, double *sensmat, double *varimat);

void Grad_Hess_Gauss_Pair(double corr, int *flag, double *gradcorr, double *gradient, 
			  double *hessian, double *hesscorr, int *npar,  double *par, 
			  double u, double v);

void Grad_Cond_Gauss(double corr, int *flag, double *gradcorr, double *gradient, int *npar, 
		     double *par, double u, double v);

void Grad_Diff_Gauss(double corr, int *flag, double *gradcorr, 
		     double *gradient, int *npar, double *par, double u, double v);

void Grad_Diff_Vario(double corr, int *flag, double *gradcorr, 
		     double *gradient, int *npar, double *par);

void Grad_Pair_Gauss(double corr, int *flag, double *gradcorr, double *gradient, int *npar,
		     double *par, double u, double v);

void Sens_Cond_Gauss(int *corrmod, double *dista, double *eps, int *flagcorr,
		     int *flagnuis, double *lags, int *nsite, double *nuisance,
		     int *npair, int *npar, int *nparc, double *parcorr, double *sensmat);

void Sens_Cond_Gauss_ij(double corr, int *flag, double *gradcorr, int *npar, 
			int *nparc, double *par, double *sensmat);

void Sens_Diff_Gauss(int *corrmod, double *dista, double *eps, int *flagcorr,
		     int *flagnuis, double *lags, int *nsite, double *nuisance,
		     int *npair, int *npar, int *nparc, double *parcorr, double *sensmat);

void Sens_Diff_Gauss_ij(double *gradient, int *npar, double *sensmat);

void Sens_Pair_Gauss(int *corrmod, double *dista, double *eps, int *flagcorr,
		     int *flagnuis, double *lags, int *nsite, double *nuisance,
		     int *npair, int *npar, int *nparc, double *parcorr, double *sensmat);

void Sens_Pair_Gauss_ij(double corr, int *flag, double *gradcorr, int *npar, 
			int *nparc, double *par, double *sensmat);

void Sensitivity(double *coordx, double *coordy, int *corrmod, double *dista, 
		 double *eps, int *flagcorr, int *flagnuis, double *lags, 
		 int *like, int *model, int *npair, int *npar, int *nparc, 
		 int *nsite, double *parcorr, double *nuisance, double *sensmat, 
		 int *type);

void Vari_SubSamp(double *coordx, double *coordy, int *corrmod, double *data, 
		  double *dista, double *eps, int *flagcorr, int *flagnuis, 
		  int *like, int *lonlat, int *npair, int *npar, int *nparc, 
		  int *nsite, double *nuisance, double *parcorr, int *type,
		  double *varimat, double *winc);

/*----------------------------------------------------------------
File name: Godambe.c
Description: procedures for the computation of the Godambe matrix
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

void Empiric_Variogram(double *bins, double *coordx, double *coordy, double *data, double *lags, 
		       double *lenbins, double *maxdist, double *moments, int *npairs, int *nsite, 
		       int *nbins);

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

void Distances(double *coordx, double *coordy, double *lags, int *nsite, int *type);

double Dist_geodesic(double lonx, double latx, double lony, double laty);

void SetSampling(double *coordx, double *coordy, double *data, int *npts, 
		 double *scoordx, double *scoordy, double *sdata, int *size, 
		 double xmax, double xmin, double ymax, double ymin);

double Maxima(double *x, int *size);

double Minima(double *x, int *size);

void Range(double *x, double *ran, int *size);

void Seq(double *x, int len, double *res);

/*----------------------------------------------------------------
File name: Utility.c
Description: procedures for the computation of useful quantities.
End
 ---------------------------------------------------------------*/
