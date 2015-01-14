#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>

#define LOW -1.0e15
#define MAXERR 1e-6

//---------START GLOBAL VARIABLES-----------
int *isst;//is a spatio-temporal random field?
int *istap;//is tapering?
double *lags;// vector of spatial distances for tapering
double *lagt;// vector of temporal distance for tapering
double **mlags;// vector of spatial distances
double **mlagt;// vector of temporal distances
double *maxdist;// the threshould of the spatial distances
double *maxtime;// the threshould of the temporal distances
// below which the pairs are considered
double *maximdista;// the maximum spatial distance
double *maximtime;// the maximum temporal distance
double *minimdista; // the minimum spatial distance
double *minimtime;// the minimum temporal distance
int *ncoord;// number of total spatial coordinates
int *ncoordx;// number of the first spatial coordinates
int *ncoordy;// number of the second spatial coordinates
int *npairs;// effective number of pairs
int *nrep;// number of iid replicates of the random field
int *ntime;// number of times
double *tapsep; // parameter separability for space time quasi taper
//int *totpairs;// total number of pairs
//---------END GLOBAL VARIABLES-------------
//void indx(double *ind, int *n);

// fortran declaration for bivariate normal cdf:
double F77_NAME(bvnmvn)(double *lower, double *upper, int* infin, double *correl);

// Internal function declaration:
// 1)
/*----------------------------------------------------------------
File name: CorrelationFunction.c
Description: procedures for computation of correlation functions
Start
---------------------------------------------------------------*/

double CheckCor(int *cormod, double *par);

void Comp_supp(double *c_supp,int *cormod, double h,double u, double *par);

double CorFct(int *cormod, double h, double u, double *par);

double CorFunCauchy(double lag, double power2, double scale);

double CorFunGenCauchy(double lag, double power1, double power2, double scale);

double CorFunSferical(double lag, double scale);

double CorFunStable(double lag, double power, double scale);

		       //double CorFunDobStable(double lag, double power_s, double power_t,
		       //double scale_s, double scale_t, double tsep);

double CorFunWitMat(double lag, double scale, double smooth);

double CorFunWend1(double lag,double scale);

double CorFunWend2(double lag,double scale);

double CorFunWend3(double lag,double scale);

void CorrelationMat(double *rho, int *cormod, double *nuis, double *par);

void CorrelationMat_tap(double *rho, int *cormod, double *nuis, double *par);

void CorrelationMat_st(double *rho, int *cormod, double *nuis, double *par);

void Corr_c(double *cc,double *coordx, double *coordy, double *coordt, int *cormod, int *grid, double *locx,  double *locy,int *ncoord, int *nloc,int*tloc,
                 int *ntime, double *par, int *spt, double *time,int *type);

void Corr_c_tap(double *cc,double *cc_tap,double *coordx, double *coordy, double *coordt, int *cormod, int *cormodtap, int *grid, double *locx,  double *locy,int *ncoord, int *nloc,int*tloc,
                 int *ntime, double *par, int *spt, double *time,int *type);

void DCorrelationMat(int *cormod,double *drho,double *eps,int *flagcor,
		     int *nparcor, double *parcor,double *rho);

void DCorrelationMat_tap(int *cormod,double *drho,double *eps,int *flagcor,
		     int *nparcor, double *parcor,double *rho);

void DCorrelationMat_st(int *cormod,double *drho,double *eps,int *flagcor,
			int *nparcor, double *parcor,double *rho);

double DCauchyPow(double power2, double scale, double rho);

double DCauchySc(double lag, double power2, double scale, double rho);

double DExpoSc(double lag, double scale, double rho);

double DGaussSc(double lag, double scale, double rho);

double DGenCauP1(double lag, double power1, double power2, double scale, double rho);

double DGenCauP2(double lag, double power1, double scale, double rho);

double DGenCauSc(double lag, double power1, double power2, double scale, double rho);

double DSferiSc(double lag, double scale);

double DStabPow(double lag, double power, double scale, double rho);

double DStabSc(double lag, double power, double scale, double rho);

double DWhMatSc(double *eps, double lag, double scale, double smooth);

double DWhMatSm(double *eps, double lag, double scale, double smooth);

double DGneiting_sc_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);

double DGneiting_sc_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);

double DGneiting_pw_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);

double DGneiting_pw_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);

double DGneiting_sep(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);

double DIaco_sc_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double power2);

double DIaco_sc_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double power2);

double DIaco_pw_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double power2);

double DIaco_pw_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double power2);

double DIaco_pw2(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double power2);

double DPorcu_sc_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);

double DPorcu_sc_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);

double DPorcu_pw_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);

double DPorcu_pw_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);

double DPorcu_sep(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);

double DStein_sc_s(double h,double u, double power_t,double scale_s,double scale_t,double smooth);

double DStein_sc_t(double h,double u, double power_t,double scale_s,double scale_t,double smooth);

double DStein_pw_s(double h,double u, double power_t,double scale_s,double scale_t,double smooth);

double DStein_pw_t(double h,double u, double power_t,double scale_s,double scale_t,double smooth);

double DStein_sm(double h,double u,   double power_t,double scale_s,double scale_t,double smooth);

double DExp_Cauchy_sc_s(double h,double u,double scale_s,double scale_t,double power2);

double DExp_Cauchy_sc_t(double h,double u,double scale_s,double scale_t,double power2);

double DExp_Cauchy_pw2 (double h,double u,double scale_s,double scale_t,double power2);

double DExp_Exp_sc_s(double h,double u,double scale_s,double scale_t);

double DExp_Exp_sc_t(double h,double u,double scale_s,double scale_t);

double DMat_Exp_sc_s(double h,double u,double scale_s,double scale_t,double smooth);

double DMat_Exp_sc_t(double h,double u,double scale_s,double scale_t,double smooth);

double DMat_Exp_sm(double h,double u,double *eps,double scale_s,double scale_t,double smooth);

double DMat_Cauchy_sc_s(double h,double u,double power2,double scale_s,double scale_t,double smooth);

double DMat_Cauchy_sc_t(double h, double u,double power2,double scale_s,double scale_t,double smooth);

double DMat_Cauchy_pw2(double h,double u,double power2,double scale_s,double scale_t,double smooth);

double DMat_Cauchy_sm(double h,double u,double *eps, double power2,double scale_s,double scale_t,double smooth);

void ExtCoeff(int *cormod, double *extc, double *lag, int *model,
	      int *nlags, double *nuis, double *par);

void GradCorrFct(double rho, int *cormod, double *eps, int *flag,
		 double *grad, double h, double u, double *par);

void GradVarioFct(double vario, int *cormod, double *eps, int *flag,
		  double *grad, double lag, double *par, double tsep);

void TapVectCorrelation(double *rho,int *cormod,double *tdists,int *ntdists,double *nuis,double *par);

double Variogram(int *cormod, double h, double u, double *nuis, double *par);

double VarioFct(int *cormod, double h, double *par, double u);

double VarioDobStable(double lag, double power_s, double power_t, double scale_s, double scale_t, double tsep);

double VarioGCauchy(double lag, double power1, double power2, double scale);

void VectCorrelation(double *rho, int *cormod, double *h, int *nlags,
		     int *nlagt, double *par, double *u);

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

double BrowResnllik(double a, double c, double x, double y);

void Comp_Cond_Gauss(int *cormod, double *data, double *nuis, double *par, double *thr, double *res);

void Comp_Cond_Gauss_st(int *cormod, double *data, double *nuis, double *par, double *thr, double *res);

void Comp_Cond_BinGauss(int *cormod, double *data, double *nuis, double *par, double *thr, double *res);

void Comp_Cond_BinGauss_st( int *cormod, double *data, double *nuis, double *par,double *thr, double *res);

void Comp_Diff_Gauss(int *cormod, double *data, double *nuis, double *par, double *thr, double *res);

void Comp_Diff_Gauss_st(int *cormod, double *data, double *nuis, double *par, double *thr, double *res);

void Comp_Diff_BinGauss( int *cormod, double *data, double *nuis, double *par, double *thr,double *res);

void Comp_Diff_BinGauss_st(int *cormod, double *data, double *nuis, double *par,double *thr, double *res);

void Comp_Pair_Gauss(int *cormod, double *data, double *nuis, double *par, double *thr, double *res);

void Comp_Pair_Gauss_st(int *cormod, double *data, double *nuis, double *par, double *thr, double *res);

void Comp_Pair_BinGauss(int *cormod, double *data, double *nuis, double *par, double *thr,double *res);

void Comp_Pair_BinGauss_st(int *cormod, double *data, double *nuis, double *par, double *thr,double *res);

void Comp_Brow_Resn(int *cormod, double *data, double *nuis, double *par, double *thr, double *res);

void Comp_Brow_Resn_st(int *cormod, double *data, double *nuis, double *par, double *thr, double *res);

void Comp_Ext_Gauss(int *cormod, double *data, double *nuis, double *par, double *thr, double *res);

void Comp_Ext_T(int *cormod, double *data, double *nuis, double *par, double *thr, double *res);


/*----------------------------------------------------------------
File name: CompositeLikelihood.c
Description: functions for composite log-likelihood evaluation
End
 ---------------------------------------------------------------*/

// 3)
/*----------------------------------------------------------------
File name: Distributions.c
Description: procedures for the computation of useful distributions
Start
 ---------------------------------------------------------------*/

double d1x_dt(double x, double df);

double d2x_dt(double x, double df);

double ddf_pt(double x, double df);

double ddf_t_dt(double x, double clc, double df, double somc2);

double ddf_t_d1x_dt(double x, double clc, double df, double somc2);

double dts(double x, double df);

void Dist2Dist(double *data, double *eloc, double *escale, double *eshape,
	       int *ndata, int *nsite, double *ploc, double *pscale,
	       double *pshape, int *type, double *res);

double dgev(double x, double loc, double scale, double shape);

double d2norm(double x, double y, double rho);

void GevLogLik(double *data, int *ndata, double *par, double *res);

double int_pt(double x, double df);

void integr_pt(double *x, int n, void *ex);

double pbnorm(int *cormod, double h, double u, double *nuis, double *par, double thr);

double pgev(double x, double loc, double scale, double shape);

double qgev(double x, double loc, double scale, double shape);

void vpbnorm(int *cormod, double *h, double *u, int *nlags, int *nlagt,
	     double *nuis, double *par, double *rho, double *thr);

/*----------------------------------------------------------------
File name: Distributions.c
Description: procedures for the computation of useful distributions
End
 ---------------------------------------------------------------*/

// 4)
/*----------------------------------------------------------------
File name: Godambe.c
Description: procedures for the computation of the Godambe matrix
Start
 ---------------------------------------------------------------*/

void God_Cond_Gauss(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
		    int *npar, int *nparc, double *parcor, double *nuis, double *score,
		    double *sensmat, double *varimat);

void God_Cond_Gauss_st(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
		       int *npar, int *nparc, double *parcor, double *nuis, double *score,
		       double *sensmat, double *varimat);

void God_Diff_Gauss(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
		    int *npar, int *nparc, double *parcor, double *nuis, double *score,
		    double *sensmat, double *varimat);

void God_Diff_Gauss_st(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
		       int *npar, int *nparc, double *parcor, double *nuis, double *score,
		       double *sensmat, double *varimat);

void God_Pair_Gauss(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
		    int *npar, int *nparc, double *parcor, double *nuis, double *score,
		    double *sensmat, double *varimat);

void God_Pair_Gauss_st(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
		       int *npar, int *nparc, double *parcor, double *nuis, double *score,
		       double *sensmat, double *varimat);

void God_BrowResn(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
		   int *npar, int *nparc, double *parcor, double *nuis, double *score,
		   double *sensmat, double *varimat);

void God_BrowResn_st(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
		     int *npar, int *nparc, double *parcor, double *nuis, double *score,
		     double *sensmat, double *varimat);

void God_Ext_Gauss(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
		   int *npar, int *nparc, double *parcor, double *nuis, double *score,
		   double *sensmat, double *varimat);

void God_Ext_T(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
	       int *npar, int *nparc, double *parcor, double *nuis, double *score,
	       double *sensmat, double *varimat);

void GodambeMat_emp(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
		    int *like, int *model, int *npar, int *nparc, double *parcor, double *nuis,
		    double *score, double *sensmat, int *spt, double *varimat, int *type);

void GodambeMat(double *coordx, double *coordy, int *cormod, double *data, int *dist,
		double *eps,int *flagcor, int *flagnuis, int *grid, int *like, int *model,
		int *npar, int *nparc, double *parcor, double *nuis, double *score,
		double *sensmat, int *spt, double *thr, int *type, double *varimat,
		int *vartype, double *winc, double *winstp);

void GodambeMat_Diff(double *coordx, double *coordy, int *cormod, double *eps, int *flagcor,
		     int *flagnuis, int *model, int *npar, int *nparc, double *parcor,
		     double *nuis, double *sensmat, double *varimat);

void Sens_Cond_Gauss(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis, double *nuis, int *np,
		     int *npar, int *nparc, double *parcor, double *score, double *sensmat);

void Sens_Cond_Gauss_st(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis, double *nuis, int *np,
			int *npar, int *nparc, double *parcor, double *score, double *sensmat);

void Sens_Cond_Gauss_ij(double rho, int *flag, double *gradcor, int *npar,
			int *nparc, double *par, double *sensmat);

void Sens_Diff_Gauss(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
		     double *nuis, int *np,int *npar, int *nparc, double *parcor,double *score, double *sensmat);

void Sens_Diff_Gauss_st(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
			double *nuis, int *np,int *npar, int *nparc, double *parcor, double *score, double *sensmat);

void Sens_Diff_Gauss_ij(double *gradient, int *npar, double *sensmat);

void Sens_Pair_Gauss(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
		     double *nuis, int *np, int *npar, int *nparc, double *parcor, double *score, double *sensmat);

void Sens_Pair_Gauss_st(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis, double *nuis,
			int *np, int *npar, int *nparc, double *parcor, double *score, double *sensmat);

void Sens_Pair_Gauss_ij(double rho, int *flag, double *gradcor, int *npar,
			int *nparc, double *par, double *sensmat);

void Sensitivity(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
		 int *like, int *model, int *npar, int *nparc, double *parcor,
		 double *nuis, int *np, double *score, double *sensmat, int *spt, int *type);

void Vari_SubSamp(double *coordx, double *coordy, int *cormod, double *data,
		  int *dist, double *eps, int *flagcor, int *flagnuis, int *grid,
            int *like,int *model, int *npar, int *nparc, double *nuis, int *np,
		  double *parcor, double *thr, int *type, double *varimat,
		  double *winc, double *winstp);

void Vari_SubSamp_st(int *cormod, double *data, int *dist, double *eps, int *flagcor,int *flagnuis,
		     int *like, int *npar, int *nparc, double *nuis, int *np,double *parcor,
		     int *type, double *varimat, double *winc, double *winstp);

/*----------------------------------------------------------------
File name: Godambe.c
Description: procedures for the computation of the Godambe matrix
End
 ---------------------------------------------------------------*/

// 5)
/*----------------------------------------------------------------
File name: Gradient.c
Description: procedures for computation of the gredients
Start
 ---------------------------------------------------------------*/

void Grad_Cond_Bin(double rho,double pij, double p,int *flag, double *gradcor,
		   double *grad, int *npar, double *nuis, double *thr, double u, double v);

void Grad_Cond_Gauss(double rho, int *flag, double *gradcor, double *grad,
		     int *npar, double *par, double u, double v);

void Grad_Diff_Bin(double rho,double pij, double p,int *flag, double *gradcor,
		   double *grad, int *npar, double *nuis, double *thr, double u, double v);

void Grad_Diff_Gauss(double rho, int *flag, double *gradcor, double *grad,
		     int *npar, double *par, double u, double v);

void Grad_Diff_Vario(double rho, int *flag, double *gradcor, double *grad,
		     int *npar, double *par);

void Grad_Pair_Bin(double rho,double pij, double p,int *flag, double *gradcor,
		   double *grad, int *npar, double *nuis, double *thr, double u, double v);

void Grad_Pair_Gauss(double rho, int *flag, double *gradcor, double *grad,
		     int *npar, double *par, double u, double v);

void Grad_Ext_Gauss(double rho, int *flag, double *gradcor, double *grad,
		    int *npar, double *par, double u, double v);

void Grad_Brow_Resn(double vario, int *flag, double *gradcor, double *grad,
		    int *npar, double *par, double x, double y);

void Grad_Ext_T(double rho, int *flag, double *gradcor, double *grad,
		int *npar, double *par, double x, double y);

double mij(double qij, double w, double pij, double p);

double nij(double qij, double w, double pij, double p);

/*----------------------------------------------------------------
File name: Gradient.c
Description: procedures for computation of the gradients
End
 ---------------------------------------------------------------*/

// 5)
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


// 6)
/*----------------------------------------------------------------
File name: WeightedLeastSquare.c
Description: procedures for the estimation of the model parameters
via the weighted least square method.
Start
 ---------------------------------------------------------------*/
void Binned_Lorelogram(double *bins, double *data, int *lbins, double *moms,int *nbins);

void Binned_Lorelogram_st(double *bins, double *bint, double *data, int *lbins,
			  int *lbinst, int *lbint, double *moms,double *momst,
			  double *momt, int *nbins, int *nbint);

void Binned_Madogram(double *bins, double *data, int *lbins, double *moments, int *nbins);

void Binned_Madogram_st(double *bins, double *bint, double *data, int *lbins,
			int *lbinst, int *lbint, double *moms,double *momst,
			double *momt, int *nbins, int *nbint);

void Binned_Variogram(double *bins, double *data, int *lbins, double *moms, int *nbins);

void Binned_Variogram_2(double *bins, double *data, int *lbins, double *moms, int *nbins);

void Binned_Variogram_st(double *bins, double *bint, double *data, int *lbins,
			 int *lbinst, int *lbint, double *moms, double *momst,
			 double *momt, int *nbins, int *nbint);

void Cloud_Madogram(double *bins, double *data, int *lbins, double *moms, int *nbins);

void Cloud_Variogram(double *bins, double *data, int *lbins, double *moms, int *nbins);

void LeastSquare_G(double *bins, double *bint, int *cormod, double *lbins, double *moms,
		      int *nbins, int *nbint, double *nuis, double *par, double *res);

void LeastSquare_MBR(double *bins, double *bint, int *cormod, double *lbins, double *moms,
		      int *nbins, int *nbint, double *nuis, double *par, double *res);

void LeastSquare_MEG(double *bins, double *bint, int *cormod, double *lbins, double *moms,
		      int *nbins, int *nbint, double *nuis, double *par, double *res);

void LeastSquare_MET(double *bins, double *bint, int *cormod, double *lbins, double *moms,
		      int *nbins, int *nbint, double *nuis, double *par, double *res);

void WLeastSquare_G(double *bins, double *bint, int *cormod, double *lbins, double *moms,
		      int *nbins, int *nbint, double *nuis, double *par, double *res);

void WLeastSquare_MEG(double *bins, double *bint, int *cormod, double *lbins, double *moms,
		      int *nbins, int *nbint, double *nuis, double *par, double *res);

/*----------------------------------------------------------------
File name: WeightedLeastSquare.c
Description: procedures for the estimation of the model parameters
via the weighted least square method.
End
 ---------------------------------------------------------------*/

// 7)
/*----------------------------------------------------------------
File name: Utility.c
Description: procedures for the computation of useful quantities.
Start
 ---------------------------------------------------------------*/

void DeleteGlobalVar();

void RangeDist(double *max, double *min);

double Dist_geodesic(double loni, double lati, double lonj, double latj);

double Dist_chordal(double loni, double lati, double lonj, double latj);

void GeoDist(double *coordx, double *coordy, int *ncoord, double *res,int *type_dist);

void ComputeMaxima(double *df, double *maxima, int *model,
		   int *nblock, int *nsite, double *sim);

int is_equal(double val1, double val2);

double Maxima(double *x, int *size);

double Minima(double *x, int *size);

void Range(double *x, double *ran, int *size);

void Seq(double *x, int len, double *res);

void SeqStep(double *x, int len, double step, double *res);

void SetSampling(double *coordx, double *coordy, double *data, int n, int *npts,
		 double *scoordx, double *scoordy, double *sdata, double xmax,
		 double xmin, double ymax, double ymin);

void SetSampling_st(double *data,double *sdata,int *ncoord,int *ntime,
		    int wint,int k,int n,int *nrep);

void SetGlobalVar(double *coordx,double *coordy,double *coordt,int *grid,int *ia,
		  int *idx,int *ismal,int *ja,int *nsite,int *nsitex,int *nsitey,
		  int *npair,int *replic,double *srange, double *sep, int *times,double *trange,
		  int *tap,int *tapmodel,int *type,int *weighted);

void Space_Dist(double *coordx,double *coordy,int grid,int *ia,int *idx,
		int *ismal,int *ja,double thres,int type);

void SpaceTime_Dist(double *coordx,double *coordy,double *coordt,int *grid,int *ia,int *idx,int *ismal,int *ja,
                    int *tapmodel,double thres,double thret,int type);

/*----------------------------------------------------------------
File name: Utility.c
Description: procedures for the computation of useful quantities.
End
---------------------------------------------------------------*/
