#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>

#define LOW -1.0e15
#define MAXERR 1e-6

//---------START GLOBAL VARIABLES-----------
extern int *istap;//is tapering?
extern int *isst;//is a spatio-temporal random field?
extern double *lags;// vector of spatial distances for tapering
extern double *lagt;// vector of temporal distance for tapering
extern double **mlags;// vector of spatial distances
extern double **mlagt;// vector of temporal distances
extern double *maxdist;// the threshould of the spatial distances
extern double *maxtime;// the threshould of the temporal distances
// below which the pairs are considered
extern double *maximdista;// the maximum spatial distance
extern double *maximtime;// the maximum temporal distance
extern double *minimdista; // the minimum spatial distance
extern double *minimtime;// the minimum temporal distance
extern int *ncoord;// number of total spatial coordinates
extern int *ncoordx;// number of the first spatial coordinates
extern int *ncoordy;// number of the second spatial coordinates
extern int *npairs;// effective number of pairs
extern int *nrep;// number of iid replicates of the random field
extern int *ntime;// number of times
extern double *tapsep; // parameter separability for space time quasi taper
//int *totpairs;// total number of pairs
//---------END GLOBAL VARIABLES-------------
//void indx(double *ind, int *n);

// fortran declaration for bivariate normal cdf:
extern double F77_NAME(bvnmvn)(double *lower, double *upper, int* infin, double *correl);

// Internal function declaration:
// 1)
/*----------------------------------------------------------------
File name: CorrelationFunction.c
Description: procedures for computation of correlation functions
Start
---------------------------------------------------------------*/

extern double CheckCor(int *cormod, double *par);
extern void Comp_supp(double *c_supp,int *cormod, double h,double u, double *par);
extern double CorFct(int *cormod, double h, double u, double *par);
extern double CorFunCauchy(double lag, double power2, double scale);
extern double CorFunGenCauchy(double lag, double power1, double power2, double scale);
extern double CorFunSferical(double lag, double scale);
extern double CorFunStable(double lag, double power, double scale);

		       //double CorFunDobStable(double lag, double power_s, double power_t,
		       //double scale_s, double scale_t, double tsep);

extern double CorFunWitMat(double lag, double scale, double smooth);
extern double CorFunWend1(double lag,double scale);
extern double CorFunWend2(double lag,double scale);
extern double CorFunWend3(double lag,double scale);
extern double DCauchyPow(double power2, double scale, double rho);
extern double DCauchySc(double lag, double power2, double scale, double rho);
extern double DExpoSc(double lag, double scale, double rho);
extern double DGaussSc(double lag, double scale, double rho);
extern double DGenCauP1(double lag, double power1, double power2, double scale, double rho);
extern double DGenCauP2(double lag, double power1, double scale, double rho);
extern double DGenCauSc(double lag, double power1, double power2, double scale, double rho);
extern double DSferiSc(double lag, double scale);
extern double DStabPow(double lag, double power, double scale, double rho);
extern double DStabSc(double lag, double power, double scale, double rho);
extern double DWhMatSc(double *eps, double lag, double scale, double smooth);
extern double DWhMatSm(double *eps, double lag, double scale, double smooth);
extern double DGneiting_sc_s(double h,double u, double power_s,double power_t,
                             double scale_s,double scale_t,double sep);
extern double DGneiting_sc_t(double h,double u, double power_s,double power_t,
                             double scale_s,double scale_t,double sep);
extern double DGneiting_pw_s(double h,double u, double power_s,double power_t,
                             double scale_s,double scale_t,double sep);
extern double DGneiting_pw_t(double h,double u, double power_s,double power_t,
                             double scale_s,double scale_t,double sep);
extern double DGneiting_sep(double h,double u, double power_s,double power_t,
                            double scale_s,double scale_t,double sep);
extern double DIaco_sc_s(double h,double u, double power_s,double power_t,
                         double scale_s,double scale_t,double power2);
extern double DIaco_sc_t(double h,double u, double power_s,double power_t,
                         double scale_s,double scale_t,double power2);
extern double DIaco_pw_s(double h,double u, double power_s,double power_t,
                         double scale_s,double scale_t,double power2);
extern double DIaco_pw_t(double h,double u, double power_s,double power_t,
                         double scale_s,double scale_t,double power2);
extern double DIaco_pw2(double h,double u, double power_s,double power_t,
                        double scale_s,double scale_t,double power2);
extern double DPorcu_sc_s(double h,double u, double power_s,double power_t,
                          double scale_s,double scale_t,double sep);
extern double DPorcu_sc_t(double h,double u, double power_s,double power_t,
                          double scale_s,double scale_t,double sep);
extern double DPorcu_pw_s(double h,double u, double power_s,double power_t,
                          double scale_s,double scale_t,double sep);
extern double DPorcu_pw_t(double h,double u, double power_s,double power_t,
                          double scale_s,double scale_t,double sep);
extern double DPorcu_sep(double h,double u, double power_s,double power_t,
                         double scale_s,double scale_t,double sep);
extern double DStein_sc_s(double h,double u, double power_t,double scale_s,
                          double scale_t,double smooth);
extern double DStein_sc_t(double h,double u, double power_t,double scale_s,
                          double scale_t,double smooth);
extern double DStein_pw_s(double h,double u, double power_t,double scale_s,
                          double scale_t,double smooth);
extern double DStein_pw_t(double h,double u, double power_t,double scale_s,
                          double scale_t,double smooth);
extern double DStein_sm(double h,double u,   double power_t,double scale_s,
                        double scale_t,double smooth);
extern double DExp_Cauchy_sc_s(double h,double u,double scale_s,double scale_t,double power2);
extern double DExp_Cauchy_sc_t(double h,double u,double scale_s,double scale_t,double power2);
extern double DExp_Cauchy_pw2 (double h,double u,double scale_s,double scale_t,double power2);
extern double DExp_Exp_sc_s(double h,double u,double scale_s,double scale_t);
extern double DExp_Exp_sc_t(double h,double u,double scale_s,double scale_t);
extern double DMat_Exp_sc_s(double h,double u,double scale_s,double scale_t,double smooth);
extern double DMat_Exp_sc_t(double h,double u,double scale_s,double scale_t,double smooth);
extern double DMat_Exp_sm(double h,double u,double *eps,double scale_s,double scale_t,double smooth);
extern double DMat_Cauchy_sc_s(double h,double u,double power2,double scale_s,
                               double scale_t,double smooth);
extern double DMat_Cauchy_sc_t(double h, double u,double power2,double scale_s,
                               double scale_t,double smooth);
extern double DMat_Cauchy_pw2(double h,double u,double power2,double scale_s,
                              double scale_t,double smooth);
extern double DMat_Cauchy_sm(double h,double u,double *eps, double power2,
                             double scale_s,double scale_t,double smooth);
extern void GradCorrFct(double rho, int *cormod, double *eps, int *flag, double *grad, double h, double u, double *par);
extern void GradVarioFct(double vario, int *cormod, double *eps, int *flag, double *grad, double lag, double *par, double tsep);
extern double Variogram(int *cormod, double h, double u, double *nuis, double *par);
extern double VarioFct(int *cormod, double h, double *par, double u);
extern double VarioDobStable(double lag, double power_s, double power_t, double scale_s,
                             double scale_t, double tsep);
extern double VarioGCauchy(double lag, double power1, double power2, double scale);

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

extern double BrowResnllik(double a, double c, double x, double y);


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

extern double d1x_dt(double x, double df);
extern double d2x_dt(double x, double df);
extern double ddf_pt(double x, double df);
extern double ddf_t_dt(double x, double clc, double df, double somc2);
extern double ddf_t_d1x_dt(double x, double clc, double df, double somc2);
extern double dts(double x, double df);
extern double dgev(double x, double loc, double scale, double shape);
extern double d2norm(double x, double y, double rho);
extern double int_pt(double x, double df);
extern void integr_pt(double *x, int n, void *ex);
extern double pbnorm(int *cormod, double h, double u, double *nuis, double *par, double thr);
extern double pgev(double x, double loc, double scale, double shape);
extern double qgev(double x, double loc, double scale, double shape);


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

extern void GodambeMat_Diff(double *coordx, double *coordy, int *cormod, double *eps, int *flagcor, int *flagnuis, int *model, int *npar, int *nparc, double *parcor, double *nuis, double *sensmat, double *varimat);
extern void Sens_Cond_Gauss(int *cormod, double *data, double *eps, int *flagcor,
                            int *flagnuis, double *nuis, int *np, int *npar, int *nparc,
                            double *parcor, double *score, double *sensmat);
extern void Sens_Cond_Gauss_st(int *cormod, double *data, double *eps, int *flagcor,
                               int *flagnuis, double *nuis, int *np, int *npar, int *nparc,
                               double *parcor, double *score, double *sensmat);
extern void Sens_Cond_Gauss_ij(double rho, int *flag, double *gradcor, int *npar,
                                     int *nparc, double *par, double *sensmat);
extern void Sens_Diff_Gauss(int *cormod, double *data, double *eps, int *flagcor,
                            int *flagnuis, double *nuis, int *np,int *npar, int *nparc,
                            double *parcor,double *score, double *sensmat);
extern void Sens_Diff_Gauss_st(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
            double *nuis, int *np,int *npar, int *nparc, double *parcor, double *score, double *sensmat);
extern void Sens_Diff_Gauss_ij(double *gradient, int *npar, double *sensmat);
extern void Sens_Pair_Gauss(int *cormod, double *data, double *eps, int *flagcor,
                            int *flagnuis, double *nuis, int *np, int *npar, int *nparc,
                            double *parcor, double *score, double *sensmat);
extern void Sens_Pair_Gauss_st(int *cormod, double *data, double *eps, int *flagcor,
                               int *flagnuis, double *nuis, int *np, int *npar, int *nparc,
                               double *parcor, double *score, double *sensmat);
extern void Sens_Pair_Gauss_ij(double rho, int *flag, double *gradcor, int *npar,
                                     int *nparc, double *par, double *sensmat);
extern void Sensitivity(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
                            int *like, int *model, int *npar, int *nparc, double *parcor,
                            double *nuis, int *np, double *score, double *sensmat, int *spt, int *type);
extern void Vari_SubSamp(double *coordx, double *coordy, int *cormod, double *data,
                             int *dist, double *eps, int *flagcor, int *flagnuis, int *grid,
                         int *like,int *model, int *npar, int *nparc, double *nuis, int *np,
                         double *parcor, double *thr, int *type, double *varimat,
                         double *winc, double *winstp);
extern void Vari_SubSamp_st(int *cormod, double *data, int *dist, double *eps,
                            int *flagcor,int *flagnuis, int *like, int *npar,
                            int *nparc, double *nuis, int *np,double *parcor,
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

extern void Grad_Cond_Bin(double rho,double pij, double p,int *flag, double *gradcor,
           double *grad, int *npar, double *nuis, double *thr, double u, double v);
extern void Grad_Cond_Gauss(double rho, int *flag, double *gradcor, double *grad,
             int *npar, double *par, double u, double v);
extern void Grad_Diff_Bin(double rho,double pij, double p,int *flag, double *gradcor,
           double *grad, int *npar, double *nuis, double *thr, double u, double v);
extern void Grad_Diff_Gauss(double rho, int *flag, double *gradcor, double *grad,
             int *npar, double *par, double u, double v);
extern void Grad_Diff_Vario(double rho, int *flag, double *gradcor, double *grad,
             int *npar, double *par);
extern void Grad_Pair_Bin(double rho,double pij, double p,int *flag, double *gradcor,
           double *grad, int *npar, double *nuis, double *thr, double u, double v);
extern void Grad_Pair_Gauss(double rho, int *flag, double *gradcor, double *grad,
             int *npar, double *par, double u, double v);
extern void Grad_Ext_Gauss(double rho, int *flag, double *gradcor, double *grad,
            int *npar, double *par, double u, double v);
extern void Grad_Brow_Resn(double vario, int *flag, double *gradcor, double *grad,
            int *npar, double *par, double x, double y);
extern void Grad_Ext_T(double rho, int *flag, double *gradcor, double *grad,
        int *npar, double *par, double x, double y);
extern double mij(double qij, double w, double pij, double p);
extern double nij(double qij, double w, double pij, double p);

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

extern double Dist_geodesic(double loni, double lati, double lonj, double latj);
extern double Dist_chordal(double loni, double lati, double lonj, double latj);
extern int is_equal(double val1, double val2);
extern double Maxima(double *x, int *size);
extern double Minima(double *x, int *size);
extern void Range(double *x, double *ran, int *size);
extern void SeqStep(double *x, int len, double step, double *res);
extern void SetSampling(double *coordx, double *coordy, double *data, int n, int *npts, double *scoordx, double *scoordy, double *sdata, double xmax, double xmin, double ymax, double ymin);
extern void SetSampling_st(double *data,double *sdata,int *ncoord,int *ntime, int wint,int k,int n,int *nrep);

/*----------------------------------------------------------------
File name: Utility.c
Description: procedures for the computation of useful quantities.
End
---------------------------------------------------------------*/
