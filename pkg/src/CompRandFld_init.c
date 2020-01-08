#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

// define all the globals in this file (and only)
#define extern
#include "header.h"
#undef extern

/*----------------------------------------------------------------
File name: CorrelationFunction.c
---------------------------------------------------------------*/

extern void CorrelationMat(double *rho, int *cormod, double *nuis, double *par);
extern void CorrelationMat_tap(double *rho, int *cormod, double *nuis, double *par);
extern void CorrelationMat_st(double *rho, int *cormod, double *nuis, double *par);
extern void Corr_c(double *cc,double *coordx, double *coordy, double *coordt,
                   int *cormod, int *grid, double *locx,  double *locy,int *ncoord,
                   int *nloc,int*tloc, int *ntime, double *par, int *spt, double *time,
                   int *type);
extern void Corr_c_tap(double *cc,double *cc_tap,double *coordx, double *coordy,
                       double *coordt, int *cormod, int *cormodtap, int *grid,
                       double *locx,  double *locy,int *ncoord, int *nloc,int*tloc,
                       int *ntime, double *par, int *spt, double *time,int *type);
extern void DCorrelationMat(int *cormod,double *drho,double *eps,int *flagcor,
                                int *nparcor, double *parcor,double *rho);
extern void DCorrelationMat_tap(int *cormod,double *drho,double *eps,int *flagcor,
                                    int *nparcor, double *parcor,double *rho);
extern void DCorrelationMat_st(int *cormod,double *drho,double *eps,int *flagcor,
                                     int *nparcor, double *parcor,double *rho);
extern void ExtCoeff(int *cormod, double *extc, double *lag, int *model,
          int *nlags, double *nuis, double *par);
/* extern void TapVectCorrelation(double *rho,int *cormod,double *tdists,int *ntdists,
                               double *nuis,double *par); */
extern void VectCorrelation(double *rho, int *cormod, double *h, int *nlags,
int *nlagt, double *par, double *u);


/*----------------------------------------------------------------
File name: CompositeLikelihood.c
 ---------------------------------------------------------------*/

extern void Comp_Cond_Gauss(int *cormod, double *data, double *nuis, double *par, double *thr, double *res);
extern void Comp_Cond_Gauss_st(int *cormod, double *data, double *nuis, double *par, double *thr, double *res);
extern void Comp_Cond_BinGauss(int *cormod, double *data, double *nuis, double *par, double *thr, double *res);
extern void Comp_Cond_BinGauss_st( int *cormod, double *data, double *nuis, double *par,double *thr, double *res);
extern void Comp_Diff_Gauss(int *cormod, double *data, double *nuis, double *par, double *thr, double *res);
extern void Comp_Diff_Gauss_st(int *cormod, double *data, double *nuis, double *par, double *thr, double *res);
extern void Comp_Diff_BinGauss( int *cormod, double *data, double *nuis, double *par, double *thr,double *res);
extern void Comp_Diff_BinGauss_st(int *cormod, double *data, double *nuis, double *par,double *thr, double *res);
extern void Comp_Pair_Gauss(int *cormod, double *data, double *nuis, double *par, double *thr, double *res);
extern void Comp_Pair_Gauss_st(int *cormod, double *data, double *nuis, double *par, double *thr, double *res);
extern void Comp_Pair_BinGauss(int *cormod, double *data, double *nuis, double *par, double *thr,double *res);
extern void Comp_Pair_BinGauss_st(int *cormod, double *data, double *nuis, double *par, double *thr,double *res);
extern void Comp_Brow_Resn(int *cormod, double *data, double *nuis, double *par, double *thr, double *res);
extern void Comp_Brow_Resn_st(int *cormod, double *data, double *nuis, double *par, double *thr, double *res);
extern void Comp_Ext_Gauss(int *cormod, double *data, double *nuis, double *par, double *thr, double *res);
extern void Comp_Ext_T(int *cormod, double *data, double *nuis, double *par, double *thr, double *res);

/*----------------------------------------------------------------
File name: Distributions.c
 ---------------------------------------------------------------*/

extern void Dist2Dist(double *data, double *eloc, double *escale, double *eshape, int *ndata, int *nsite, double *ploc, double *pscale, double *pshape, int *type, double *res);
extern void GevLogLik(double *data, int *ndata, double *par, double *res);
extern void vpbnorm(int *cormod, double *h, double *u, int *nlags, int *nlagt, double *nuis, double *par, double *rho, double *thr);

/*----------------------------------------------------------------
File name: Godambe.c
 ---------------------------------------------------------------*/

extern void God_Cond_Gauss(int *cormod, double *data, double *eps, int *flagcor,
                           int *flagnuis, int *npar, int *nparc, double *parcor,
                           double *nuis, double *score, double *sensmat, double *varimat);
extern void God_Cond_Gauss_st(int *cormod, double *data, double *eps, int *flagcor,
                              int *flagnuis,int *npar, int *nparc, double *parcor,
                              double *nuis, double *score, double *sensmat, double *varimat);
extern void God_Diff_Gauss(int *cormod, double *data, double *eps, int *flagcor,
                          int *flagnuis, int *npar, int *nparc, double *parcor,
                          double *nuis, double *score, double *sensmat, double *varimat);
extern void God_Diff_Gauss_st(int *cormod, double *data, double *eps, int *flagcor,
                              int *flagnuis, int *npar, int *nparc, double *parcor,
                              double *nuis, double *score,double *sensmat, double *varimat);
extern void God_Pair_Gauss(int *cormod, double *data, double *eps, int *flagcor,
                           int *flagnuis, int *npar, int *nparc, double *parcor,
                           double *nuis, double *score, double *sensmat, double *varimat);
extern void God_Pair_Gauss_st(int *cormod, double *data, double *eps, int *flagcor,
                              int *flagnuis, int *npar, int *nparc, double *parcor,
                              double *nuis, double *score, double *sensmat, double *varimat);
extern void God_BrowResn(int *cormod, double *data, double *eps, int *flagcor,
                         int *flagnuis, int *npar, int *nparc, double *parcor,
                         double *nuis, double *score, double *sensmat, double *varimat);
extern void God_BrowResn_st(int *cormod, double *data, double *eps, int *flagcor,
                            int *flagnuis, int *npar, int *nparc, double *parcor,
                            double *nuis, double *score, double *sensmat, double *varimat);
extern void God_Ext_Gauss(int *cormod, double *data, double *eps, int *flagcor,
                          int *flagnuis, int *npar, int *nparc, double *parcor,
                          double *nuis, double *score, double *sensmat, double *varimat);
extern void God_Ext_T(int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
                        int *npar, int *nparc, double *parcor, double *nuis, double *score,
                        double *sensmat, double *varimat);
extern void GodambeMat_emp(int *cormod, double *data, double *eps, int *flagcor,
                           int *flagnuis, int *like, int *model, int *npar, int *nparc,
                           double *parcor, double *nuis, double *score, double *sensmat,
                           int *spt, double *varimat, int *type);
extern void GodambeMat(double *coordx, double *coordy, int *cormod, double *data,
                       int *dist, double *eps,int *flagcor, int *flagnuis, int *grid,
                       int *like, int *model, int *npar, int *nparc, double *parcor,
                       double *nuis, double *score, double *sensmat, int *spt,
                       double *thr, int *type, double *varimat, int *vartype,
                       double *winc, double *winstp);

/*----------------------------------------------------------------
File name: Likelihood.c
 ---------------------------------------------------------------*/

extern void Binned_Lorelogram(double *bins, double *data, int *lbins, double *moms,int *nbins);
extern void Binned_Lorelogram_st(double *bins, double *bint, double *data, int *lbins,
              int *lbinst, int *lbint, double *moms,double *momst,
              double *momt, int *nbins, int *nbint);
extern void Binned_Madogram(double *bins, double *data, int *lbins, double *moments, int *nbins);
extern void Binned_Madogram_st(double *bins, double *bint, double *data, int *lbins,
            int *lbinst, int *lbint, double *moms,double *momst,
            double *momt, int *nbins, int *nbint);
extern void Binned_Variogram(double *bins, double *data, int *lbins, double *moms, int *nbins);
extern void Binned_Variogram_2(double *bins, double *data, int *lbins, double *moms, int *nbins);
extern void Binned_Variogram_st(double *bins, double *bint, double *data, int *lbins,
             int *lbinst, int *lbint, double *moms, double *momst,
             double *momt, int *nbins, int *nbint);
extern void Cloud_Madogram(double *bins, double *data, int *lbins, double *moms, int *nbins);
extern void Cloud_Variogram(double *bins, double *data, int *lbins, double *moms, int *nbins);
extern void LeastSquare_G(double *bins, double *bint, int *cormod, double *lbins, double *moms,
              int *nbins, int *nbint, double *nuis, double *par, double *res);
extern void LeastSquare_MBR(double *bins, double *bint, int *cormod, double *lbins, double *moms,
              int *nbins, int *nbint, double *nuis, double *par, double *res);
extern void LeastSquare_MEG(double *bins, double *bint, int *cormod, double *lbins, double *moms,
              int *nbins, int *nbint, double *nuis, double *par, double *res);
extern void LeastSquare_MET(double *bins, double *bint, int *cormod, double *lbins, double *moms,
              int *nbins, int *nbint, double *nuis, double *par, double *res);
extern void WLeastSquare_G(double *bins, double *bint, int *cormod, double *lbins, double *moms,
              int *nbins, int *nbint, double *nuis, double *par, double *res);
extern void WLeastSquare_MEG(double *bins, double *bint, int *cormod, double *lbins, double *moms,
              int *nbins, int *nbint, double *nuis, double *par, double *res);

/*----------------------------------------------------------------
File name: Utility.c
 ---------------------------------------------------------------*/

extern void DeleteGlobalVar();
extern void RangeDist(double *max, double *min);
extern void GeoDist(double *coordx, double *coordy, int *ncoord, double *res,int *type_dist);
extern void ComputeMaxima(double *df, double *maxima, int *model, int *nblock, int *nsite, double *sim);
extern void Seq(double *x, int len, double *res);
extern void SetGlobalVar(double *coordx,double *coordy,double *coordt,int *grid,int *ia,
                             int *idx,int *ismal,int *ja,int *nsite,int *nsitex,int *nsitey,
                             int *npair,int *replic,double *srange, double *sep, int *times,
                         double *trange,int *tap,int *tapmodel,int *type,int *weighted);
extern void Space_Dist(double *coordx,double *coordy,int grid,int *ia,int *idx,
                            int *ismal,int *ja,double thres,int type);
extern void SpaceTime_Dist(double *coordx,double *coordy,double *coordt,int *grid,
                           int *ia,int *idx,int *ismal,int *ja, int *tapmodel,
                           double thres,double thret,int type);


static const R_CMethodDef CEntries[] = {
    /* -------------------------- CorrelationFunction.c -----------*/
    {"CorrelationMat",        (DL_FUNC) &CorrelationMat,          4},
    {"CorrelationMat_tap",    (DL_FUNC) &CorrelationMat_tap,      4},
    {"CorrelationMat_st",     (DL_FUNC) &CorrelationMat_st,       4},
    {"Corr_c",                (DL_FUNC) &Corr_c,                 16},
    {"Corr_c_tap",            (DL_FUNC) &Corr_c_tap,             18},
    {"DCorrelationMat",       (DL_FUNC) &DCorrelationMat,         7},
    {"DCorrelationMat_tap",   (DL_FUNC) &DCorrelationMat_tap,     7},
    {"DCorrelationMat_st",    (DL_FUNC) &DCorrelationMat_st,      7},
    {"ExtCoeff",              (DL_FUNC) &ExtCoeff,                7},
    {"VectCorrelation",       (DL_FUNC) &VectCorrelation,         7},
    /* -------------------------- CompositeLikelihood.c -----------*/
    {"Comp_Cond_Gauss",       (DL_FUNC) &Comp_Cond_Gauss,         6},
    {"Comp_Cond_Gauss_st",    (DL_FUNC) &Comp_Cond_Gauss_st,      6},
    {"Comp_Cond_BinGauss",    (DL_FUNC) &Comp_Cond_BinGauss,      6},
    {"Comp_Cond_BinGauss_st", (DL_FUNC) &Comp_Cond_BinGauss_st,   6},
    {"Comp_Diff_Gauss",       (DL_FUNC) &Comp_Diff_Gauss,         6},
    {"Comp_Diff_Gauss_st",    (DL_FUNC) &Comp_Diff_Gauss_st,      6},
    {"Comp_Diff_BinGauss_st", (DL_FUNC) &Comp_Diff_BinGauss_st,   6},
    {"Comp_Pair_Gauss",       (DL_FUNC) &Comp_Pair_Gauss,         6},
    {"Comp_Pair_BinGauss",    (DL_FUNC) &Comp_Pair_BinGauss,      6},
    {"Comp_Pair_BinGauss_st", (DL_FUNC) &Comp_Pair_BinGauss_st,   6},
    {"Comp_Brow_Resn",        (DL_FUNC) &Comp_Brow_Resn,          6},
    {"Comp_Brow_Resn_st",     (DL_FUNC) &Comp_Brow_Resn_st,       6},
    {"Comp_Ext_Gauss",        (DL_FUNC) &Comp_Ext_Gauss,          6},
    {"Comp_Ext_T",            (DL_FUNC) &Comp_Ext_T,              6},
    /* -------------------------- Distributions.c -----------------*/
    {"Dist2Dist",             (DL_FUNC) &Dist2Dist,              11},
    {"GevLogLik",             (DL_FUNC) &GevLogLik,               4},
    {"vpbnorm",               (DL_FUNC) &vpbnorm,                 9},
    /* -------------------------- Godambe.c -----------------------*/
    {"God_Cond_Gauss",        (DL_FUNC) &God_Cond_Gauss,         12},
    {"God_Cond_Gauss_st",     (DL_FUNC) &God_Cond_Gauss_st,      12},
    {"God_Diff_Gauss",        (DL_FUNC) &God_Diff_Gauss,         12},
    {"God_Diff_Gauss_st",     (DL_FUNC) &God_Diff_Gauss_st,      12},
    {"God_Pair_Gauss",        (DL_FUNC) &God_Pair_Gauss,         12},
    {"God_Pair_Gauss_st",     (DL_FUNC) &God_Pair_Gauss_st,      12},
    {"God_BrowResn",          (DL_FUNC) &God_BrowResn,           12},
    {"God_BrowResn_st",       (DL_FUNC) &God_BrowResn_st,        12},
    {"God_Ext_Gauss",         (DL_FUNC) &God_Ext_Gauss,          12},
    {"God_Ext_T",             (DL_FUNC) &God_Ext_T,              12},
    {"GodambeMat_emp",        (DL_FUNC) &GodambeMat_emp,         16},
    {"GodambeMat",            (DL_FUNC) &GodambeMat,             24},
    /* -------------------------- Likelihood.c --------------------*/
    {"Binned_Lorelogram",     (DL_FUNC) &Binned_Lorelogram,       5},
    {"Binned_Lorelogram_st",  (DL_FUNC) &Binned_Lorelogram_st,   11},
    {"Binned_Madogram",       (DL_FUNC) &Binned_Madogram,         5},
    {"Binned_Madogram_st",    (DL_FUNC) &Binned_Madogram_st,     11},
    {"Binned_Variogram",      (DL_FUNC) &Binned_Variogram,        5},
    {"Binned_Variogram_2",    (DL_FUNC) &Binned_Variogram_2,      5},
    {"Binned_Variogram_st",   (DL_FUNC) &Binned_Variogram_st,    11},
    {"Cloud_Madogram",        (DL_FUNC) &Cloud_Madogram,          5},
    {"Cloud_Variogram",       (DL_FUNC) &Cloud_Variogram,         5},
    {"LeastSquare_G",         (DL_FUNC) &LeastSquare_G,          10},
    {"LeastSquare_MBR",       (DL_FUNC) &LeastSquare_MBR,        10},
    {"LeastSquare_MEG",       (DL_FUNC) &LeastSquare_MEG,        10},
    {"LeastSquare_MET",       (DL_FUNC) &LeastSquare_MET,        10},
    {"WLeastSquare_G",        (DL_FUNC) &WLeastSquare_G,         10},
    {"WLeastSquare_MEG",      (DL_FUNC) &WLeastSquare_MEG,       10},
    /* -------------------------- Utility.c -----------------------*/
    {"DeleteGlobalVar",       (DL_FUNC) &DeleteGlobalVar,         0},
    {"RangeDist",             (DL_FUNC) &RangeDist,               2},
    {"GeoDist",               (DL_FUNC) &GeoDist,                 5},
    {"ComputeMaxima",         (DL_FUNC) &ComputeMaxima,           6},
    {"Seq",                   (DL_FUNC) &Seq,                     3},
    {"SetGlobalVar",          (DL_FUNC) &SetGlobalVar,           21},
    {"Space_Dist",            (DL_FUNC) &Space_Dist,              9},
    {"SpaceTime_Dist",        (DL_FUNC) &SpaceTime_Dist,         12},
    {NULL, NULL, 0}
};

void R_init_CompoRandFld(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

void R_unload_CompRndFld(DllInfo *info) {
  // just to avoid warning from compiler
}
