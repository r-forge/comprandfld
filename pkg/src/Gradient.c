#include "header.h"

/*###################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Email: simone.padoan@epfl.ch.
### Institute: EPFL.
### File name: Gradient.c
### Description:
### This file contains a set of procedures
### for the computation of the composite likelihood
### gradients.
### Last change: 05/12/2010.
##################################################*/

// Compute the gradient vector of the conditional pairwise log likelihood for a Gaussian model:
void Grad_Cond_Gauss(double rho, int *flag, double *gradcor, double *grad,
		     int *npar, double *par, double u, double v)
{
  // Initialization variables:
  double mean=par[0], nugget=par[1], sill=par[2];
  double C=0.0, L=0.0, M=0.0, Q=0.0, R=0.0, S=0.0;
  double a=0.0, b=0.0, c=0.0, d=0.0, da=0.0, e=0.0;
  double pu=0.0, pv=0.0, su=0.0, sv=0.0;
  double suv=0.0, uv=0.0;
  int h=0, i=0, j=0;
  //defines useful quantities:
  a=nugget+sill;
  b=sill*rho;
  c=pow(a,2)-pow(b,2);
  d=a/c;
  da=2*a;
  e=b/c;
  u=u-mean;
  v=v-mean;
  uv=u+v;
  pu=pow(u,2);
  pv=pow(v,2);
  R=pu+pv;
  S=0.5*R/c;
  L=u*v;
  M=L/c;
  Q=d*R-2*e*L-1;
  su=(-1+pu/a)/da; //first statistics: first component
  sv=(-1+pv/a)/da; //second statistics: second component
  suv=su+sv;
  // Derivatives of the conditional respect with the mean
  if(flag[0]==1){ grad[i]=uv*(2/(a+b)-1/a); i++; }
  // Derivative of the conditional respect with the nugget
  if(flag[1]==1){ grad[i]=2*(d*Q-S)-suv; i++; }
  // Derivative of the conditional respect with the sill
  if(flag[2]==1){ grad[i]=2*((d-e)*Q-S-rho*M)-suv; i++; }
  // Derivatives with respect to the correlation parameters
  h=0;
  C=sill*(M-e*Q);
  for(j=i;j<*npar;j++) { grad[j]=2*C*gradcor[h]; h++; }
  return;
}

// Compute the gradient vector of the conditional log likelihood for a Gaussian-Binary model :
void Grad_Cond_Bin(double rho,double pij, double p,int *flag, double *gradcor, double *grad,
		   int *npar, double *nuis, double *thr, double u, double v)
{
  // Initialization variables:
  double dpij=0.0, dij=0.0, rvar=0.0, dpdm=0.0, f=0.0;
  double q1=0.0, q2=0.0, q3=0.0, sh=0.0, vario=0.0, z=0;
  int h=0, i=0, j=0;
  //init variables:
  z=(nuis[0]-*thr)/sqrt(nuis[2]+nuis[1]);
  rvar=nuis[2]/(nuis[2]+nuis[1]);
  //set derivatives components:
  q1=dnorm(z,0,1,0);//stand normal pdf
  q2=pnorm(z*sqrt((1-rvar*rho)/(1+rvar*rho)),0,1,1,0);// stand norm cdf
  q3=d2norm(z,z,rvar*rho);// biv stand norm pdf
  //derivatives:
  dpdm=q1/sqrt(nuis[2]+nuis[1]);/*dp/dmu*/
  dpij=2*dpdm*q2;/*dpij/dmu*/
  f=-(0.5*(nuis[0]-*thr)*dpdm)/(nuis[2]+nuis[1]);/* dp/dsill*/
  dij=2*f*q2;/* dpij/dsill*/
  vario=2*(p-pij);//variogramma binario!!!
  sh=1/(1-2*p+pij);
  // Derivative of the difference respect with the mean
  if(flag[0]==1) { grad[i]=(dpij-2*dpdm)*(1-((u+v)*nij(dpij,dpdm,pij,p)+
					   (u*v)*mij(dpij,dpdm,pij,p)))*sh+dpdm*(1-(u+v)/(2*p))/(1-p); i++; }
  // Derivative of the difference respect with the nugget
  if(flag[1]==1) { grad[i]=1; i++; }
  // Derivative of the difference respect with the sill
  if(flag[2]==1) { grad[i]=(dij-2*f)*(1-((u+v)*nij(dij,f,pij,p)+
					 (u*v)*mij(dij,f,pij,p)))*sh+f*(1-(u+v)/(2*p))/(1-p); i++; }
  // Derivatives with respect to the correlation parameters
  for(j=i;j<*npar;j++) { grad[j]=gradcor[j]*q3*rvar*(1-((u+v)*2*(p-1)/vario +
                                                (u*v)*2*(pij-2*pow(p,2)+p)/(vario*pij)))*sh; h++; }

  return;
}

// Compute the gradient vector of the difference log likelihood for a Gaussian model :
void Grad_Diff_Gauss(double rho, int *flag, double *gradcor, double *grad,
		     int *npar, double *par, double u, double v)
{
  // Initialization variables:
  double nugget=par[1], sill=par[2];
  double vario=0.0, sh=0.0;
  int h=0, i=0, j=0;
  //variogram:
  vario=nugget+sill*(1-rho);
  sh=0.5*(0.5*pow(u-v,2)/vario-1)/vario;
  // Derivative of the conditional respect with the nugget
  if(flag[1]==1) { grad[i]=sh; i++; }
  // Derivative of the conditional respect with the sill
  if(flag[2]==1) { grad[i]=(1-rho)*sh; i++; }
  // Derivatives with respect to the correlation parameters
  for(j=i;j<*npar;j++) { grad[j]=-sill*gradcor[h]*sh; h++; }
  return;
}

// Compute the gradient vector of the difference log likelihood for a Gaussian-Binary model :
void Grad_Diff_Bin(double rho,double pij, double p,int *flag, double *gradcor,  double *grad,
		   int *npar, double *nuis, double *thr, double u, double v)
{
  // Initialization variables:
  double dpij=0.0, dij=0.0, rvar=0.0, dpdm=0.0, f=0.0;
  double q1=0.0, q2=0.0, q3=0.0, sh=0.0, vario=0.0, z=0;
  int h=0, i=0, j=0;
  //init variables:
  z=(nuis[0]-*thr)/sqrt(nuis[2]+nuis[1]);
  rvar=nuis[2]/(nuis[2]+nuis[1]);
  //set derivatives components:
  q1=dnorm(z,0,1,0);//stand normal pdf
  q2=pnorm(z*sqrt((1-rvar*rho)/(1+rvar*rho)),0,1,1,0);// stand norm cdf
  q3=d2norm(z,z,rvar*rho);// biv stand norm pdf
  //derivatives:
  dpdm=q1/sqrt(nuis[2]+nuis[1]);/*dp/dmu*/
  dpij=2*dpdm*q2;/*dpij/dmu*/
  f=-(0.5*(nuis[0]-*thr)*dpdm)/(nuis[2]+nuis[1]);/* dp/dsill*/
  dij=2*f*q2;/* dpij/dsill*/
  vario=2*(p-pij);//variogramma binario!!!
  sh=2*(1-pow(u-v,2)/vario)/(1-vario);
  // Derivative of the difference respect with the mean
  if(flag[0]==1) { grad[i]=(dpij-dpdm)*sh; i++; }
  // Derivative of the difference respect with the nugget
  if(flag[1]==1) { grad[i]=1; i++; }
  // Derivative of the difference respect with the sill
  if(flag[2]==1) { grad[i]=(dij-f)*sh; i++; }
  // Derivatives with respect to the correlation parameters
  for(j=i;j<*npar;j++) { grad[j]=gradcor[j]*q3*rvar*sh; h++; }
  return;
}


// Compute the gradient vector of the variogram for a Gaussian model :
void Grad_Diff_Vario(double rho, int *flag, double *gradcor,
		     double *grad, int *npar, double *par)
{
  // Initialization variables:
  double nugget=par[1], sill=par[2];
  double vario=0.0, sh=0.0;
  int h=0, i=0, j=0;
  //variogram:
  vario=nugget+sill*(1-rho);
  sh=1/vario;
  // Derivative of the conditional respect with the nugget
  if(flag[1]==1){ grad[i]=sh; i++; }
  // Derivative of the conditional respect with the sill
  if(flag[2]==1){ grad[i]=(1-rho)*sh; i++; }
  // Derivatives with respect to the correlation parameters
  for(j=i;j<*npar;j++){ grad[j]=-sill*gradcor[h]*sh; h++;}
  return;
}

// Compute the gradient vector of the pairwise log likelihood for a Gaussian-Binary model :
void Grad_Pair_Bin(double rho,double pij, double p,int *flag, double *gradcor, double *grad,
		   int *npar, double *nuis, double *thr, double u, double v)
{
  // Initialization variables:
  double dpij=0.0, dij=0.0, rvar=0.0, dpdm=0.0, f=0.0;
  double q1=0.0, q2=0.0, q3=0.0, sh=0.0, vario=0.0, z=0;
  int h=0, i=0, j=0;
  //init variables:
  z=(nuis[0]-*thr)/sqrt(nuis[2]+nuis[1]);
  rvar=nuis[2]/(nuis[2]+nuis[1]);
  //set derivatives components:
  q1=dnorm(z,0,1,0);//stand normal pdf
  q2=pnorm(z*sqrt((1-rvar*rho)/(1+rvar*rho)),0,1,1,0);// stand norm cdf
  q3=d2norm(z,z,rvar*rho);// biv stand norm pdf
  //derivatives:
  dpdm=q1/sqrt(nuis[2]+nuis[1]);/*dp/dmu*/
  dpij=2*dpdm*q2;/*dpij/dmu*/
  f=-(0.5*(nuis[0]-*thr)*dpdm)/(nuis[2]+nuis[1]);/* dp/dsill*/
  dij=2*f*q2;/* dpij/dsill*/
  vario=2*(p-pij);//variogramma binario!!!
  sh=1/(1-2*p+pij);
  // Derivative of the difference respect with the mean
  if(flag[0]==1) { grad[i]=(dpij-2*dpdm)*(1-((u+v)*nij(dpij,dpdm,pij,p)+
					 (u*v)*mij(dpij,dpdm,pij,p)))*sh; i++; }
  // Derivative of the difference respect with the nugget
  if(flag[1]==1) { grad[i]=1; i++; }
  // Derivative of the difference respect with the sill
  if(flag[2]==1) { grad[i]=(dij-2*f)*(1-((u+v)*nij(dij,f,pij,p)+
					 (u*v)*mij(dij,f,pij,p)))*sh; i++; }
  // Derivatives with respect to the correlation parameters
  for(j=i;j<*npar;j++) { grad[j]=gradcor[j]*q3*rvar*(1-((u+v)*2*(p-1)/vario +
                                                (u*v)*2*(pij-2*pow(p,2)+p)/(vario*pij)))*sh; h++; }
  return;
}


// Compute the gradient vector of the pairwise log likelihood for a Gaussian model:
void Grad_Pair_Gauss(double rho, int *flag, double *gradcor, double *grad, int *npar,
		      double *par, double u, double v)
{
  // Initialization variables:
  double mean=par[0], nugget=par[1], sill=par[2];
  double C=0.0, L=0.0, M=0.0, Q=0.0, R=0.0, S=0.0;
  double a=0.0, b=0.0, c=0.0, d=0.0, e=0.0;
  int h=0, i=0, j=0;
  // defines useful quantities:
  a=nugget+sill;
  b=sill*rho;
  c=pow(a,2)-pow(b,2);
  d=a/c;
  e=b/c;
  u=u-mean;
  v=v-mean;
  R=pow(u,2)+pow(v,2);
  S=0.5*R/c;
  L=u*v;
  M=L/c;
  Q=d*R-2*e*L-1;
  // Derivative of the pairwise respect with the mean
  if(flag[0]==1) { grad[i]=(u+v)/(a+b); i++; }
  // Derivative of the pairwise respect with the nugget
  if(flag[1]==1) { grad[i]=d*Q-S; i++; }
  // Derivative of the pairwise respect with the sill
  if(flag[2]==1) { grad[i]=(d-e)*Q-S-rho*M; i++; }
  // Derivatives with respect to the correlation parameters
  C=sill*(M-e*Q);
  for(j=i;j<*npar;j++) { grad[j]=C*gradcor[h]; h++;}
  return;
}
// Gradient of the max-stable Extremal Gaussian model:
void Grad_Ext_Gauss(double rho, int *flag, double *gradcor, double *grad,
		    int *npar, double *par, double u, double v)
{
  // Initialization variables:
  double a=0.0, a3=0.0, C=0.0, omr2=0.0, d2V=0.0, drV=0.0;
  double drd2V=0.0, drduV=0.0, drdvV=0.0, duV=0.0, dvV=0.0;
  double u2=0.0, v2=0.0;
  int h=0, i=0, j=0;
  // defines useful quantities:
  u2=pow(u,2);
  v2=pow(v,2);
  rho=par[0]*rho;//rho=sill*corr
  omr2=1-pow(rho,2);
  //sqrt of the quadratic form
  a=sqrt(u2+v2-2*rho*u*v);
  a3=pow(a,3);
  duV=-0.5*(u*rho-a-v)/(u2*a);
  dvV=-0.5*(v*rho-a-u)/(v2*a);
  d2V=0.5*omr2/a3;
  drV=0.5/a;
  drduV=-0.5*(u-v*rho)/a3;
  drdvV=-0.5*(v-u*rho)/a3;
  drd2V=-rho/a3+(3*u*v*omr2)/pow(a,5);
  C=drV+(drd2V+drduV*dvV+duV*drdvV)/(d2V+duV*dvV);
  // Derivative of the pairwise respect with the sill
  if(flag[0]==1) { grad[i]=C*rho/par[0]; i++; }
  // Derivatives with respect to the correlation parameters
  for(j=i;j<*npar;j++) { grad[j]=C*par[0]*gradcor[h]; h++;}
  return;
}
// Gradient of the max-stable extremal-t model:
void Grad_Ext_T(double rho, int *flag, double *gradcor, double *grad,
		int *npar, double *par, double x, double y)
{
  // Initialization variables:
  double a=0.0, ac=0.0, acb=0.0, acd=0.0, aci=0.0, acib=0.0, acie=0.0;
  double A=0.0, b=0.0, B=0.0, c=0.0, clc=0.0, cilci=0.0, C=0.0, d=0.0;
  double d1ty=0.0, d2tx=0.0, d2ty=0.0, d2V=0.0, d1tx=0.0, df12=0.0;
  double dd1tx=0.0, dd1ty=0.0, dd2V=0.0, ddfa=0.0, ddfac=0.0, ddfaci=0.0;
  double ddfc=0.0, ddfci=0.0, ddfd2V=0.0, ddfdxV=0.0, ddfdyV=0.0, ddfV=0.0;
  double ddfw=0.0, ddxV=0.0, ddyV=0.0, ddfz=0.0, ddtx=0.0, ddty=0.0, den=0.0;
  double df=par[0], df1=par[0]+1, dptx=0.0, dpty=0.0, dtx=0.0, dtxd=0.0;
  double dty=0.0, dtye=0.0, dxV=0.0, dyV=0.0, dV=0.0, D=0.0, e=0.0, E=0.0;
  double F=0.0, G=0.0, omr2=0.0, opdf=1+1/par[0], ptx=0.0, pty=0.0, somr2=0.0;
  double x2=0.0, x2d=0.0, xy=0.0, xy2d=0.0, xyd=0.0, y2=0.0, y2d=0.0;
  double yx2d=0.0, w=0.0, z=0.0;
  int i=0, j=0, h=0;
  // defines useful quantities:
  rho=par[1]*rho;//rho=sill*corr
  omr2=1-pow(rho, 2);
  somr2 = sqrt(omr2);
  a=sqrt(df1/omr2);
  c=pow(y/x,1/df);
  ac=a*c;
  aci=a/c;
  z=(c-rho)*a;
  w=(1/c-rho)*a;
  x2=pow(x,2);
  y2=pow(y,2);
  x2d=x2*df;
  y2d=y2*df;
  yx2d=x2d*y;
  xy2d=y2d*x;
  xy=x*y;
  xyd=xy*df;
  ptx=pt(z,df1,1,0);
  pty=pt(w,df1,1,0);
  dtx=dt(z,df1,0);
  dty=dt(w,df1,0);
  d1tx=d1x_dt(z,df1);
  d1ty=d1x_dt(w,df1);
  //defines the log-likelihood components:
  A=dtx*opdf+d1tx*ac/df;
  B=dty*opdf+d1ty*aci/df;
  dxV=ptx/x2+dtx*ac/x2d-dty*aci/xyd;
  dyV=pty/y2+dty*aci/y2d-dtx*ac/xyd;
  d2V=ac*A/x2d/y+aci*B/y2d/x;
  den=d2V+dxV*dyV;
  // Derives with respect to the degree of freedom:
  if(flag[0]==1)
    {
      clc=c*log(c);
      cilci=log(1/c)/c;
      df12=2*df1;
      dptx=ddf_pt(z,df1);
      dpty=ddf_pt(w,df1);
      ddtx=ddf_t_dt(z,clc,df1,somr2);
      ddty=ddf_t_dt(w,cilci,df1,somr2);
      dd1tx=ddf_t_d1x_dt(z,clc,df1,somr2);
      dd1ty=ddf_t_d1x_dt(w,cilci,df1,somr2);
      ddfc=-clc/df;
      ddfci=log(c)/c/df;
      ddfa=1/(2*a*omr2);
      ddfac=ddfa*c+a*ddfc;
      ddfaci=ddfa/c+a*ddfci;
      ddfz=z/df12+ddfc*a;
      ddfw=w/df12+ddfci*a;
      E=dptx+dtx*ddfz;
      F=dpty+dty*ddfw;
      // Defines the derivatives of the log-likelihood:
      ddfV=-E/x-F/y;
      ddfdxV=E/x2+(ddtx*ac+dtx*(ddfac-ac/df))/x2d-
      	(ddty*aci+dty*(ddfaci-aci/df))/xyd;
      ddfdyV=F/y2+(ddty*aci+dty*(ddfaci-aci/df))/y2d-
	(ddtx*ac+dtx*(ddfac-ac/df))/xyd;
      ddfd2V=((ddfac-ac/df)*A+
	      ac*(dd1tx*ac+ddtx*df1+d1tx*(ddfac-ac/df)-dtx/df)/df)/yx2d+
	((ddfaci-aci/df)*B+
	 aci*(dd1ty*aci+ddty*df1+d1ty*(ddfaci-aci/df)-dty/df)/df)/xy2d;
      grad[i]=ddfV+(ddfd2V+ddfdxV*dyV+dxV*ddfdyV)/den;
      i++;
    }
  // Derivatives with the respect to the sill and correlation parameters:
  b=rho/omr2;
  d=-a+b*z;
  e=-a+b*w;
  acd=ac*d;
  acb=ac*b;
  acib=aci*b;
  acie=aci*e;
  d2tx=d2x_dt(z,df1);
  d2ty=d2x_dt(w,df1);
  dtxd=dtx*d;
  dtye=dty*e;
  C=d1tx*acd+dtx*acb;
  D=d1ty*acie+dty*acib;
  // Defines the derivatives of the log-likelihood:
  dV=-dtxd/x-dtye/y;
  ddxV=dtxd/x2+C/x2d-D/xyd;
  ddyV=dtye/y2+D/y2d-C/xyd;
  dd2V=acb*A/yx2d+ac*(d1tx*d*opdf+d2tx*acd/df+d1tx*acb/df)/yx2d+
    acib*B/xy2d+aci*(d1ty*e*opdf+d2ty*acie/df+d1ty*acib/df)/xy2d;
  G=dV+(dd2V+ddxV*dyV+dxV*ddyV)/den;
  // Derivatives with respect to the sill parameter
  if(flag[1]==1){ grad[i]=G*rho/par[1]; i++;}
  // Derivatives with respect to the correlation parameters
  for(j=i;j<*npar;j++){ grad[j]=G*par[1]*gradcor[h]; h++;}
  return;
}
// Gradient of the max-stable Brown-Resnick model:
void Grad_Brow_Resn(double vario, int *flag, double *gradcor, double *grad,
		    int *npar, double *par, double x, double y)
{
  // Initialization variables:
  double a=0.0, a2x=0.0, a2y=0.0, ao2=0.0, ax=0.0, ay=0.0;
  double axy=0.0, ax2=0.0, ay2=0.0, dx=0.0, dy=0.0, C=0.0;
  double d2V=0.0, daV=0.0, dad2V=0.0, dadxV=0.0, dadyV=0.0;
  double dxV=0.0, dyV=0.0, lyx=0.0, omz2=0.0, omw2=0.0;
  double opzw=0.0, px=0.0, py=0.0, x2=0.0, y2=0.0, w=0.0, z=0.0;
  int i=0;
  // defines useful quantities:
  a=sqrt(vario);// Husler-Reiss coefficient (lambda)
  ao2=0.5*a;
  ax=a*x;
  ay=a*y;
  a2x=a*ax;
  a2y=a*ay;
  axy=ax*y;
  x2=pow(x,2);
  y2=pow(y,2);
  ax2=a*x2;
  ay2=a*y2;
  lyx=log(y/x)/a;
  z=ao2+lyx;
  w=ao2-lyx;
  opzw=1+z*w;
  omz2=1-pow(z,2);
  omw2=1-pow(w,2);
  px=pnorm(z,0,1,1,0);
  py=pnorm(w,0,1,1,0);
  dx=dnorm(z,0,1,0);
  dy=dnorm(w,0,1,0);
  dxV=px/x2+dx/ax2-dy/axy;
  dyV=py/y2+dy/ay2-dx/axy;
  d2V=(w*dx*y+z*dy*x)/(ax2*ay2);
  daV=-w*dx/ax-z*dy/ay;
  dadxV=(opzw*dy/y-omw2*dx/x)/a2x;
  dadyV=(opzw*dx/x-omz2*dy/y)/a2y;
  dad2V=(z*omw2-2*w)*dx/(ax2*a2y)+
    (w*omz2-2*z)*dy/(ay2*a2x);
  dad2V=(z-z*w*w-2*w)*dx/ax2/a2y+(w-z*z*w-2*z)*dy/ay2/a2x;
  C=daV+(dad2V+dadxV*dyV+dxV*dadyV)/(d2V+dxV*dyV);
  // Derivatives with respect to the variogram parameters
  for(i=0;i<*npar;i++) grad[i]=0.5*C*gradcor[i]/a;
  return;
}


double mij(double qij, double w, double pij, double p)
{
  double val=0.0;
  val=(2*pij*w*(p-1)+qij*(p+pij-2*pow(p,2)))/(pij*(p-pij)*(qij-2*w));
  return(val);
}


double nij(double qij, double w, double pij, double p)
{
  double val=0.0;
  val=(w*(1-pij)-qij*(1-p))/((p-pij)*(qij-2*w));
  return(val);
}
