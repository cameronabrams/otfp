#include "cvs.h"

// global variables

enum {ZSDCIRCLE, 
      ZSDXRANGE, 
      ZSDRING, 
      RMSD, 
      LINE, 
      BOND,
      BONDS,
      HALFBOND,
      S,
      ANGLE,
      DIHED,
      COGX,
      COGY,
      COGZ,
      CARTESIAN_X,
      CARTESIAN_Y,
      CARTESIAN_Z, 
      NULL_CV};

char * CVSTRINGS[NULL_CV] = {
      "ZSDCIRCLE",
      "ZSDXRANGE",
      "ZSDRING",
      "RMSD",
      "LINE",
      "BOND",
      "BONDS",
      "HALFBOND",
      "S",
      "ANGLE",
      "DIHED",
      "COGX",
      "COGY",
      "COGZ",
      "CARTESIAN_X",
      "CARTESIAN_Y",
      "CARTESIAN_Z"};

// zsd circle
double zsdc_x,zsdc_y;
double zsdc_d; //radius of the cilinder
double zsdc_s; //radius of the lipids
//cv * zsdc_c=NULL;

// zsd ring
double zsdr_x,zsdr_y;
double zsdr_r1; //internal radius of the cilinder
double zsdr_r2; //external radius of the cilinder
double zsdr_s; //radius of the lipids


// Constructores

cv * cv_init ( char * typ, int nC, int * ind, 
    double zmin, double zmax, char * boundstr, double boundk,
    char * outfile, int outputFreq )  {
  int i;
  cv * c=malloc(sizeof(cv));

  c->nC=nC;
  c->val=0.0;

  c->ind=calloc(nC,sizeof(int));
  for (i=0;i<nC;i++) c->ind[i]=ind[i];

  c->gr=(double**)malloc(nC*sizeof(double*));
  for (i=0;i<nC;i++) c->gr[i]=(double*)malloc(3*sizeof(double));

  c->typ=cv_getityp(typ);
  switch(c->typ) {
    case ZSDCIRCLE:   c->calc = calccv_zsd_circle; break;
    case ZSDXRANGE:   c->calc = calccv_zsd_xrange; break;
    case ZSDRING:     c->calc = calccv_zsd_ring; break;
    case RMSD:        c->calc = calccv_rmsd; 
                      c->ref = (double**)malloc(nC*sizeof(double*));
                      for (i=0;i<nC;i++) c->ref[i] = calloc(3,sizeof(double));
                      break;
    case LINE:    c->calc = calccv_line; 
                      c->ref = (double**)malloc(nC*sizeof(double*));
                      for (i=0;i<nC;i++) c->ref[i] = calloc(3,sizeof(double));
                      c->ref2 = (double**)malloc(nC*sizeof(double*));
                      for (i=0;i<nC;i++) c->ref2[i] = calloc(3,sizeof(double));
                      break;
    case BOND:        c->calc = calccv_bond; break;
    case BONDS:       c->calc = calccv_bonds; break;
    case HALFBOND:    c->calc = calccv_halfbond; break;
    case S:           c->calc = calccv_s; break;
    case ANGLE:       c->calc = calccv_angle; break;
    case DIHED:       c->calc = calccv_dihed; break;
    case COGX:	      c->calc = calccv_cogx; break;
    case COGY:	      c->calc = calccv_cogy; break;
    case COGZ:	      c->calc = calccv_cogz; break;
    case CARTESIAN_X: c->calc = calccv_x; break;
    case CARTESIAN_Y: c->calc = calccv_y; break;
    case CARTESIAN_Z: c->calc = calccv_z; break;
    default: 
      fprintf(stderr,"ERROR, CV not recognized");
      fflush(stderr);exit(-1);break;
  }
                       
  // boundary function
  c->f=0.;
  c->u=0.;
  c->min=zmin;
  c->max=zmax;
  c->half_domain=0.5*(zmax-zmin);
  c->boundk=boundk;
  if     (!strcmp(boundstr,"SOFTUPPER")) {c->boundFunc = cv_SoftUpperWall;}
  else if(!strcmp(boundstr,"SOFTLOWER")) {c->boundFunc = cv_SoftLowerWall;} 
  else if(!strcmp(boundstr,"SOFT"     )) {c->boundFunc = cv_SoftWalls    ;} 
  else if(!strcmp(boundstr,"NADA"     )) {c->boundFunc = cv_nada         ;} 
  else {
    fprintf(stderr, "Error: boundary type not recognized");
    exit(1);
  }
 

  // output
  c->boutput=(outputFreq>0);
  if (c->boutput) {
    c->boutput=1;
    c->outputFreq=outputFreq;
    c->ofp=fopen(outfile,"w");
  }
   
  return c;
}

// Get addresses
double * cv_access_ref ( cv * c, int i ) {
  if (!c)         {fprintf(stderr,"CVS) null argument\n"); exit(-1);}
  if (i>=c->nC)   {fprintf(stderr,"CVS) out of size\n"); exit(-1);}
  return c->ref[i];
}                   

double * cv_access_ref2 ( cv * c, int i ) {
  if (!c)         {fprintf(stderr,"CVS) null argument\n"); exit(-1);}
  if (i>=c->nC)   {fprintf(stderr,"CVS) out of size\n"); exit(-1);}
  return c->ref2[i];
}                   
 
// Get properties
int cv_getityp ( char * typ ) {
  int i;
  for (i=0;i<NULL_CV&&strcmp(typ,CVSTRINGS[i]);i++);
  if (i<NULL_CV) return i;
  else return -1;
}

char * cv_getstyp ( int ityp ) {
  if (ityp<NULL_CV) return CVSTRINGS[ityp];
  else return "NOT_FOUND";
}
 

// CV computing 

int calccv_rmsd (cv * c, double ** R ) {
  int l,i,d;

  double aux, sum2;

  //  optimal rotation ? 
  //
  //    E. A. Coutsias, C. Seok, and K. A. Dill.
  //    Using quaternions to calculate RMSD.
  //    J. Comput. Chem., 25(15):1849-1857, 2004. 
  //

  c->val=0.;
  for (l=0;l<c->nC;l++) {
    for (d=0;d<3;d++) {
      i=c->ind[l];
      aux=R[i][d]-c->ref[i][d];
      c->val+=aux*aux;
    }
  }

  for (l=0;l<c->nC;l++) {
    for (d=0;d<3;d++) {
      i=c->ind[l];
      c->gr[i][d]=2.0*(R[i][d]-c->ref[i][d]);
    }
  } 

  c->val=c->val/c->nC;

  return 0;

}

int set_line (cv * c) {
  int l,i,d;
  double aux;

  // // Corrects COG of ref
  // for (d=0;d<3;d++) {
  //   aux=0.;
  //   for (l=0;l<c->nC;l++) aux+=c->ref[l][d];
  //   for (l=0;l<c->nC;l++) c->ref[l][d]-=aux/c->nC;
  // }
  //
  // // Corrects COG of ref2
  // for (d=0;d<3;d++) {
  //   aux=0.;
  //   for (l=0;l<c->nC;l++) aux+=c->ref2[l][d];
  //   for (l=0;l<c->nC;l++) c->ref2[l][d]-=aux/c->nC;
  // }
  //
  //Optimal rotation ? 
  
  // Computes ref2-ref and store in gr
  c->refmod=0.;
  for (d=0;d<3;d++) {
    for (l=0;l<c->nC;l++) {
      c->gr[l][d]=c->ref2[l][d]-c->ref[l][d];
      c->refmod+=c->gr[l][d]*c->gr[l][d];
    }
  }

  // 0 is ref and 1 is ref2 so refmod go to the square
  // c->refmod=sqrt(c->refmod);

  // Normalize
  for (d=0;d<3;d++) {
    for (l=0;l<c->nC;l++) c->gr[l][d]=c->gr[l][d]/c->refmod;
  }
   
  return 0;

}

int calccv_line (cv * c, double ** R ) {
  int l,i,d;
  double cog[3];
  double aux;

  //Correct COG

  // // Computes COG of ref
  // for (d=0;d<3;d++) {
  //   aux=0.;
  //   for (l=0;l<c->nC;l++) {
  //     i=c->ind[l];
  //     aux+=R[i][d];
  //   }
  //   cog[d] = aux/c->nC;
  // }
                    
  // Computes cv value
  c->val=0.;
  for (l=0;l<c->nC;l++) {
    for (d=0;d<3;d++) {
      i=c->ind[l];
      // c->val+=(R[i][d]-cog[d]-c->ref[l][d])*c->gr[l][d];
      c->val+=(R[i][d]-c->ref[l][d])*c->gr[l][d];
    }
  }
               
  //The gradient is constant
     
  return 0;

}
 
int calccv_s ( cv * c, double ** R ) {
  /* S is the CV of bond networks (see \cite{Barducci2006}) */
  int j,k,l;
  double r,aux;

  c->val=0.;

  for (j=0;j<c->nC;j+=2) {
    k=c->ind[j];
    l=c->ind[j+1];
    r=my_getbond( R[k], R[l], c->gr[j], c->gr[j+1]);
    aux=pow((r/2.5),6);
    c->val+=1./(1.+aux);
    aux=(6.*aux/r)/((1.+aux)*(1.+aux));

    for (k=0;k<3;k++) {
      c->gr[j  ][k]=-aux*c->gr[j  ][k];
      c->gr[j+1][k]=-aux*c->gr[j+1][k];
    }
  }

  return 0;
}

int set_zsd_circle ( double x,double y, double xy, double s  ) {

  zsdc_x=x;
  zsdc_y=y;
  zsdc_d=xy; //radius of the cilinder
  zsdc_s=s; //radius of the lipids

  return 0;
}

int calccv_zsd_circle ( cv * c, double ** R ) {
  int i,l;
  double zsdc_v,zsdc_z;
  double aux,aux1,aux2,d,norm;


  // Compute the average z
  zsdc_z=0.;
  norm=0.;
  for (l=0;l<c->nC;l++) {
    i=c->ind[l];

    // Compute the weight
    aux1= (R[i][0]-zsdc_x);
    aux2= (R[i][1]-zsdc_y);
    d=sqrt(aux1*aux1+aux2*aux2);
    aux1= zsdc_d-d;
    aux2=-zsdc_d-d;
    aux=cdf(aux1/zsdc_s)-cdf(aux2/zsdc_s);

    // Accumulate and save the unnormalized weight
    norm+=aux;
    zsdc_z+=aux*R[i][2];
    c->gr[l][2]=aux;
 
    // Save the derivatvies of the unnormalized weight
    // (cdf(aux1/s))'=gauss(aux1,s,0)*(x-x0)/d
    aux1=exp(-aux1*aux1/(2.*zsdc_s*zsdc_s));
    aux2=exp(-aux2*aux2/(2.*zsdc_s*zsdc_s));

    // I need a versor for the distance but d could be zero...
    if (d<.1e-15) { 
      aux=0.;  
    } else {
     aux=-(aux1-aux2)/(zsdc_s*sqrt(2*M_PI)*d);
    }
    c->gr[l][0]=aux*(R[i][0]-zsdc_x);
    c->gr[l][1]=aux*(R[i][1]-zsdc_y);
  }
  zsdc_z=zsdc_z/norm;

  // Compute the variance of z
  zsdc_v=0.;
  for (l=0;l<c->nC;l++) {
    i=c->ind[l];

    // Accumulate
    aux=(R[i][2]-zsdc_z);
    zsdc_v+=aux*aux*c->gr[l][2];
  }

  zsdc_v=zsdc_v/norm;
  //fprintf(stdout,"BLP) %.5f %.5f\n",aux*aux,c->gr[l][2]);

  c->val=zsdc_v;

  // Force on each lipid
  for (l=0;l<c->nC;l++) {
    i=c->ind[l];

    aux=(R[i][2]-zsdc_z);

    //c->gr[l][0]=c->gr[l][0]*(aux*aux-zsdc_v)/norm;
    //c->gr[l][1]=c->gr[l][1]*(aux*aux-zsdc_v)/norm;
    c->gr[l][0]=c->gr[l][0]*(aux*aux-3*zsdc_v)/norm;
    c->gr[l][1]=c->gr[l][1]*(aux*aux-3*zsdc_v)/norm;
    c->gr[l][2]=2*c->gr[l][2]*aux/norm;

  }

  return 0;

}

int calccv_zsd_xrange ( cv * c, double ** R ) {
  // set_zsd_cricle set enough number of parameters for this cv
  // therefore we will use:
  //   zsdc_d as the widht of the x range
  //   zsdc_x as the position of the x range
  //   zsdc_y discarded
  int i,l;
  double zsdc_v,zsdc_z;
  double aux,aux1,aux2,d,norm;


  // Compute the average z
  zsdc_z=0.;
  norm=0.;
  for (l=0;l<c->nC;l++) {
    i=c->ind[l];

    // Compute the weight
    d= abs(R[i][0]-zsdc_x);
    aux1= zsdc_d-d;
    aux2=-zsdc_d-d;
    aux=cdf(aux1/zsdc_s)-cdf(aux2/zsdc_s);

    // Accumulate and save the unnormalized weight
    norm+=aux;
    zsdc_z+=aux*R[i][2];
    c->gr[l][2]=aux;
 
    // Save the derivatvies of the unnormalized weight
    // (cdf(aux1/s))'=gauss(aux1,s,0)*(x-x0)/d
    aux1=exp(-aux1*aux1/(2.*zsdc_s*zsdc_s));
    aux2=exp(-aux2*aux2/(2.*zsdc_s*zsdc_s));

    // I need a versor for the distance but d could be zero...
    if (d<.1e-15) { 
      aux=0.;  
    } else {
      aux=-(aux1-aux2)/(zsdc_s*sqrt(2*M_PI)*d);
    }
    c->gr[l][0]=aux*(R[i][0]-zsdc_x);
    c->gr[l][1]=0.0;
  }
  zsdc_z=zsdc_z/norm;

  // Compute the variance of z
  zsdc_v=0.;
  for (l=0;l<c->nC;l++) {
    i=c->ind[l];

    // Accumulate
    aux=(R[i][2]-zsdc_z);
    zsdc_v+=aux*aux*c->gr[l][2];
  }

  zsdc_v=zsdc_v/norm;
  //fprintf(stdout,"BLP) %.5f %.5f\n",aux*aux,c->gr[l][2]);

  c->val=zsdc_v;

  // Force on each lipid
  for (l=0;l<c->nC;l++) {
    i=c->ind[l];

    aux=(R[i][2]-zsdc_z);

    //c->gr[l][0]=c->gr[l][0]*(aux*aux-zsdc_v)/norm;
    //c->gr[l][1]=c->gr[l][1]*(aux*aux-zsdc_v)/norm;
    c->gr[l][0]=c->gr[l][0]*(aux*aux-3*zsdc_v)/norm;
    c->gr[l][1]=0.0;
    c->gr[l][2]=2*c->gr[l][2]*aux/norm;

  }

  return 0;

}
 
int set_zsd_ring ( double x,double y, double r1, double r2, double s  ) {

  zsdr_x=x;
  zsdr_y=y;
  zsdr_r1=r1;
  zsdr_r2=r2;
  zsdr_s=s; 

  return 0;
}

int calccv_zsd_ring ( cv * c, double ** R ) {
  int i,l;
  double zsdr_v,zsdr_z;
  double aux,aux1,aux2,aux3,aux4,d,norm;


  // Compute the average z
  zsdr_z=0.;
  norm=0.;
  for (l=0;l<c->nC;l++) {
    i=c->ind[l];

    // Compute the weight
    aux1= (R[i][0]-zsdr_x);
    aux2= (R[i][1]-zsdr_y);
    d=sqrt(aux1*aux1+aux2*aux2);
    aux1= zsdr_r2-d;
    aux2= zsdr_r1-d;
    aux=cdf(aux1/zsdr_s)-cdf(aux2/zsdr_s);
    aux3= -zsdr_r1-d;
    aux4= -zsdr_r2-d;
    aux+=cdf(aux3/zsdr_s)-cdf(aux4/zsdr_s);

    // Accumulate and save the unnormalized weight
    norm+=aux;
    zsdr_z+=aux*R[i][2];
    c->gr[l][2]=aux;
 
    // Save the derivatvies of the unnormalized weight
    // (cdf(aux1/s))'=gauss(aux1,s,0)*(x-x0)/d
    aux1=exp(-aux1*aux1/(2.*zsdr_s*zsdr_s));
    aux2=exp(-aux2*aux2/(2.*zsdr_s*zsdr_s));
    aux3=exp(-aux3*aux3/(2.*zsdr_s*zsdr_s));
    aux4=exp(-aux4*aux4/(2.*zsdr_s*zsdr_s));

    // d can't be zero here
    //// I need a versor for the distance but d could be zero...
    //if (d<.1e-15) { 
    //  aux=0.;  
    //} else {
    //  aux=-(aux1-aux2)/(zsdc_s*sqrt(2*M_PI)*d);
    aux=-(aux1-aux2-aux1+aux2)/(zsdr_s*sqrt(2*M_PI)*d);
    //} 
    c->gr[l][0]=aux*(R[i][0]-zsdr_x);
    c->gr[l][1]=aux*(R[i][1]-zsdr_y);
  }
  zsdr_z=zsdr_z/norm;

  // Compute the variance of z
  zsdr_v=0.;
  for (l=0;l<c->nC;l++) {
    i=c->ind[l];

    // Accumulate
    aux=(R[i][2]-zsdr_z);
    zsdr_v+=aux*aux*c->gr[l][2];
  }

  zsdr_v=zsdr_v/norm;
  //fprintf(stdout,"BLP) %.5f %.5f\n",aux*aux,c->gr[l][2]);

  c->val=zsdr_v;

  // Force on each lipid
  for (l=0;l<c->nC;l++) {
    i=c->ind[l];

    aux=(R[i][2]-zsdr_z);

    //c->gr[l][0]=c->gr[l][0]*(aux*aux-zsdr_v)/norm;
    //c->gr[l][1]=c->gr[l][1]*(aux*aux-zsdr_v)/norm;
    c->gr[l][0]=c->gr[l][0]*(aux*aux-3*zsdr_v)/norm;
    c->gr[l][1]=c->gr[l][1]*(aux*aux-3*zsdr_v)/norm;
    c->gr[l][2]=2*c->gr[l][2]*aux/norm;

  }

  return 0;

}
     
int calccv_bonds ( cv * c, double ** R ) {
  /* BONDS is a bond network not normalized */
  int j,k,l;
  double r,aux;

  c->val=0.;
  for (j=0;j<c->nC;j+=2) {
    k=c->ind[j];
    l=c->ind[j+1];
    r=my_getbond( R[k], R[l], c->gr[j], c->gr[j+1]);
    c->val+=r;
  }

  return 0;
}
     
int calccv_bond ( cv * c, double ** R ) {
  c->val=my_getbond(R[c->ind[0]],R[c->ind[1]],c->gr[0],c->gr[1]);
  return 0;
}
 
int calccv_halfbond ( cv * c, double ** R ) {
  c->val=my_getbond(R[c->ind[0]],R[c->ind[1]],c->gr[0],c->gr[1]);
  c->gr[1][0]=0.;
  c->gr[1][1]=0.;
  c->gr[1][2]=0.;
  return 0;
}
 
int calccv_angle ( cv * c, double ** R ) {
  c->val=my_getangle(R[c->ind[0]],R[c->ind[1]],R[c->ind[2]],
			    c->gr[0],        c->gr[1],        c->gr[2]);
  return 0;
}

int calccv_dihed ( cv * c, double ** R ) {
  c->val=my_getdihed(R[c->ind[0]],R[c->ind[1]],R[c->ind[2]],R[c->ind[3]],
  		             c->gr[0],        c->gr[1],       c->gr[2],        c->gr[3]);
#ifdef _PARANOIA_
	if (_PARANOIA_) {
	  if (c->val!=c->val) {
	    fprintf(stderr,"CVS/C/PARANOIA) Tripped at dihed cvi->val %.5f\n",cvi->val);
	    fprintf(stderr,"Program exits.\n");
	    fflush(stderr);
	    exit(-1);
	  }
	}
#endif
  return 0;
}

int calccv_x ( cv * c, double ** R ) {
  c->val=R[c->ind[0]][0];
  c->gr[0][0]=1.0;
  c->gr[0][1]=0.0;
  c->gr[0][2]=0.0; 
  return 0;
}

int calccv_y ( cv * c, double ** R ) {
  c->val=R[c->ind[0]][1];
  c->gr[0][0]=0.0;
  c->gr[0][1]=1.0;
  c->gr[0][2]=0.0; 
  return 0;
}

int calccv_z ( cv * c, double ** R ) {
  c->val=R[c->ind[0]][2];
  c->gr[0][0]=0.0;
  c->gr[0][1]=0.0;
  c->gr[0][2]=1.0; 
  return 0;
}

int calccv_cogx ( cv * c, double ** R ) {
  double aux,aux2;
  int l,i;

  aux=0.;
  aux2=1./c->nC;
  for (l=0;l<c->nC;l++) {
    i=c->ind[l];
    c->gr[l][0]=aux2;
    c->gr[l][1]=0.;
    c->gr[l][2]=0.;
    aux+=R[i][0];
  }

  c->val = aux/c->nC;
  return 0;
}

int calccv_cogy ( cv * c, double ** R ) {
  double aux,aux2;
  int l,i;

  aux=0.;
  aux2=1./c->nC;
  for (l=0;l<c->nC;l++) {
    i=c->ind[l];
    c->gr[l][0]=0.;
    c->gr[l][1]=aux2;
    c->gr[l][2]=0.;
    aux+=R[i][1];
  }

  c->val = aux/c->nC;
  return 0;
}
 
int calccv_cogz ( cv * c, double ** R ) {
  double aux,aux2;
  int l,i;

  aux=0.;
  aux2=1./c->nC;
  for (l=0;l<c->nC;l++) {
    i=c->ind[l];
    c->gr[l][0]=0.;
    c->gr[l][1]=0.;
    c->gr[l][2]=aux2;
    aux+=R[i][2];
  }  

  c->val = aux/c->nC;
  return 0;
}
 

double cdf(double x)
//cumulative density function (CDF) of a standard normal (Gaussian) random variable
//cdf, thanks to John D. Cook. http://www.johndcook.com/blog/cpp_phi/
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
 
    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);
 
    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
 
    return 0.5*(1.0 + sign*y);
}


//double zsdc_cut(double r,double cu,double cl,double b,double dzsdc_cut){
//  double r1,r2,r3,r4,r5;
//  
//  r1=cu+d/2
//
//
//  if(r<r0) {
//    fc=1.0;
//    dfc=0.0;
//  } else if(r>r2) {
//    fc=0.0;
//    dfc=0.0;
//  } else {
//    aux=M_PI*(r-(r1+r2)*0.5)/(r2-r1);
//    fc =0.5-0.5*sin(aux);
//    dfc=-cos(aux)*M_PI*r/(r2-r1)*0.5;
//  }
//
//  return fc;
//}    
//double fcut(double r,double r1,double r2,double dfcut){
//  double fc,dfc,aux;
//  
//  if(r<r1) {
//    fc=1.0;
//    dfc=0.0;
//  } else if(r>r2) {
//    fc=0.0;
//    dfc=0.0;
//  } else {
//    aux=M_PI*(r-(r1+r2)*0.5)/(r2-r1);
//    fc =0.5-0.5*sin(aux);
//    dfc=-cos(aux)*M_PI*r/(r2-r1)*0.5;
//  }
//
//  return fc;
//}
//

// BOUNDARIES

int cv_nada ( cv * c ) {
  return 0;
}

int cv_SoftWalls ( cv * c ) {
  cv_SoftUpperWall(c);
  cv_SoftLowerWall(c);
}
 
int cv_SoftLowerWall ( cv * c ) {
  double aux;
  aux=c->val-c->min;
  if (aux>0.) return 0;
  c->f-=c->boundk*aux;
  c->u+=.5*c->boundk*aux*aux;
  return 0;
}

int cv_SoftUpperWall ( cv * c ) {
  double aux;
  aux=c->val-c->max;
  if (aux<0.) return 0;
  c->f-=c->boundk*aux;
  c->u+=.5*c->boundk*aux*aux;
  return 0;
}

void cv_output ( cv * c ) {
  fprintf(c->ofp,"%11.5f",c->val);
  fprintf(c->ofp," %11.5f\n",c->f);
  fflush(c->ofp);
}
            
