#include "cvs.h"

// global variables

enum {BOND, ANGLE, DIHED, CARTESIAN_X, CARTESIAN_Y, CARTESIAN_Z, S,BILAYP,COGX, COGY, COGZ, NULL_CV};
char * CVSTRINGS[NULL_CV] = {
  "BOND", "ANGLE", "DIHED", 
  "CARTESIAN_X", "CARTESIAN_Y", "CARTESIAN_Z", 
  "COGX", "COGY", "COGZ", 
  "S","BILAYP"};

// bylayer
double blpx,blpy,blpz;
double blpdo2,blpd2; //diameter of the cilinder
//cvStruct * blpc=NULL;


int cv_dimension ( cvStruct * c ) {
  // I shuld remove this
  int d;
  d=0;
  switch(c->typ) {
    case CARTESIAN_X: d=0; break;
    case CARTESIAN_Y: d=1; break;
    case CARTESIAN_Z: d=2; break;
    case COGX: d=0; break;
    case COGY: d=1; break;
    case COGZ: d=2; break;
  }
  return d;
}
 
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
 
cvStruct * New_cvStruct ( int typ, int nC, int * ind ) {
  int i;
  cvStruct * c=malloc(sizeof(cvStruct));

  c->typ=typ;
  c->calc = calccv_s;
  switch(c->typ) {
    case CARTESIAN_X: c->calc = calccv_x; break;
    case CARTESIAN_Y: c->calc = calccv_y; break;
    case CARTESIAN_Z: c->calc = calccv_z; break;
    case COGX: c->calc = calccv_cogx; break;
    case COGY: c->calc = calccv_cogy; break;
    case COGZ: c->calc = calccv_cogz; break;
    case S:           c->calc = calccv_s; break;
    case BOND:        c->calc = calccv_bond; break;
    case DIHED:       c->calc = calccv_dihed; break;
    case ANGLE:       c->calc = calccv_angle; break;
    case BILAYP:
      //if (blpc) {
      //  fprintf(stderr, "Error: only 1 bilayerp allowed for now");
      //  exit(1);
      //}
      c->calc = calccv_bilayerpoint;
      //blpc=c;
      break;
  }

  c->nC=nC;
  c->val=0.0;

  c->ind=calloc(nC,sizeof(int));
  if (ind) for (i=0;i<nC;i++) c->ind[i]=ind[i];

  c->gr=(double**)malloc(nC*sizeof(double*));
  for (i=0;i<nC;i++) c->gr[i]=(double*)malloc(3*sizeof(double));

  return c;
}


int set_bilayerpoint ( double x,double y, double xy ) {

  blpx=x;
  blpy=y;
  blpdo2=xy*.5; //diameter square of the cilinder
  blpd2=xy*xy;

  return 0;
}

int calccv_s ( cvStruct * c, DataSpace * ds ) {
  /* S is the CV of bond networks (see \cite{Barducci2006}) */
  int j,k,l;
  double r,aux;

  c->val=0.;

  for (j=0;j<c->nC;j+=2) {
    k=c->ind[j];
    l=c->ind[j+1];
    r=my_getbond( ds->R[k], ds->R[l], c->gr[j], c->gr[j+1]);
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

int calccv_bilayerpoint ( cvStruct * c, DataSpace * ds ) {
  int i,j,k,l;
  double aux,cl,cu,d,normu,norml;

  blpz=0.;

  j=0;

  // loop over each lipid
  for (l=0;l<c->nC;l++) {
    i=c->ind[l];
    c->gr[i][0]=0.;
    c->gr[i][1]=0.;
    c->gr[i][2]=0.;

    // Compute the distance to the point
    d =(ds->R[i][0]-blpx)*(ds->R[i][0]-blpx)
      +(ds->R[i][1]-blpy)*(ds->R[i][1]-blpy);
    
    if(d-blpd2<0.) {
      blpz+=ds->R[i][2];
      j++;
    }

  }
  
  blpz=blpz/j;

  cu=0.; //upper position
  cl=0.; //lower position
  normu=0.;
  norml=0.;
 


  // Compute the CV and the normfactor in each layer
  for (l=0;l<c->nC;l++) {
    i=c->ind[l];

    // Compute the distance to the point
    d =(ds->R[i][0]-blpx)*(ds->R[i][0]-blpx)
      +(ds->R[i][1]-blpy)*(ds->R[i][1]-blpy);

    aux=sqrt(d);
    aux=cdf(-aux+blpdo2)-cdf(-aux-blpdo2);

    if(ds->R[i][2]>blpz) {
      cu+=ds->R[i][2]*aux;
      c->gr[i][2]=aux;
      normu+=aux;
    } else {
      cl+=ds->R[i][2]*aux;
      c->gr[i][2]=-aux;
      norml+=aux;
    }

  }

  c->val=cu/normu-cl/norml;

  
  // Normalize the force on each lipid
  for (l=0;l<c->nC;l++) {
    i=c->ind[l];

    if(ds->R[i][2]>blpz) {
      c->gr[i][2]=c->gr[i][2]/normu;
    } else {
      c->gr[i][2]=c->gr[i][2]/norml;
    }

    //fprintf(stderr,"aa0 %f\n",c->gr[i][0]);
    //fprintf(stderr,"aa1 %f\n",c->gr[i][1]);
    //fprintf(stderr,"aa2 %f\n",c->gr[i][2]);

  }

  //fprintf(stderr,"aaaaaaaaaaaaa %f\n",cu/normu);
  //fprintf(stderr,"aaaaaaaaaaaaa %f\n",cl/norml);
  //fprintf(stderr,"aaaaaaaaaaaaa %f\n",c->val);
  //exit(1);

}

int calccv_bond ( cvStruct * c, DataSpace * ds ) {
  c->val=my_getbond(ds->R[c->ind[0]],ds->R[c->ind[1]],c->gr[0],c->gr[1]);
  return 0;
}
 
int calccv_angle ( cvStruct * c, DataSpace * ds ) {
  c->val=my_getangle(ds->R[c->ind[0]],ds->R[c->ind[1]],ds->R[c->ind[2]],
			    c->gr[0],        c->gr[1],        c->gr[2]);
  return 0;
}

int calccv_dihed ( cvStruct * c, DataSpace * ds ) {
  c->val=my_getdihed(ds->R[c->ind[0]],ds->R[c->ind[1]],ds->R[c->ind[2]],ds->R[c->ind[3]],
  		             c->gr[0],        c->gr[1],       c->gr[2],        c->gr[3]);
#ifdef _PARANOIA_
	if (_PARANOIA_) {
	  if (c->val!=c->val) {
	    fprintf(stderr,"CFACV/C/PARANOIA) Tripped at dihed cvi->val %.5f\n",cvi->val);
	    fprintf(stderr,"Program exits.\n");
	    fflush(stderr);
	    exit(-1);
	  }
	}
#endif
  return 0;
}

int calccv_x ( cvStruct * c, DataSpace * ds ) {
  c->val=ds->R[c->ind[0]][0];
  c->gr[0][0]=1.0;
  c->gr[0][1]=0.0;
  c->gr[0][2]=0.0; 
  return 0;
}

int calccv_y ( cvStruct * c, DataSpace * ds ) {
  c->val=ds->R[c->ind[0]][1];
  c->gr[0][0]=0.0;
  c->gr[0][1]=1.0;
  c->gr[0][2]=0.0; 
  return 0;
}

int calccv_z ( cvStruct * c, DataSpace * ds ) {
  c->val=ds->R[c->ind[0]][2];
  c->gr[0][0]=0.0;
  c->gr[0][1]=0.0;
  c->gr[0][2]=1.0; 
  return 0;
}

int calccv_cogx ( cvStruct * c, DataSpace * ds ) {
  double aux,aux2;
  int l,i;

  aux=0.;
  aux2=1./c->nC;
  for (l=0;l<c->nC;l++) {
    i=c->ind[l];
    c->gr[i][0]=aux2;
    c->gr[i][1]=0.;
    c->gr[i][2]=0.;
    aux+=ds->R[i][0];
  }

  c->val = aux/c->nC;
  return 0;
}

int calccv_cogy ( cvStruct * c, DataSpace * ds ) {
  double aux,aux2;
  int l,i;

  aux=0.;
  aux2=1./c->nC;
  for (l=0;l<c->nC;l++) {
    i=c->ind[l];
    c->gr[i][0]=0.;
    c->gr[i][1]=aux2;
    c->gr[i][2]=0.;
    aux+=ds->R[i][1];
  }

  c->val = aux/c->nC;
  return 0;
}
 
int calccv_cogz ( cvStruct * c, DataSpace * ds ) {
  double aux,aux2;
  int l,i;

  aux=0.;
  aux2=1./c->nC;
  for (l=0;l<c->nC;l++) {
    i=c->ind[l];
    c->gr[i][0]=0.;
    c->gr[i][1]=0.;
    c->gr[i][2]=aux2;
    aux+=ds->R[i][2];
  }

  c->val = aux/c->nC;
  return 0;
}
 

double cdf(double x)
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


//double blpcut(double r,double cu,double cl,double b,double dblpcut){
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
