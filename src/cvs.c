
#include "cvs.h"

// global variables
enum {BOND, ANGLE, DIHED, CARTESIAN_X, CARTESIAN_Y, CARTESIAN_Z, S,BILAYP, NULL_CV};
char * CVSTRINGS[NULL_CV] = {"BOND", "ANGLE", "DIHED", "CARTESIAN_X", "CARTESIAN_Y", "CARTESIAN_Z", "S","BILAYP"};
double blpx,blpy,blpz;
double blpd=3.; //diameter of the cilinder
 
int cv_dimension ( cvStruct * c ) {
  // I shuld remove this
  int d;
  d=0;
  switch(c->typ) {
    case CARTESIAN_X: d=0; break;
    case CARTESIAN_Y: d=1; break;
    case CARTESIAN_Z: d=2; break;
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
    case S:           c->calc = calccv_s; break;
    case BOND:        c->calc = calccv_bond; break;
    case DIHED:       c->calc = calccv_dihed; break;
    case ANGLE:       c->calc = calccv_angle; break;
    case BILAYP:      c->calc = calccv_bilayerpoint; break;
  }

  c->nC=nC;
  c->val=0.0;

  c->ind=calloc(nC,sizeof(int));
  if (ind) for (i=0;i<nC;i++) c->ind[i]=ind[i];

  c->gr=(double**)malloc(nC*sizeof(double*));
  for (i=0;i<nC;i++) c->gr[i]=(double*)malloc(3*sizeof(double));

  return c;
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
  double aux,cl,cu,d;
  double r[3];

  cu=0.; //upper position
  cl=0.; //lower position
  c->val=0.;

  // loop over each lipid
  for (i=0;i<c->nC;i+=3) {

    // Compute the COM of the lipid
    for (k=0;i<3.;k++) {
      r[k]+=ds->R[i  ][k];
      r[k]+=ds->R[i+1][k];
      r[k]+=ds->R[i+2][k];
      r[k]=r[k]/3.;
    }

    //if (pbc) {}

    d =(r[0]-blpx)*(r[0]-blpx)
      +(r[1]-blpy)*(r[1]-blpy);

    
    aux=(cdf(d)-cdf(d-blpd));
    

    if(r[2]>blpz) {
      cu=cu+r[2]*aux;
      c->gr[i  ][2]+=aux;
      c->gr[i+1][2]+=aux;
      c->gr[i+2][2]+=aux;
    } else {
      cl=cl+r[2]*aux;
      c->gr[i  ][2]-=aux;
      c->gr[i+1][2]-=aux;
      c->gr[i+2][2]-=aux;
    }

  }

  c->val=cu-cl;

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
   
