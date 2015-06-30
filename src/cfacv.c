#include "cfacv.h"

// Handling SIGSEGV or SIGINT signal to make a call stack using backtrace.
// The call stack can be used to profile the program
// See for instance:
// http://stackoverflow.com/a/378024/1342186
// http://stackoverflow.com/a/77336/1342186

#include <execinfo.h>
#include <signal.h>
#include <unistd.h>

//#define _PARANOIA_ 1

void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}

void cfacvBanner ( void ) {
  printf("CFACV/C) Version (Git Hash): %s\n",VERSION);
}

FILE * my_fopen ( char * name, char * code ) {
  if (strcmp(name,"stdout")) return fopen(name,code);
  else return stdout;
}

FILE * my_binfopen ( char * name, char * code, unsigned int outputLevel, DataSpace * ds ) {
  FILE * fp=fopen(name,code);
  fwrite(&outputLevel,sizeof(unsigned int),1,fp);
  fwrite(&(ds->iK),sizeof(int),1,fp);
  return fp;
}


char * RFSTRINGS[NULL_RF] = {"HARMONIC","HARMCUTO"};
int rf_getityp ( char * typ ) {
  int i;
  for (i=0;i<NULL_RF&&strcmp(typ,RFSTRINGS[i]);i++);
  if (i<NULL_RF) return i;
  else return -1;
}

char * BADCVSTRINGMESSAGE = "NOT FOUND";
char * rf_getstyp ( int ityp ) {
  if (ityp<NULL_RF) return RFSTRINGS[ityp];
  else return BADCVSTRINGMESSAGE;
}




// RESTRAIN STRUCTURES

restrStruct * New_restrStruct ( 
    double k, double z, int nCV, 
    double * cvc, char * rftypstr, 
    double zmin, double zmax,
    char * boundstr, double boundk
    ) {

  restrStruct * newr=malloc(sizeof(restrStruct));
  int i,b_pbc,b_peri;
  b_pbc=0;
  b_peri=0;

  newr->rfityp=rf_getityp(rftypstr);
  newr->k=k;
  newr->z=z;
  newr->u=newr->f=0.0;
  newr->nCV=nCV;
  newr->cvc=(double*)malloc(nCV*sizeof(double));
  if (cvc) for (i=0;i<nCV;i++) newr->cvc[i]=cvc[i];
  else for (i=0;i<nCV;i++) newr->cvc[i]=0;

  newr->min=zmin;
  newr->max=zmax;
  newr->half_domain=0.5*(zmax-zmin);
  
  // evolve function
  newr->tamd_noise=0.0;
  newr->tamd_restraint=0.0;
  newr->tamdOpt=NULL;
  newr->smdOpt=NULL;
  newr->evolve=1;
  newr->evolveFunc = NULL;

  // boundary function
  newr->boundk=boundk;
  if     (!strcmp(boundstr,"SOFTUPPER")) {newr->boundFunc = SoftUpperWall      ;}
  else if(!strcmp(boundstr,"SOFTLOWER")) {newr->boundFunc = SoftLowerWall      ;} 
  else if(!strcmp(boundstr,"SOFT"     )) {newr->boundFunc = SoftWalls          ;} 
  else if(!strcmp(boundstr,"PBC"))       {newr->boundFunc = pbc; b_pbc=1       ;} 
  else if(!strcmp(boundstr,"PERIODIC" )) {newr->boundFunc = Periodic; b_peri=1 ;} 
  else if(!strcmp(boundstr,"NADA"     )) {newr->boundFunc = nada               ;} 
  else {
    fprintf(stderr, "Error: boundary type not recognized");
    exit(1);
  }

  // energy function
  if     (b_pbc) {
    if      (newr->rfityp==HARMONIC) {newr->energyFunc = HarmonicCart_pbc;}
    else if (newr->rfityp==HARMCUTO) {newr->energyFunc = HarmonicCart_cutoff_pbc;}
  } else if (b_peri) {
    //TODO
    //if      (newr->rfityp==HARMONIC) {newr->energyFunc = HarmonicCart_periodic;}
    //else if (newr->rfityp==HARMCUTO) {newr->energyFunc = HarmonicCart_cutoff_periodic;}
  } else {
    if      (newr->rfityp==HARMONIC) {newr->energyFunc = HarmonicCart;}
    else if (newr->rfityp==HARMCUTO) {newr->energyFunc = HarmonicCart_cutoff;}
  }


  return newr;
}

tamdOptStruct * New_tamdOptStruct ( double g, double kt, double dt, int riftyp) {
  tamdOptStruct * tamd=malloc(sizeof(tamdOptStruct));
  tamd->kT=kt;
  tamd->gamma=g;
  tamd->dt=dt;
  if (g>0.0) tamd->ginv=1.0/g;
  else tamd->ginv=0.0;
  tamd->noise = sqrt(2.0*kt*dt*tamd->ginv);
  return tamd;
}

smdOptStruct * New_smdOptStruct ( double target, int t0, int t1 ) {
  smdOptStruct * smd=malloc(sizeof(smdOptStruct));
  smd->t0=t0;
  smd->t1=t1;
  if (!(t1-t0)) smd->invinterval=0.;
  else smd->invinterval = 1.0/(t1-t0);
  smd->target=target;
  return smd;
}

int smdOptInit ( smdOptStruct * smd, double initval) {
  double dist=0.0;
  if (!smd) return -1;
  
  smd->initval=initval;
  dist = smd->target-smd->initval;
  //FIXME: Inquire if the restrain is linked to periodic instead.
  //if (periodic) {
  //  if (dist<-M_PI) {
  //    dist+=2*M_PI;
  //  }
  //  else if (dist > M_PI) {
  //    dist-=2*M_PI;
  //  } 
  //}
  smd->increment=dist*smd->invinterval;
  printf("CFACV/C) smd increment is %.5f based on target %.5f and initval %.5f\n",
         smd->increment,smd->target,smd->initval);
  return 0;
  
}

void ds_saverestrains ( DataSpace * ds, int timestep ) {
  int i;
  restrStruct * r;

  if (!(timestep%ds->restrsavefreq)) {
    rewind(ds->restrofs);
    for (i=0;i<ds->iK;i++) {
      r=ds->restr[i];
      fwrite(&(r->z),sizeof(double),1,ds->restrofs);
      fwrite(&(r->val),sizeof(double),1,ds->restrofs);
    }
    fflush(ds->restrofs);
  }
}


//ENERGY FUNCTIONS

int HarmonicCart ( restrStruct * r ) {
  double tmp=r->val-r->z;
  r->f=-r->k*tmp;
  r->u=0.5*r->k*tmp*tmp;
  return 0;
}

int HarmonicCart_cutoff ( restrStruct * r ) {
  //double r1=10.0,r2=12.0;
 
  // Non tie
  if (r->z>r->max) {
    r->f=0.0;
    r->u=0.0;
    r->evolve=0;
    return 0;
  } 


  //// Set up initial position above the CV
  //if(!r->evolve) {
  //  r->z=r->val;
  //  r->evolve=1;
  //  r->f=0.0;
  //  r->u=0.0;
  //  return 0;
  //}

  // Set up initial acording to a distribution around the CV
  if(!r->evolve) {
    r->z=r->val+sqrt(r->tamdOpt->kT/r->k)*my_whitenoise(Xi);
    r->evolve=1;
  }

  // Abrupt cutoff
  HarmonicCart(r);

  //// Force an energy with a smooth cutoff
  //if(r->z<r1) {
  //  aux=0.0;
  //  aux1=1.0;
  //} else {
  //tmp==r->val-r->z;
  //  if (r->pbc) {
  //    if (tmp>r->half_domain) tmp-=2*r->half_domain;
  //    if (tmp<-r->half_domain) tmp+=2*r->half_domain;
  //  }
  //  aux=M_PI*(tmp-(r1+r2)*0.5)/(r2-r1);
  //  aux1=0.5-0.5*sin(aux);
  //  aux=-cos(aux)*M_PI*tmp/(r2-r1)*0.5;
  //  //dfcut=-0.5_dp*cos(aux)*pio2*tmp/(r2-r1);
  //}
  //aux2=0.5*r->k*tmp*tmp;
  //r->f=-r->k*tmp*aux1-aux2*aux;
  //r->u=aux2*aux1;

  return 0;
}

int HarmonicCart_pbc ( restrStruct * r ) {
  double tmp=r->val-r->z;
  if (tmp>r->half_domain) {tmp-=2*r->half_domain;}
  if (tmp<-r->half_domain) {tmp+=2*r->half_domain;}
  r->f=-r->k*tmp;
  r->u=0.5*r->k*tmp*tmp;
  return 0;
}

int HarmonicCart_cutoff_pbc ( restrStruct * r ) {

  // Release tie
  if (r->z>r->max) {
    r->f=0.0;
    r->u=0.0;
    r->evolve=0;
    return 0;
  } 

  // Set up initial acording to a distribution around the CV
  if(!r->evolve) {
    r->z=r->val+sqrt(r->tamdOpt->kT/r->k)*my_whitenoise(Xi);
    r->evolve=1;
  }

  // Abrupt cutoff
  HarmonicCart_pbc(r);

  return 0;
}
      



//BOUNDARIES
 
int nada ( restrStruct * r ) {
}

int SoftWalls ( restrStruct * r ) {
  SoftUpperWall(r);
  SoftLowerWall(r);
}
 
int SoftLowerWall ( restrStruct * r ) {
  double aux;
  aux=r->z-r->min;
  if (aux>0.) return 0;
  // r->f is minus the force on z!
  r->f+=r->boundk*aux;
  r->u+=.5*r->boundk*aux*aux;
  return 0;
}

int SoftUpperWall ( restrStruct * r ) {
  double aux;
  aux=r->z-r->max;
  if (aux<0.) return 0;
  // r->f is minus the force on z!
  r->f+=r->boundk*aux;
  r->u+=.5*r->boundk*aux*aux;
  return 0;
}
 
int Periodic ( restrStruct * r ) {
  // TODO: Periodic is the same that PBC and can be unified.
  if (r->z>M_PI) r->z-=2*M_PI;
  else if (r->z<-M_PI) r->z+=2*M_PI;
  return 0;
}
 
int pbc ( restrStruct * r ) {
  if (r->z>r->max) r->z-=2*r->half_domain;
  if (r->z<-r->min) r->z+=2*r->half_domain;
  return 0;
}


// EVOLUTION OF RESTRIANS

/* Brownian dynamics (TAMD)*/
int cbd ( restrStruct * r, double f ) {
  tamdOptStruct * tamd=r->tamdOpt;
  double dd = tamd->ginv*tamd->dt*f;
  double rd = tamd->noise*my_whitenoise(Xi);

  if (!r->evolve) return 0;
  r->tamd_noise=rd;
  r->tamd_restraint=dd;
  r->z=r->z+dd+rd;
  return 0;
}

/* Constant Velocity (SMD)*/
int uniformvelocity ( restrStruct * r, double f ) {
  if (!r->evolve) return 0;
  r->z=r->z+f;
  return 0;
}

//double smdUpdate ( double z, double increment, int OK ) {
//  return z+OK*increment;
//}
//
//double smdUpdate_Periodic ( double z, double increment, int OK ) {
//  double newz=z+increment;
//  if (newz>M_PI) newz-=2*OK*M_PI;
//  else if (newz<-M_PI) newz+=2*OK*M_PI;
//  return newz;
//}







// DATA SPACE


DataSpace * NewDataSpace ( int N, int M, int K, long int seed ) {
  DataSpace * ds = malloc(sizeof(DataSpace));
  int i;
  ds->N=N;
  ds->M=M;
  ds->K=K;
  ds->iN=ds->iM=ds->iK=0;

  // Random seed. I deatached this from the data space a do it global
  // for all the file (see for example HarmonicCart_cutoff function)
  Xi=(unsigned short *)malloc(3*sizeof(unsigned short));
  Xi[0]=(seed & 0xffff00000000) >> 32;
  Xi[1]=(seed & 0x0000ffff0000) >> 16;
  Xi[2]= 0x330e;

  ds->R=(double**)malloc(N*sizeof(double*));
  for (i=0;i<N;i++) ds->R[i]=calloc(3,sizeof(double));

  ds->ac = (atomCenterStruct**)malloc(N*sizeof(atomCenterStruct*));

  ds->cv=(cvStruct**)malloc(M*sizeof(cvStruct*));
  ds->restr=(restrStruct**)malloc(K*sizeof(restrStruct*));

  ds->doAnalyticalCalc=0;
  ds->O[0]=ds->O[1]=ds->O[2]=0.0;
  ds->L[0]=ds->L[1]=ds->L[2]=0.0;
  ds->hL[0]=ds->hL[1]=ds->hL[2]=0.0;
  ds->evolveAnalyticalParameters=0;
  ds->beginaccum=0;

  return ds;
}

int DataSpace_SetupPBC ( DataSpace * ds, int pbc, double Ox, double Oy, double Oz, \
			      double Lx, double Ly, double Lz ) {
  ds->pbc=pbc;
  ds->O[0]=Ox;
  ds->O[1]=Oy;
  ds->O[2]=Oz;
  ds->L[0]=Lx;  ds->hL[0]=Lx/2.;
  ds->L[1]=Ly;  ds->hL[1]=Ly/2.;
  ds->L[2]=Lz;  ds->hL[2]=Lz/2.;
  ds->Max[0]=ds->L[0]-ds->O[0]; ds->Min[0]=ds->Max[0]-ds->L[0];
  ds->Max[1]=ds->L[1]-ds->O[1]; ds->Min[1]=ds->Max[1]-ds->L[1];
  ds->Max[2]=ds->L[2]-ds->O[2]; ds->Min[2]=ds->Max[2]-ds->L[2];
 
  return 0;
}

chapeau * DataSpace_get_chapeauadress ( DataSpace * ds, int i ) {
  return ds->ch[i];
}
  
int DataSpace_Setup1Dchapeau ( DataSpace * ds, int numrep, double min, int nKnots,
    double max, int beginaccum, int useTAMDforces, char * outfile, 
    int outfreq, int outlevel, int nupdate) {
  int i,j;
  int N=ds->N; // number of centers
  char buf[80];

  fprintf(stdout,"CFACV/C) DEBUG DataSpace_Setup %i\n",N);fflush(stdout);
  ds->doAnalyticalCalc=1;
  ds->useTAMDforces=useTAMDforces;
  if (ds->useTAMDforces) fprintf(stdout,"CFACV/C) INFO: using TAMD forces instead of analytical forces.\n");

  ds->evolveAnalyticalParameters=1;

  // This is to let the system equilibrate
  ds->beginaccum=beginaccum;

  //FIXME: This code was here to allow compute chapeau functions separatedly
  //for different pair types of particles. For instance, this allow to
  //recover SOD SOD, CLA CLA and SOD CLA pair potentials in 1 TAMD
  //simulation. Each index has a number in ch_id which allow to sort the pair
  //in the different chapeau objects on the c code.  From the studies with
  //SOD CLA, this pair potentials will be OK only if the ficticius
  //temperature is the same that the real one.  On the other hand, a better
  //way to achive this is needed (without saving a lot of numbers in ch_id).
  //For understand how this worked, see the previous versions of the code. 
  //i=N*(N-1)/2; // The upper triangular half of a rectangular matrix
  //ds->ch_id=(int*)malloc(i*sizeof(int));
  //ds->ch_num=chnum;

  ds->ch_num=numrep;
  ds->ch=(chapeau**)malloc(ds->ch_num*sizeof(chapeau*));

  for (i=0;i<ds->ch_num;i++) {
    ds->ch[i]=chapeau_alloc(nKnots,min,max,N);

    //No needed, calloc in chapeau_alloc
    //gsl_matrix_set_zero(ds->ch[i]->A);
    //gsl_vector_set_zero(ds->ch[i]->b);
    sprintf(buf, "%s_ch%d.bsp", outfile,i);
    chapeau_setupoutput(ds->ch[i],buf,outfreq,outlevel);
    ds->ch[i]->nupdate=nupdate;
  }

  ds->restrsavefreq=outfreq;

  //TODO: this two lines should be before the restart
  sprintf(buf, "%s_restr.restart",outfile);
  ds->restrofs=fopen(buf,"w");

  ds->lamfric=0.0;
  ds->lamdt=0.0;

  return 0;
}

int DataSpace_getN ( DataSpace * ds ) {
  if (ds) return ds->N;
  else return -1;
}

//int * DataSpace_chid ( DataSpace * ds ) {
//FIXME: This code was here to allow compute chapeau functions separatedly
//for different pair types of particles. For instance, this allow to
//recover SOD SOD, CLA CLA and SOD CLA pair potentials in 1 TAMD
//simulation. Each index has a number in ch_id which allow to sort the pair
//in the different chapeau objects on the c code.  From the studies with
//SOD CLA, this pair potentials will be OK only if the ficticius
//temperature is the same that the real one.  On the other hand, a better
//way to achive this is needed (without saving a lot of numbers in ch_id).
//For understand how this worked, see the previous versions of the code.
//  if (ds) {return ds->ch_id;}
//  return NULL;
//}
 
double * DataSpace_centerPos ( DataSpace * ds, int i ) {
  if (ds) {
    if (i<ds->N) {
      return ds->R[i];
    }
    return NULL;
  }
  return NULL;
}

int DataSpace_checkdata ( DataSpace * ds ) {
  if (ds) {
    int i,j;
    for (i=0;i<ds->N;i++) {
      for (j=0;j<3;j++) {
	if (ds->R[i][j]!=ds->R[i][j]) {
	  return 1;
	}
      }
    }
    return 0;
  }
  return 0;
}

int DataSpace_dump ( DataSpace * ds ) {
  if (ds) {
    int i,j;
    for (i=0;i<ds->N;i++) {
      fprintf(stdout,"%i ",i);
      for (j=0;j<3;j++)
	fprintf(stdout,"%.5f ",ds->R[i][j]);
      fprintf(stdout,"\n");
      fflush(stdout);
    }
  }
  return 0;
}

int DataSpace_AddAtomCenter ( DataSpace * ds, int n, int * ind, double * m ) {
  if (ds) {
    if (ds->iN<ds->N) {
      int i,j;
      ds->ac[ds->iN++]=New_atomCenterStruct(n);
      i=ds->iN-1;
      ds->ac[i]->ind=calloc(n,sizeof(int));
      ds->ac[i]->m=calloc(n,sizeof(double));
      for (j=0;j<n;j++) {
	ds->ac[i]->ind[j]=ind[j];
	ds->ac[i]->m[j]=m[j];
	ds->ac[i]->M+=m[j];
      }
      return (ds->iN-1);
    }
    return -1;
  }
  return -1;
}

int DataSpace_AddCV ( DataSpace * ds, char * typ, int nind, int * ind ) {
  int ityp=cv_getityp(typ);

  if (!ds) return -1;
  
  if (ds->iM<ds->M) {
    ds->cv[ds->iM++]=New_cvStruct(ityp,nind,ind);
    return (ds->iM-1);
  }
  return -1;
  
}

restrStruct * DataSpace_AddRestr ( DataSpace * ds, double k, double z, int nCV, double * cvc, char * rftypstr, 
			 double zmin, double zmax,char * boundf, double boundk ) {

  if (!ds) return NULL;
  
  if (ds->iK<ds->K) {
    ds->restr[ds->iK]=New_restrStruct(k,z,nCV,cvc,rftypstr,zmin,zmax,boundf,boundk);
    ds->iK++;
    return ds->restr[ds->iK-1];
  }

  return NULL;
}

int restr_set_rchid ( restrStruct * r, DataSpace * ds, int chid) {
  if (!r) return -1;
  if (!ds) return -1;
  if (chid>=ds->ch_num) return -1;
  r->chid = chid;
  return 0;
}
 
int restr_AddTamdOpt ( restrStruct * r, double g, double kt, double dt ) {
  
  if (!r) return -1;

  r->evolveFunc = cbd;
  r->tamdOpt=New_tamdOptStruct(g,kt,dt,r->rfityp);
  r->chid = 0;
  return 0;
}

int restr_UpdateTamdOpt ( restrStruct * r, double g, double kt, double dt ) {
  tamdOptStruct * tamd;
  if (!r) return -1;
  tamd = r->tamdOpt;
  tamd->kT = kt;
  tamd->gamma = g;
  tamd->dt = dt;
  tamd->ginv = 1.0/g;
  tamd->noise = sqrt(2.0*kt*dt*tamd->ginv);
  // rescalevels [expr sqrt(1.0*$NEWTEMP/$OLDTEMP)] no needed in cbd
  return 0;
}
 
int restr_AddSmdOpt  ( restrStruct * r, double target, int t0, int t1 ) {
  if (!r) return -1;
  r->evolveFunc = uniformvelocity;
  r->smdOpt = New_smdOptStruct(target,t0,t1);
  r->chid = 0;
  return 0;
}

double restr_getz ( restrStruct * r ) {
  if (!r) return -1;
  return r->z;
}

double restr_getu ( restrStruct * r ) {
  if (!r) return -1;
  return r->u;
}
  
int DataSpace_SetRestraints ( DataSpace * ds, double * rval ) {
  int i;
  restrStruct * r;
  for (i=0;i<ds->iK;i++) {
    r=ds->restr[i];
    r->z=rval[i];
    printf("CFACV/C) Reinitialize Z[%i] = %.5f\n",i,r->z);
  }
  return 0;
}

int DataSpace_ComputeCVs ( DataSpace * ds ) {
  int i;
  cvStruct * c;

 
  if (!ds) return -1;
  for (i=0;i<ds->iM;i++){
    c=ds->cv[i];
    c->calc(c,ds);
  }
  return 0;
}

int DataSpace_InitKnots ( DataSpace * ds, char * filename, int j) {
  if (ds) {
    chapeau * ch = ds->ch[j];
    FILE * fp = fopen(filename,"r");
    double * knots=(double*)calloc(ch->m,sizeof(double));
    int i=0;
    char ln[255];
    while (fgets(ln,255,fp)) {
      sscanf(ln,"%lf",&knots[i++]);
    }
    fclose(fp);
    chapeau_setPeaks(ch,knots);
    free(knots);
    fprintf(stdout,"INFO) Read %i knot values from %s\n",i,filename);
    fflush(stdout);
  }
  return 0;
}

int DataSpace_RestrainingForces ( DataSpace * ds, int first, int timestep ) {
  int ii,i,j,jj,k;
  double * cvc,fi;
  int d=0;
  int N=ds->N,K=ds->iK;
  int ic,aux;
  cvStruct * cv;
  restrStruct * r;
  chapeau * ch;
      
  // install SIGINT handler
  signal(SIGINT, handler);   
  if (!ds) return -1;

  /* clear out the position arrays to hold forces to be communicated
     back to tclForces */
  for (i=0;i<ds->N;i++) for (d=0;d<3;d++) ds->R[i][d]=0.0;

  /* Loop over all restraints and compute restraint values from their CV's */
  for (i=0;i<ds->iK;i++) {
    r=ds->restr[i];
    cvc=r->cvc;
    r->val=0.0;
    for (j=0;j<r->nCV;j++) {
      cv=ds->cv[j];
      r->val+=cvc[j]*cv->val;
#ifdef _PARANOIA_
      if (_PARANOIA_) {
        if (r->val!=r->val) {
          fprintf(stderr,"CFACV/C/PARANOIA) Tripped at r->val %.5lf\n",r->val);
          fprintf(stderr,"CFACV/C/PARANOIA) cvc[%i] %.5lf cv->val %.5lf\n",j,cv->val);
          fprintf(stderr,"Program exits.\n");
          fflush(stderr);
          exit(-1);
        }
      }
#endif
    }
    
    if (first) {
      if (r->tamdOpt) {
        r->z=r->val;
      }
      if (r->smdOpt) {
        r->z=r->val;
        smdOptInit(r->smdOpt,r->val);
      }
    }
  }

  //printf("CFACV/C) %i restraint %i val %.4f targ %.4f\n",timestep,i,r->val,r->z);

  /* compute forces on all restraints */
  for (i=0;i<ds->iK;i++) {
    r=ds->restr[i];
    cvc=r->cvc;

    // if this is a cartesian single-cv restraint, get the dimension
    for (j=0;j<r->nCV;j++) {
      if (cvc[j]) {
          d=cv_dimension(ds->cv[j]);
      }
    }

    // Compute the energy and the force of the restrain
    r->energyFunc(r);

    /* accumulate force increments in ds->R for each real particle*/
    /* cv->gr[jj][kk] = \partial CV / \partial (kk-coord of center jj of CV) */
    for (j=0;j<r->nCV;j++) {
      cv=ds->cv[j];
      /* if the j'th cv referenced by this restraint contributes to this restraint: */
      if (r->cvc[j]) {
        /* ... for each center in this CV ... */
        for (jj=0;jj<cv->nC;jj++) {
          /* increment the forces of the jj'th cartesian center on the j'th cv of this restraint */
          for (d=0;d<3;d++) ds->R[ cv->ind[jj] ][d]+=r->cvc[j]*cv->gr[jj][d]*r->f;
#ifdef _PARANOIA_
          if (_PARANOIA_) {
            if ((ds->R[cv->ind[jj]][0]!=ds->R[cv->ind[jj]][0])
      	  ||(ds->R[cv->ind[jj]][1]!=ds->R[cv->ind[jj]][1])
      	  ||(ds->R[cv->ind[jj]][2]!=ds->R[cv->ind[jj]][2])) {
      	fprintf(stderr,"CFACV/C/PARANOIA) Tripped.\n");
      	fprintf(stderr,"CFACV/C/PARANOIA) r->cvc[%i] %.5lf\n",j,r->cvc[j]);
      	fprintf(stderr,"CFACV/C/PARANOIA) cv->gr[%i][0] %.5lf\n",jj,cv->gr[jj][0]);
      	fprintf(stderr,"CFACV/C/PARANOIA) cv->gr[%i][1] %.5lf\n",jj,cv->gr[jj][1]);
      	fprintf(stderr,"CFACV/C/PARANOIA) cv->gr[%i][2] %.5lf\n",jj,cv->gr[jj][2]);
      	fprintf(stderr,"CFACV/C/PARANOIA) r->f %.5lf\n",r->f);
      	fprintf(stderr,"Program exits\n");
      	fflush(stderr);
      	exit(-1);
            }
          }
#endif
        }
      }
    }
  }


  // Evolution of the z variables
  // TODO: Im wondering if the loop over restrain can be saftely
  // place outside the if
  if (ds->doAnalyticalCalc) {

    // Add statistic to b and A matrix (tridiagonal)
    if (timestep>ds->beginaccum) {
      for (i=0;i<ds->iK;i++) {
        r=ds->restr[i];
        if (r->tamdOpt) fes1D(ds,r);
      }
    }

    // Moving the z variables
    if (!ds->useTAMDforces) {

      // Free energy gradient to evolve auxiliariy variables

      // map the interparticle gradients back into the restraints
      // compute gradients to be applied to z's from an analytical
      // function rather than from the atomic forces in this case,
      // each restraint must have one and only one CV

      fprintf(stderr,"The needed code was temporal removed. \n");
      fprintf(stderr,"Please useTAMDforces or see the original version of the code \n");
      exit(1);

      /*for (i=0;i<N;i++) {ds->F[i][0]=ds->F[i][1]=ds->F[i][2]=0.0;}
        for (i=0;i<K;i++) {
          r=ds->restr[i];
          // each restraint must have one and only one CV
          for (j=0;!r->cvc[j];j++);
          d=i%3;
          ii = ds->cv[j]->ind[0]; // index of particle "i"
          nl=ds->nl[ii];
          inlc=ds->nlc[ii];
          for (jj=0;jj<inlc;jj++) {
            k=nl[jj];
            if ( ds->rr[ii][k] < ds->squaredPairCutoff ) {
              fi=ds->gr[ii][k]*ds->RR[ii][k][d];
              ds->F[ii][d]-=fi;
              ds->F[k][d]+=fi;
            }
          }
          r->evolveFunc(r,ds->F[ii][d]);
        }*/
    } else {

      // TAMD forces to evolve auxilary variables
      for (i=0;i<K;i++) {
        r=ds->restr[i];

        // Compute the boundary forces not included in the atom forces!!
        r->boundFunc(r);

        if (r->tamdOpt) r->evolveFunc(r,-r->f);
      }
    }


    // Solving the chapeau equations
    if (ds->evolveAnalyticalParameters) {

      for (ic=0;ic<ds->ch_num;ic++) {
        ch=ds->ch[ic];

        if (timestep>ds->beginaccum && !(timestep % ch->nupdate)) {
          chapeau_update_peaks(ch);
          //aux=ch->nsample/ch->nupdate
          //if(!(aux%ch->outputFreq)) chapeau_output(ch);
          chapeau_output(ch,timestep);
        }
      }
    }

  } else { // update z's using the measured forces
    for (i=0;i<ds->iK;i++) {
      r=ds->restr[i];

      // Compute the boundary forces not included in the atom forces!!
      r->boundFunc(r);

      if (r->tamdOpt) r->evolveFunc(r,-r->f);
    }
  }

  //Evolve smd restrains (no matter the value of ds->doAnalyticalCalc)
  for (i=0;i<ds->iK;i++) {
    r=ds->restr[i];
  
    if (r->smdOpt) {
      
      // Compute the boundary forces not included in the atom forces!!
      r->boundFunc(r);

      r->evolve=(int)(r->smdOpt->t0<=timestep)&&(r->smdOpt->t1>=timestep);
      r->evolveFunc(r,r->smdOpt->increment);
    }
  }

  //Save restr trajectory for future restart
  if (!(timestep % ds->restrsavefreq)) ds_saverestrains(ds,timestep);

  return 0;
}

double DataSpace_RestraintEnergy ( DataSpace * ds ) {
  if (ds) {
    double u=0.0;
    int i;
    restrStruct * r;
    for (i=0;i<ds->iK;i++) {
      r=ds->restr[i];
      u+=r->u;
    }
    return u;
  }
  return 0.0;    
}

void DataSpace_ReportCV ( DataSpace * ds, int * active, double * res ) {
  if (ds) {
    int i;
    for (i=0;i<ds->iM;i++) {
      if (active[i]) res[i]=ds->cv[i]->val;
      else res[i]=0.0;
    }
  }
}

void DataSpace_ReportAll ( DataSpace * ds ) {
  if (ds) {
    int i,j,k;
    printf("CFACV_OTFP/C) DataSpace: %i Centers, %i CV's, %i Restraints\n",
	   ds->N, ds->iM, ds->iK);
    for (i=0;i<ds->N;i++) {
      printf("CFACV_OTFP/C)  Center %i : %.5f %.5f %.5f\n",
	     i,ds->R[i][0],ds->R[i][1],ds->R[i][2]);
    }
    for (i=0;i<ds->iM;i++) {
      printf("CFACV_OTFP/C)      CV %i : typ %s val %.5f ind ",i,cv_getstyp(ds->cv[i]->typ),
	     ds->cv[i]->val);
      for (j=0;j<ds->cv[i]->nC;j++) printf("%i ",ds->cv[i]->ind[j]);
      printf(" grads: ");
      for (j=0;j<ds->cv[i]->nC;j++) 
	for (k=0;k<3;k++) printf("%.5f ",ds->cv[i]->gr[j][k]);
      printf("\n");
    }
    for (i=0;i<ds->iK;i++) {
      printf("CFACV_OTFP/C)       R %i : k %.3f z %.3f cvc ",i,ds->restr[i]->k,ds->restr[i]->z);
      for (j=0;j<ds->restr[i]->nCV;j++) printf("%.2f ",ds->restr[i]->cvc[j]);
      if (ds->restr[i]->tamdOpt) {
	printf("TAMD(kT=%.2f,gamma=%.2f,dt=%.3f) ",
	       ds->restr[i]->tamdOpt->kT,
	       ds->restr[i]->tamdOpt->gamma,
	       ds->restr[i]->tamdOpt->dt);
      }
      printf(" typ %s min %.5f max %.5f\n",rf_getstyp(ds->restr[i]->rfityp),ds->restr[i]->min,ds->restr[i]->max);
      printf("\n");
    }
    
  }
}

void DataSpace_ReportRestraints ( DataSpace * ds, int step, int outputlevel, FILE * fp ) {
  int i;
  
  if (outputlevel & 1) {
    fprintf(fp,"CFACV_OTFP/C) Z  % 10i ",step);
    for (i=0;i<ds->iK;i++) {
      fprintf(fp,"% 11.5f",ds->restr[i]->z);
    }
    fprintf(fp,"\n");
  }
  if (outputlevel & 2) {
    fprintf(fp,"CFACV_OTFP/C) Th % 10i ",step);
    for (i=0;i<ds->iK;i++) {
      fprintf(fp,"% 11.5f",ds->restr[i]->val);
    }
    fprintf(fp,"\n");
  }
  if (outputlevel & 4) {
    fprintf(fp,"CFACV_OTFP/C) FD % 10i ",step);
    for (i=0;i<ds->iK;i++) {
      fprintf(fp,"% 11.5f",ds->restr[i]->tamd_restraint);
    }
    fprintf(fp,"\n");
  }
  if (outputlevel & 8) {
    fprintf(fp,"CFACV_OTFP/C) ND % 10i ",step);
    for (i=0;i<ds->iK;i++) {
      fprintf(fp,"% 11.5f",ds->restr[i]->tamd_noise);
    }
    fprintf(fp,"\n");
  }
  fflush(fp);
/*   printf("CFACV_OTFP/C) INFO %i ",step); */
/*   for (i=0;i<ds->iK;i++) { */
/*     printf(" | typ %s min %.5lf max %.5lf ",rf_getstyp(ds->restr[i]->rfityp),ds->restr[i]->min,ds->restr[i]->max); */
/*   } */
/*   printf("\n"); */
}

void DataSpace_BinaryReportRestraints ( DataSpace * ds, int step, int outputlevel, FILE * fp ) {
  int i;
  
  if (outputlevel & 1) {
    for (i=0;i<ds->iK;i++) {
      fwrite(&(ds->restr[i]->z),sizeof(double),1,fp);
    }
  }
  if (outputlevel & 2) {
    for (i=0;i<ds->iK;i++) {
      fwrite(&(ds->restr[i]->val),sizeof(double),1,fp);
    }
  }
  if (outputlevel & 4) {
    for (i=0;i<ds->iK;i++) {
      fwrite(&(ds->restr[i]->tamd_restraint),sizeof(double),1,fp);
    }
  }
  if (outputlevel & 8) {
    for (i=0;i<ds->iK;i++) {
      fwrite(&(ds->restr[i]->tamd_noise),sizeof(double),1,fp);
    }
  }
  fflush(fp);
}

int fes1D( DataSpace * ds, restrStruct * r ) { 
  int i,j,k,m;
  double R,F;
  chapeau * ch;

  R=r->z;
  F=-r->f; // F=k(\theta-z)

  // Following the paper E definition A and b should be
  //
  // A=1./nsteps \sum_{steps} \sum_{i} dphi_m(z)/dz dphi_n(z)/dz
  // b=-1./nsteps \sum_{steps} \sum_{i} dphi_m(z)/dz F
  //
  // But here we compute -b*nsteps and A*nsteps and then A and b are scaled
  // by the proper factor when equation A\lambda=b is solved (see the
  // chapeau_update_peaks procedure)
  
  // TODO: Repair the support to evolve trough Free Energy Gradient has was
  // done before
  
  /* TODO: This should be generalized by putting and index in the restrain to
           refer the asociated chapeu*/
  //ch=ds->ch[ ds->ch_id[m] ];
  ch=ds->ch[r->chid];

  if ( R > ch->rmax ) return;

  //FIXME: Early return if no chapeau object is asociated with this pair
  //if ( ds->ch_id[m] == -1 ) continue;

  // Early return to avoid interpolations below rmin
  if ( R <= ch->rmin ) return;

  //ch->nsample++;

  /* Interpolation R->m
        m R
        | |
    o---o-x-o---o--*/
  m=(int)((R-ch->rmin)*ch->idr);
  ch->hits[m]=1;
  //ch->hits[m+1]=1;

  /* 1/dz  /| 
          / |n=m+1
         /  | 
    o---o-x-o---o--*/ 
  ch->b->data[m+1]+=ch->idr*F;
  ch->A->data[(m+1)*ch->A->tda+(m+1)]+=ch->idr*ch->idr;

  /*    |\   -1/dz
     n=m| \
        |  \
    o---o-x-o---o--*/ 
  ch->b->data[m]-=ch->idr*F;
  ch->A->data[m*ch->A->tda+m]+=ch->idr*ch->idr;

  /*    |\ /| 
       m| X |n=m+1
        |/ \| 
    o---o-x-o---o--*/ 
  ch->A->data[(m+1)*ch->A->tda+m]-=ch->idr*ch->idr;
  ch->A->data[m*ch->A->tda+m+1]-=ch->idr*ch->idr;

            

  return 0;
}
 
