#include "cfacv.h"

//#define _PARANOIA_ 1

void cfacvBanner ( void ) {
  printf("CFACV/C) Version %.2f\n",CFACVC_VERSION);
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

char * CVSTRINGS[NULL_CV] = {"BOND", "ANGLE", "DIHED", "CARTESIAN_X", "CARTESIAN_Y", "CARTESIAN_Z"};
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

char * RFSTRINGS[NULL_RF] = {"HARMONIC", "PERIODIC"};
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

cvStruct * New_cvStruct ( int typ, int nC, int * ind ) {
  cvStruct * newcv=malloc(sizeof(cvStruct));
  int i;
  newcv->typ=typ;
  newcv->nC=nC;
  newcv->val=0.0;
  newcv->ind=calloc(nC,sizeof(int));
  if (ind) for (i=0;i<nC;i++) newcv->ind[i]=ind[i];
  newcv->gr=(double**)malloc(nC*sizeof(double*));
  for (i=0;i<nC;i++) newcv->gr[i]=(double*)malloc(3*sizeof(double));
  return newcv;
}

double rf_Harmonic ( double k, double v, double z, double half_domain ) {
  return -k*(v-z);
}

double rf_HarmonicCartPBC ( double k, double v, double z, double half_domain ) {
  double tmp=v-z;
  if (tmp>half_domain) {
    //    fprintf(stderr,"CFACV/C) DEBUG: -MICing %.5lf > %.5lf\n",tmp,half_domain);fflush(stderr);
    tmp-=2*half_domain;
  }
  if (tmp<-half_domain) {
    //    fprintf(stderr,"CFACV/C) DEBUG: +MICing %.5lf > %.5lf\n",tmp,half_domain);fflush(stderr);
    tmp+=2*half_domain;
  }
  return -k*tmp;
}

double re_HarmonicCartPBC ( double k, double v, double z, double half_domain ) {
  double tmp=v-z;
  if (tmp>half_domain) tmp-=2*half_domain;
  if (tmp<-half_domain) tmp+=2*half_domain;
  return 0.5*k*tmp*tmp;
}

double re_Harmonic ( double k, double v, double z, double half_domain ) {
  return 0.5*k*(v-z)*(v-z);
}

double rf_Periodic ( double k, double v, double z, double half_domain ) {
  double dz=v-z;
#ifdef _PARANOIA_
  if (_PARANOIA_) {
    if ((dz!=dz)||(k!=k)||(half_domain!=half_domain)) {
      fprintf(stderr,"CFACV/C/PARANOIA) Tripped in rf_Periodic\n");
      fprintf(stderr,"CFACV/C/PARANOIA) k %.5lf v %.5lf z %.5lf dz %.5lf half_domain %.5lf\n",
	      k,v,z,dz,half_domain);
      fprintf(stderr,"Program exits\n");
      fflush(stderr);
      exit(-1);
    }
  }
#endif
  if (dz<-half_domain) dz+=2*half_domain;
  else if (dz>half_domain) dz-=2*half_domain;
  return -k*dz;
}

double re_Periodic ( double k, double v, double z, double half_domain ) {
  double dz=v-z;
  if (dz<-half_domain) dz+=2*half_domain;
  else if (dz>half_domain) dz-=2*half_domain;
  return 0.5*k*dz*dz;
}

restrStruct * New_restrStruct ( double k, double z, int nCV, double * cvc, char * rftypstr, double zmin, double zmax ) {
  restrStruct * newr=malloc(sizeof(restrStruct));
  int i;
  newr->rfityp=rf_getityp(rftypstr);
  newr->k=k;
  newr->z=z;
  newr->u=newr->f=0.0;
  newr->nCV=nCV;
  newr->cvc=(double*)malloc(nCV*sizeof(double));
  if (cvc) for (i=0;i<nCV;i++) newr->cvc[i]=cvc[i];
  else for (i=0;i<nCV;i++) newr->cvc[i]=0;
  newr->tamdOpt=NULL;
  newr->tamd_noise=0.0;
  newr->tamd_restraint=0.0;
  newr->smdOpt=NULL;
  newr->min=zmin;
  newr->max=zmax;
  newr->half_domain=0.5*(zmax-zmin);
  if (newr->rfityp==HARMONIC) {
    newr->forceFunc = rf_HarmonicCartPBC;
    newr->energyFunc = re_HarmonicCartPBC;
  }
  else if (newr->rfityp==PERIODIC) {
    newr->forceFunc = rf_Periodic;
    newr->energyFunc = re_Periodic;
  }
  
  return newr;
}

tamdOptStruct * New_tamdOptStruct ( double g, double kt, double dt, int riftyp, double half_domain ) {
  tamdOptStruct * tamd=malloc(sizeof(tamdOptStruct));
  tamd->kT=kt;
  tamd->gamma=g;
  tamd->dt=dt;
  if (g>0.0) tamd->ginv=1.0/g;
  else tamd->ginv=0.0;
  tamd->noise = sqrt(2.0*kt*dt*tamd->ginv);
  tamd->periodic=(riftyp==PERIODIC);
  tamd->half_domain=half_domain;
  return tamd;
}

/* Diffusive dynamics! */
double tamdUpdate ( DataSpace * ds, double z, double f, tamdOptStruct * tamd, double * noise, double * res ) {
  double dd = tamd->ginv*tamd->dt*f;
  double rd = tamd->noise*my_whitenoise(ds->Xi);
  *noise=rd;
  *res=dd;
  if (tamd->periodic) {
    double newrawz=z+dd+rd;
    if (newrawz>tamd->half_domain) newrawz-=2*tamd->half_domain;
    if (newrawz<-tamd->half_domain) newrawz+=2*tamd->half_domain;
    return newrawz;
  }
  return z+dd+rd; 
}

double smdUpdate ( double z, double increment, int OK ) {
  return z+OK*increment;
}

double smdUpdate_Periodic ( double z, double increment, int OK ) {
  double newz=z+increment;
  if (newz>M_PI) newz-=2*OK*M_PI;
  else if (newz<-M_PI) newz+=2*OK*M_PI;
  return newz;
}

smdOptStruct * New_smdOptStruct ( double target, int t0, int t1, int periodic ) {
  smdOptStruct * smd=malloc(sizeof(smdOptStruct));
  smd->t0=t0;
  smd->t1=t1;
  smd->invinterval = 1.0/(t1-t0);
  smd->target=target;
  if (periodic) smd->update=smdUpdate_Periodic;
  else smd->update=smdUpdate;
  return smd;
}

int smdOptInit ( smdOptStruct * smd, double initval, int periodic ) {
  double dist=0.0;
  if (smd) {
    smd->initval=initval;
    dist = smd->target-smd->initval;
    if (periodic) {
      if (dist<-M_PI) {
	dist+=2*M_PI;
      }
      else if (dist > M_PI) {
	dist-=2*M_PI;
      } 
    }
    smd->increment=dist*smd->invinterval;
    printf("CFACV/C) smd increment is %.5lf based on target %.5lf and initval %.5lf\n",
	   smd->increment,smd->target,smd->initval);
    return 0;
  }
  return -1;
}

DataSpace * NewDataSpace ( int N, int M, int K, long int seed ) {
  DataSpace * ds = malloc(sizeof(DataSpace));
  int i;
  ds->N=N;
  ds->M=M;
  ds->K=K;
  ds->iN=ds->iM=ds->iK=0;
  ds->Xi=(unsigned short *)malloc(3*sizeof(unsigned short));
  ds->Xi[0]=(seed & 0xffff00000000) >> 32;
  ds->Xi[1]=(seed & 0x0000ffff0000) >> 16;
  ds->Xi[2]= 0x330e;

  ds->R=(double**)malloc(N*sizeof(double*));
  for (i=0;i<N;i++) ds->R[i]=calloc(3,sizeof(double));

  ds->ac = (atomCenterStruct**)malloc(N*sizeof(atomCenterStruct*));

  ds->cv=(cvStruct**)malloc(M*sizeof(cvStruct*));
  ds->restr=(restrStruct**)malloc(K*sizeof(restrStruct*));

  /* metric tensor */
  ds->mt=(double**)malloc(M*sizeof(double*));
  for (i=0;i<M;i++) ds->mt[i]=(double*)malloc(M*sizeof(double));

  /* inter-center distances */
  ds->doAnalyticalCalc=0;
  ds->F=NULL;
  ds->Z=NULL;
  ds->RR=NULL;
  ds->rr=NULL;
  ds->gr=NULL;
  ds->O[0]=ds->O[1]=ds->O[2]=0.0;
  ds->L[0]=ds->L[1]=ds->L[2]=0.0;
  ds->hL[0]=ds->hL[1]=ds->hL[2]=0.0;
  ds->evolveAnalyticalParameters=0;
  ds->beginEvolveParameters=0;
  ds->reportParamFreq=0;

  return ds;
}

void LJPair ( double * par, double r2, double * gr ) {
  double ep=par[0];
  double sig=par[1];
  if (r2>0) {
    double r6i=1.0/(r2*r2*r2);
    double s6=sig*sig*sig*sig*sig*sig;
    double s12=s6*s6;
    double r12i=r6i*r6i;
    *gr=-12*ep/r2*((s12*r12i)-(s6*r6i));
/*     fprintf(stderr,"CFACV/C) DEBUG: LJPair: ep %.5lf sig %.5lf r2 %.5lf r6i %.5lf r12i %.5lf s6 %.5lf s12 %.5lf gr %.5lf\n", */
/* 	    ep,sig,r2,r6i,r12i,s6,s12,*gr);fflush(stderr); */
  } else {
    *gr=0.0;
  }
}

int DataSpace_SetupPairCalc ( DataSpace * ds, double Ox, double Oy, double Oz, \
			      double Lx, double Ly, double Lz, double cutoff, double nlcutoff,
			      int beginEvolve, int useTAMDforces, int reportParamFreq, double spline_min, int nKnots, 
			      char * splineoutputfile, int splineoutputfreq, int splineoutputlevel, 
			      int lamupdateinterval ) {
  int i,j;
  int N=ds->N; // number of centers
  fprintf(stderr,"CFACV/C) DEBUG DataSpace_SetupPairCalc N %i\n",N);fflush(stderr);
  ds->doAnalyticalCalc=1;
  ds->useTAMDforces=useTAMDforces;
  if (ds->useTAMDforces) fprintf(stdout,"CFACV/C) INFO: using TAMD forces instead of analytical forces.\n");
  ds->squaredPairCutoff = cutoff*cutoff;
  ds->nl_squaredPairCutoff = nlcutoff*nlcutoff;
  ds->nlDispThresh = 0.5*(nlcutoff-cutoff);
  fprintf(stderr,"CFACV/C) DEBUG DataSpace_SetupPairCalc nlDispThresh %g\n",ds->nlDispThresh);fflush(stderr);
  ds->evolveAnalyticalParameters=1;
  ds->beginEvolveParameters=beginEvolve;
  ds->reportParamFreq=reportParamFreq;

  ds->Z=(double**)malloc(N*sizeof(double*));
  for (i=0;i<N;i++) ds->Z[i]=calloc(3,sizeof(double));
  ds->accumZ=(double**)malloc(N*sizeof(double*));
  for (i=0;i<N;i++) ds->accumZ[i]=calloc(3,sizeof(double));
  ds->nlc=(int*)malloc(N*sizeof(int));
  ds->nl=(int**)malloc(N*sizeof(int*));
  for (i=0;i<N;i++) ds->nl[i]=(int*)malloc(MAXNBOR*sizeof(int));
  ds->nlTrig=0;

  ds->F=(double**)malloc(N*sizeof(double*));
  for (i=0;i<N;i++) ds->F[i]=calloc(3,sizeof(double));

  ds->RR=(double***)malloc(N*sizeof(double**));
  ds->rr=(double**)malloc(N*sizeof(double*));
  ds->gr=(double**)malloc(N*sizeof(double*));
  for (i=0;i<N;i++) {
    ds->RR[i]=(double**)malloc(N*sizeof(double*));
    ds->rr[i]=(double*)malloc(N*sizeof(double));
    ds->gr[i]=(double*)malloc(N*sizeof(double));
    for (j=0;j<N;j++) ds->RR[i][j]=(double*)malloc(3*sizeof(double));
  }
  ds->O[0]=Ox;
  ds->O[1]=Oy;
  ds->O[2]=Oz;
  ds->L[0]=Lx;  ds->hL[0]=Lx/2.;
  ds->L[1]=Ly;  ds->hL[1]=Ly/2.;
  ds->L[2]=Lz;  ds->hL[2]=Lz/2.;
  ds->Max[0]=ds->L[0]-ds->O[0]; ds->Min[0]=ds->Max[0]-ds->L[0];
  ds->Max[1]=ds->L[1]-ds->O[1]; ds->Min[1]=ds->Max[1]-ds->L[1];
  ds->Max[2]=ds->L[2]-ds->O[2]; ds->Min[2]=ds->Max[2]-ds->L[2];
 
  ds->nsamples=0;

  ds->ch=chapeau_alloc(nKnots,spline_min,cutoff,N);

  chapeau_init_global_accumulators(ds->ch);
  chapeau_setupoutput(ds->ch,splineoutputfile,splineoutputfreq,splineoutputlevel);
  chapeau_setUpdateInterval(ds->ch,lamupdateinterval);

  ds->lamfric=0.0;
  ds->lamdt=0.0;
  fprintf(stderr,"CFACV_OTFP/DEBUG) setup pair calculation\n");
  fflush(stderr);


  return 0;
}

int DataSpace_getN ( DataSpace * ds ) {
  if (ds) return ds->N;
  else return -1;
}

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
      fprintf(stderr,"%i ",i);
      for (j=0;j<3;j++)
	fprintf(stderr,"%.5lf ",ds->R[i][j]);
      fprintf(stderr,"\n");
      fflush(stderr);
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
  if (ds) {
    int ityp=cv_getityp(typ);
    if (ds->iM<ds->M) {
      ds->cv[ds->iM++]=New_cvStruct(ityp,nind,ind);
      return (ds->iM-1);
    }
    return -1;
  }
  return -1;
}

int DataSpace_AddRestr ( DataSpace * ds, double k, double z, int nCV, double * cvc, char * rftypstr, 
			 double zmin, double zmax ) {
  if (ds) {
    if (ds->iK<ds->K) {
      ds->restr[ds->iK++]=New_restrStruct(k,z,nCV,cvc,rftypstr,zmin,zmax);
      return (ds->iK-1);
    }
    return -1;
  }
  return -1;
}

int DataSpace_AddTamdOpt ( DataSpace * ds, int ir, double g, double kt, double dt ) {
  if (ds) {
/*     printf("DB: adding tamd options (% 7.3f% 7.3f% 7.3f) to r %i out of %i...\n", */
/* 	   g,kt,dt,ir,ds->iK); */
    if (ir < ds->iK ) {
      ds->restr[ir]->tamdOpt=New_tamdOptStruct(g,kt,dt,ds->restr[ir]->rfityp,ds->restr[ir]->half_domain);
      return ir;
    }
    return -1;
  }
  return -1;
}

int DataSpace_AddSmdOpt  ( DataSpace * ds, int ir, double target, int t0, int t1 ) {
  if (ds) {
    if (ir < ds->iK ) {
      ds->restr[ir]->smdOpt=New_smdOptStruct(target,t0,t1,(int)(ds->restr[ir]->forceFunc==rf_Periodic));
      return ir;
    }
    return -1;
  }
  return -1;
}

int DataSpace_SetRestraints ( DataSpace * ds, double * rval ) {
  int i;
  restrStruct * r;
  for (i=0;i<ds->iK;i++) {
    r=ds->restr[i];
    r->z=rval[i];
    printf("CFACV/C) Reinitialize Z[%i] = %.5lf\n",i,r->z);
  }
  return 0;
}

int DataSpace_ComputeCVs ( DataSpace * ds ) {
  if (ds) {
    int i,j,k,dimension,ii;
    double tmp;
    cvStruct * cvi, * cvj, * cv;
    restrStruct * r;
    /* clear out the metric tensor */
    for (i=0;i<ds->M;i++) {
      for (j=0;j<ds->M;j++) {
	ds->mt[i][j]=0.0;
      }
    }

    for (i=0;i<ds->iM;i++) {
      cvi=ds->cv[i];
      if (cvi->typ==BOND) {
	cvi->val=my_getbond(ds->R[cvi->ind[0]],ds->R[cvi->ind[1]],cvi->gr[0],cvi->gr[1]);
      } 
      else if (cvi->typ==ANGLE) {
	cvi->val=my_getangle(ds->R[cvi->ind[0]],ds->R[cvi->ind[1]],ds->R[cvi->ind[2]],
			    cvi->gr[0],        cvi->gr[1],        cvi->gr[2]);
      }
      else if (cvi->typ==DIHED) {
	cvi->val=my_getdihed(ds->R[cvi->ind[0]],ds->R[cvi->ind[1]],ds->R[cvi->ind[2]],ds->R[cvi->ind[3]],
			    cvi->gr[0],        cvi->gr[1],        cvi->gr[2],        cvi->gr[3]);
#ifdef _PARANOIA_
	if (_PARANOIA_) {
	  if (cvi->val!=cvi->val) {
	    fprintf(stderr,"CFACV/C/PARANOIA) Tripped at dihed cvi->val %.5lf\n",cvi->val);
	    fprintf(stderr,"Program exits.\n");
	    fflush(stderr);
	    exit(-1);
	  }
	}
#endif
      }
      else if (cvi->typ==CARTESIAN_X) {
	cvi->val=ds->R[cvi->ind[0]][0];
	cvi->gr[0][0]=1.0;
	cvi->gr[0][1]=0.0;
	cvi->gr[0][2]=0.0;
      }
      else if (cvi->typ==CARTESIAN_Y) {
	cvi->val=ds->R[cvi->ind[0]][1];
	cvi->gr[0][0]=0.0;
	cvi->gr[0][1]=1.0;
	cvi->gr[0][2]=0.0;
      }
      else if (cvi->typ==CARTESIAN_Z) {
	cvi->val=ds->R[cvi->ind[0]][2];
	cvi->gr[0][0]=0.0;
	cvi->gr[0][1]=0.0;
	cvi->gr[0][2]=1.0;
      }
    }
    /* metric tensor: innermost loop is loop over all centers; problem is
       two aribitrary cv's have different number of centers. */
    for (i=0;i<ds->iM;i++) {
      for (j=0;j<ds->iM;j++) {
        
      }
    }

    return 0;
  }
  return -1;
}


int DataSpace_InitKnots ( DataSpace * ds, char * filename ) {
  if (ds) {
    chapeau * ch = ds->ch;
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
}

double handle_pair ( DataSpace * ds, int i, int j ) {
  int k;
  double tmp,t2,R,u,gr_r,r2;
  chapeau * ch=ds->ch;
  r2=0.0;
  for (k=0;k<3;k++) {
    tmp=ds->Z[i][k]-ds->Z[j][k];
    // minimum image convention
    while (tmp>ds->hL[k]) tmp-=ds->L[k];
    while (tmp<-ds->hL[k]) tmp+=ds->L[k];
    ds->RR[i][j][k]=tmp;
    ds->RR[j][i][k]=-tmp;
    t2=tmp*tmp;
    r2+=t2;
  }
  ds->gr[i][j]=ds->gr[j][i]=0.0;
  if ( r2 < ds->squaredPairCutoff ) {
    R=sqrt(r2);
    chapeau_pair_eval_g(ch,R,&u,&gr_r);
    ds->gr[i][j]=gr_r;
    ds->gr[j][i]=gr_r;
    chapeau_increment_particle_sum(ch,i,j,&ds->RR[i][j][0],R);
  }
  ds->rr[i][j]=r2;
  ds->rr[j][i]=r2;
  return r2;
}

int DataSpace_RestrainingForces ( DataSpace * ds, int first, int timestep ) {
  if (ds) {
    int i,j,jj,k;
    double v, u, f, kr, z, dz;
    cvStruct * cv;
    restrStruct * r;
    double * cvc;
    int dimension;

    //   printf("CFACV/C) RestrainingForces begins.\n");

    /* clear out the position arrays to hold forces to be communicated
       back to tclForces */
    for (i=0;i<ds->N;i++) {
      ds->R[i][0]=0.0;
      ds->R[i][1]=0.0;
      ds->R[i][2]=0.0;
    }

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
	  smdOptInit(r->smdOpt,r->val,(int)(r->forceFunc==rf_Periodic));	  
	}
      }
    }
      //printf("CFACV/C) %i restraint %i val %.4lf targ %.4lf\n",timestep,i,r->val,r->z);

    /* compute forces on all restraints */
    for (i=0;i<ds->iK;i++) {
      r=ds->restr[i];
      cvc=r->cvc;
      // if this is a cartesian single-cv restraint, get the dimension
      for (j=0;j<r->nCV;j++) {
	if (cvc[j]) {
	  cv=ds->cv[j];
	  switch(cv->typ) {
	  case CARTESIAN_X: 
	    dimension=0;
	    break;
	  case CARTESIAN_Y:
	    dimension=1;
	    break;
	  case CARTESIAN_Z:
	    dimension=2;
	    break;
	  default:
	    fprintf(stderr,"CFACV/C) ERROR: can't handle cv typ %s in analytics (yet)\n",CVSTRINGS[cv->typ]);
	    exit(-1);
	  }
	}
      }

      r->f=r->forceFunc(r->k,r->val,r->z,ds->hL[dimension]);
      
      r->u=r->energyFunc(r->k,r->val,r->z,ds->hL[dimension]);

      //printf("CFACV/C) restraint %i: f %.5lf u %.5lf\n",i,r->f,r->u);

      /* accumulate force increments */
      /* cv->gr[jj][kk] = \partial CV / \partial (kk-coord of center jj of CV) */
      for (j=0;j<r->nCV;j++) {
	cv=ds->cv[j];
	/* if the j'th cv referenced by this restraint contributes to this restraint: */
	if (r->cvc[j]) {
	  /* ... for each center in this CV ... */
	  for (jj=0;jj<cv->nC;jj++) {
	    /* increment the forces of the jj'th cartesian center on the j'th cv of this restraint */
	    ds->R[cv->ind[jj]][0]+=r->cvc[j]*cv->gr[jj][0]*r->f;
	    ds->R[cv->ind[jj]][1]+=r->cvc[j]*cv->gr[jj][1]*r->f;
	    ds->R[cv->ind[jj]][2]+=r->cvc[j]*cv->gr[jj][2]*r->f;
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

    if (ds->doAnalyticalCalc) {
      int nCV,ii,d,N=ds->N,K=ds->iK;
      double this_force;
      int dimension;
      double r2,t2,R,gr,gr_r;
      //double total_force[3];
      double tmp,oldz,*ff,fi;
      chapeau * ch=ds->ch;
      int * nl, * nlc, inlc;

      ds->nlTrig=0;

      // copy r->z from each restraint/CV into ds->Z[][] and 
      // r->f into ds->F[][]...
      // restraint-index into particle-dimension index...
      for (i=0;i<K;i++) {
	r=ds->restr[i];     
	for (j=0;!r->cvc[j];j++);
	cv=ds->cv[j];
	ii=cv->ind[0]; // ii = particle-index 
	dimension = i%3;  // 0 == CARTESIAN_X, 1 = CARTESIAN_Y, 2 = CARTESIAN_Z
	//fprintf(stderr,"#### ts %i ptcl %i d %i z %g f %g\n",timestep,ii,dimension,r->z,-r->f);fflush(stderr);
	oldz=ds->Z[ii][dimension];
	ds->accumZ[ii][dimension]+=r->z-oldz;
	if (!first && fabs(ds->accumZ[ii][dimension])>ds->nlDispThresh) ds->nlTrig=1;
	ds->Z[ii][dimension]=r->z;
	ds->F[ii][dimension]=-r->f; // k(theta-z)
      }
     
      // begin the N2 loop that fills accumulators and computes analytical gradient
      chapeau_init_particle_sums(ch);
      if (first || ds->nlTrig) {
	fprintf(stderr,"CFACV_OTFP/INFO) resetting neighborlists @ %i\n",timestep);
	for (i=0;i<N;i++) {
	  ds->accumZ[i][0]=ds->accumZ[i][1]=ds->accumZ[i][2]=0.0;
	  nl=ds->nl[i];
	  nlc=ds->nlc;
	  nlc[i]=0;
	  for (j=i+1;j<N;j++) {
	    r2=handle_pair(ds,i,j);
	    // add to i's neighborlist
	    if ( r2 < ds->nl_squaredPairCutoff ) {
	      nl[nlc[i]++]=j;
	      if (nlc[i]==MAXNBOR) {
		fprintf(stderr,"# ERROR: particle %i has too many neighbors %i\n",i, MAXNBOR);
		exit(-1);
	      }
	    }
	  }
	  // FINALIZE this time-step's particle-i tallies
	  chapeau_increment_global_accumulators(ch,i,&ds->F[i][0]);
	  //	  for (d=0;d<3;d++) ds->F[i][d]=0.0;
	}
      } else {
	for (i=0;i<N;i++) {  
	  nl=ds->nl[i];
	  inlc=ds->nlc[i];
	  for (jj=0;jj<inlc;jj++) {
	    j=nl[jj];
	    r2=handle_pair(ds,i,j);
	  }
	  // FINALIZE this time-step's particle-i tallies
	  chapeau_increment_global_accumulators(ch,i,&ds->F[i][0]);
	  //	  for (d=0;d<3;d++) ds->F[i][d]=0.0;
	}
      } // end of pair loop with neighborlist reconstruction
      //fprintf(stderr,"CFACV/C) DEBUG paircalc
      //ends.\n");fflush(stderr); fflush(stderr); 

      // map the interparticle gradients back into the restraints
      // compute gradients to be applied to z's from an analytical
      // function rather than from the atomic forces in this case,
      // each restraint must have one and only one CV
      for (i=0;i<N;i++) {
	ds->F[i][0]=ds->F[i][1]=ds->F[i][2]=0.0;
      }

      for (i=0;i<K;i++) {
	r=ds->restr[i];
	for (j=0;!r->cvc[j];j++);
	dimension=i%3;
	ii = ds->cv[j]->ind[0]; // index of particle "i"
	ff=ds->F[ii];
	//fprintf(stderr,"###### ts %i ptcl %i d %i z %g f %g\n",timestep,ii,dimension,r->z,-r->f);fflush(stderr);
	nl=ds->nl[ii];
	inlc=ds->nlc[ii];
	for (jj=0;jj<inlc;jj++) {
	  k=nl[jj];
	  r2=ds->rr[ii][k];
	  if ( r2 < ds->squaredPairCutoff ) {
	    fi=ds->gr[ii][k]*ds->RR[ii][k][dimension];
	    ff[dimension]-= fi;
	    ds->F[k][dimension]+=fi;
	  }
	}
	r->z=tamdUpdate(ds,r->z,ds->useTAMDforces?-r->f:ff[dimension],r->tamdOpt,&r->tamd_noise,&r->tamd_restraint);	
      }
      
      if (ds->evolveAnalyticalParameters) ds->nsamples++;
      if (ds->evolveAnalyticalParameters&&(timestep>ds->beginEvolveParameters)&&(!(timestep%ch->updateinterval))) {
	chapeau_update_peaks(ch,ds->nsamples,timestep);
	chapeau_output(ch,timestep);
	if (ds->reportParamFreq && !(timestep%ds->reportParamFreq)) {
	  fprintf(stdout,"LAMBDA: % 8i",timestep);fflush(stderr);
	  for (i=0;i<ds->ch->m;i++) {
	    fprintf(stdout,"% 17.5lf (%i) ",gsl_vector_get(ds->ch->lam,i),ds->ch->hits[i]);
	  }
	  fprintf(stdout,"\n");
	}
      }
    } // doAnalyticalCalc
    else { // update z's using the measured forces
      for (i=0;i<ds->iK;i++) {
	r=ds->restr[i];
	if (r->tamdOpt) {
	  r->z=tamdUpdate(ds,r->z,-r->f,r->tamdOpt,&r->tamd_noise,&r->tamd_restraint);
	} else if (r->smdOpt) {
	  r->z=r->smdOpt->update(r->z,r->smdOpt->increment,
				 (r->smdOpt->t0<=timestep)&&(r->smdOpt->t1>=timestep));
	}
      }
    }
    return 0;
  }
  return -1;
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
      printf("CFACV_OTFP/C)  Center %i : %.5lf %.5lf %.5lf\n",
	     i,ds->R[i][0],ds->R[i][1],ds->R[i][2]);
    }
    for (i=0;i<ds->iM;i++) {
      printf("CFACV_OTFP/C)      CV %i : typ %s val %.5lf ind ",i,cv_getstyp(ds->cv[i]->typ),
	     ds->cv[i]->val);
      for (j=0;j<ds->cv[i]->nC;j++) printf("%i ",ds->cv[i]->ind[j]);
      printf(" grads: ");
      for (j=0;j<ds->cv[i]->nC;j++) 
	for (k=0;k<3;k++) printf("%.5lf ",ds->cv[i]->gr[j][k]);
      printf("\n");
    }
    for (i=0;i<ds->iK;i++) {
      printf("CFACV_OTFP/C)       R %i : k %.3lf z %.3lf cvc ",i,ds->restr[i]->k,ds->restr[i]->z);
      for (j=0;j<ds->restr[i]->nCV;j++) printf("%.2lf ",ds->restr[i]->cvc[j]);
      if (ds->restr[i]->tamdOpt) {
	printf("TAMD(kT=%.2f,gamma=%.2f,dt=%.3f) ",
	       ds->restr[i]->tamdOpt->kT,
	       ds->restr[i]->tamdOpt->gamma,
	       ds->restr[i]->tamdOpt->dt);
      }
      printf(" typ %s min %.5lf max %.5lf\n",rf_getstyp(ds->restr[i]->rfityp),ds->restr[i]->min,ds->restr[i]->max);
      printf("\n");
    }
    
  }
}

void DataSpace_ReportRestraints ( DataSpace * ds, int step, int outputlevel, FILE * fp ) {
  int i;
  
  if (outputlevel & 1) {
    fprintf(fp,"CFACV_OTFP/C) Z  % 10i ",step);
    for (i=0;i<ds->iK;i++) {
      fprintf(fp,"% 11.5lf",ds->restr[i]->z);
    }
    fprintf(fp,"\n");
  }
  if (outputlevel & 2) {
    fprintf(fp,"CFACV_OTFP/C) Th % 10i ",step);
    for (i=0;i<ds->iK;i++) {
      fprintf(fp,"% 11.5lf",ds->restr[i]->val);
    }
    fprintf(fp,"\n");
  }
  if (outputlevel & 4) {
    fprintf(fp,"CFACV_OTFP/C) FD % 10i ",step);
    for (i=0;i<ds->iK;i++) {
      fprintf(fp,"% 11.5lf",ds->restr[i]->tamd_restraint);
    }
    fprintf(fp,"\n");
  }
  if (outputlevel & 8) {
    fprintf(fp,"CFACV_OTFP/C) ND % 10i ",step);
    for (i=0;i<ds->iK;i++) {
      fprintf(fp,"% 11.5lf",ds->restr[i]->tamd_noise);
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

