#include "chapeau.h"

chapeau * chapeau_alloc ( int m, double rmin, double rmax, int npart ) {
  chapeau * ch = (chapeau*)malloc(sizeof(chapeau));
  int d,i;

  ch->m=m;
  ch->N=npart;
  ch->rmin=rmin;
  ch->rmax=rmax;
  ch->dr=(rmax-rmin)/(m-1);
  ch->idr=1.0/ch->dr;
  ch->lam=gsl_vector_calloc(m);
  ch->lambar=gsl_vector_calloc(m);
  ch->newlam=gsl_vector_calloc(m);
  ch->hits=(int*)calloc(m,sizeof(int));

  ch->s=(double***) calloc(npart,sizeof(double**));
  for (i=0;i<npart;i++) {
    ch->s[i]=(double**)calloc(3,sizeof(double*));
    for (d=0;d<3;d++) {
      ch->s[i][d]=(double*)calloc(m,sizeof(double));
    }
  }
  
  ch->b=gsl_vector_calloc(m);
  ch->A=gsl_matrix_calloc(m,m);

  ch->bbar=gsl_vector_calloc(m);
  ch->Abar=gsl_matrix_calloc(m,m);

  ch->Permutation=gsl_permutation_alloc(m);

  fprintf(stderr,"CFACG4/DEBUG) allocated chapeau with %i knots on [%.5lf,%.5lf] dr %g \n",m,rmin,rmax,ch->dr);
  fflush(stderr);

  ch->alpha=0;

  return ch;
}

void chapeau_setUpdateInterval ( chapeau * ch, int i ) {
  if (ch) {
    ch->updateinterval=i;
  }
}

void chapeau_setPeaks ( chapeau * ch, double * peaks ) {
  if (ch) {
    if (peaks) {
      int i;
      for (i=0;i<ch->m;i++) gsl_vector_set(ch->lam,i,peaks[i]);
    }
  }
}

double right ( chapeau * ch, int m, double z, double zmin ) {
  return 1.0-ch->idr*(z-(m*ch->dr+zmin));
}

double left ( chapeau * ch, int m, double z, double zmin ) {
  return ch->idr*(z-((m-1)*ch->dr+zmin));
}

void chapeau_pair_eval_g ( chapeau * ch, double z, double * u, double * g_r ) {
  int m=(int)((z-ch->rmin)/ch->dr);
/*   fprintf(stdout,"### z %g z[%i] %g : interpolating between [%i](%g,%g)(%g) and [%i](%g,%g)(%g)\n", */
/* 	  z,m,gsl_vector_get(ch->lam,m),m,ch->rmin+m*ch->dr, */
/* 	  gsl_vector_get(ch->lam,m),gsl_vector_get(ch->lam,m)*right(ch,m,z,ch->rmin), */
/* 	  m+1,ch->rmin+(m+1)*ch->dr,gsl_vector_get(ch->lam,m+1),gsl_vector_get(ch->lam,m+1)*left(ch,m+1,z,ch->rmin)); */
  *u=gsl_vector_get(ch->lam,m+1)*left(ch,m+1,z,ch->rmin)+gsl_vector_get(ch->lam,m)*right(ch,m,z,ch->rmin);
  *g_r=(gsl_vector_get(ch->lam,m+1)-gsl_vector_get(ch->lam,m))/(z*ch->dr);
  ch->hits[m]++;
}

void chapeau_setupoutput ( chapeau * ch, char * filename, int outputFreq, int outputLevel ) {
  if (ch) {
    char * flag="CxCx";
    ch->ofp=fopen(filename,"w");
    fwrite(flag,sizeof(char),4,ch->ofp);
    fwrite(&outputLevel,sizeof(int),1,ch->ofp);
    fwrite(&ch->m,sizeof(int),1,ch->ofp);
    ch->outputFreq=outputFreq;
    ch->outputLevel=outputLevel;
  }
}

void chapeau_output ( chapeau * ch, int timestep ) {
  if (ch&&ch->ofp&&!(timestep%ch->outputFreq)) {
    int outputlevel=ch->outputLevel;
    int i;
    //fprintf(stdout,"## chapeau_output %i\n",timestep);
    fwrite(&timestep,sizeof(int),1,ch->ofp);
    if (outputlevel & 1) { // 0th bit = output knots as y
      for (i=0;i<ch->m;i++) {
	fwrite(gsl_vector_ptr(ch->lam,i),sizeof(double),1,ch->ofp);
	//fprintf(stderr,"### %i %g\n",i,*gsl_vector_ptr(ch->lam,i));
      }
    }
    if (outputlevel & 2) {
      for (i=0;i<ch->m;i++) {
	fwrite(&(ch->hits[i]),sizeof(int),1,ch->ofp);
      }
    }
    fflush(ch->ofp);
  }
}

void chapeau_init_global_accumulators ( chapeau * ch ) {
  int i,j,p,q;
  int m=ch->m;
  gsl_matrix_set_zero(ch->A);
  gsl_vector_set_zero(ch->b);
}

void chapeau_init_particle_sums ( chapeau * ch ) {
  int d,i,j;
  int m=ch->m;
  for (j=0;j<ch->N;j++) for (d=0;d<3;d++) for (i=0;i<m;i++) ch->s[j][d][i]=0.0;
}

void chapeau_increment_particle_sum ( chapeau * ch, int i, int j, double * Zij, double zij ) {
  int d;
  int m=(int)((zij-ch->rmin)*ch->idr);
  double ** si = ch->s[i];
  double ** sj = ch->s[j];
  double zijinv=1.0/(zij*ch->dr),zz;
  for (d=0;d<3;d++) {
    zz=zijinv*Zij[d];
    si[d][m]-=zz;
    si[d][m+1]+=zz;
    sj[d][m]+=zz;
    sj[d][m+1]-=zz;
  }
}

// increment the global accumulators by what is currently in the single-particle accumulators
void chapeau_increment_global_accumulators ( chapeau * ch, int i, double * F ) {
  int m,n;//,d;
  int N=ch->m;
  double ** s=ch->s[i];
  gsl_vector * b = ch->b;
  gsl_matrix * A = ch->A;
  //  for (d=0;d<3;d++) {
    for (m=0;m<N;m++) {
      b->data[m]+=F[0]*s[0][m]+F[1]*s[1][m]+F[2]*s[2][m];
      //gsl_vector_set(ch->b,m,gsl_vector_get(ch->b,m)+F[0]*s[0][m]+F[1]*s[1][m]+F[2]*s[2][m]);
      for (n=0;n<N;n++) {
	A->data[m*A->tda+n]+=s[0][m]*s[0][n]+s[1][m]*s[1][n]+s[2][m]*s[2][n];
	//gsl_matrix_set(ch->A,m,n,gsl_matrix_get(ch->A,m,n)+s[0][m]*s[0][n]+s[1][m]*s[1][n]+s[2][m]*s[2][n]);
	//fprintf(stderr,"*** iga %i %i %g\n",m,n,gsl_matrix_get(ch->A,m,n));
      }
    }
    //}
}

void chapeau_update_peaks ( chapeau * ch, int nsamples, int timestep ) {
  int i,j,I,ii,jj;
  int N=ch->m;
  int s;
  double ninv;
  double lo,lb,alpha;

  if (!(timestep%ch->updateinterval)) {
    
    gsl_matrix * Abar;
    gsl_vector * bbar;
    gsl_vector * lambar;
    gsl_permutation * p;
    int nred;

    // parameters that are allowed to evolve lie between indices for which ch->hits[] is non-zero
    // so extract the proper subspace
    for (i=0;i<ch->m&&(!ch->hits[i]);i++);
    I=i;
    nred=N-1-I;
    Abar=gsl_matrix_alloc(nred,nred);
    bbar=gsl_vector_alloc(nred);
    lambar=gsl_vector_alloc(nred);
    p=gsl_permutation_alloc(nred);
    ii=0;
    for (i=I;i<ch->m-1;i++) {
      gsl_vector_set(bbar,ii,gsl_vector_get(ch->b,i));
      jj=0;
      for (j=I;j<ch->m-1;j++) {
	gsl_matrix_set(Abar,ii,jj,gsl_matrix_get(ch->A,i,j));
	jj++;
      }
      ii++;
    }

    ninv=1.0/nsamples;
    gsl_matrix_scale(Abar,ninv);
    gsl_vector_scale(bbar,-ninv);
    
    gsl_linalg_LU_decomp(Abar,p,&s);
    gsl_linalg_LU_solve(Abar,p,bbar,lambar);
    
    alpha=0.0;//1.0-exp(-1.0e4*ninv);
    
    //fprintf(stderr,"### update %i subspace is %i x %i alpha %g\n",timestep,nred,nred,alpha);fflush(stderr);

    // update the vector of coefficients
    for (i=0;i<nred;i++) {
      lo=gsl_vector_get(ch->lam,i+I);
      lb=gsl_vector_get(lambar,i);
      gsl_vector_set(ch->lam,i+I,alpha*lo+(1-alpha)*lb);
    }
    gsl_matrix_free(Abar);
    gsl_vector_free(bbar);
    gsl_vector_free(lambar);
    gsl_permutation_free(p);
  }
}
