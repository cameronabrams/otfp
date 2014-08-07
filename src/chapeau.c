#include "chapeau.h"

void chapeau_sum ( chapeau * ch1, chapeau * ch2 ) {
  //Overwrite ch1 with ch1+ch2
  int i,j,k;

  if (sizeof(ch1)!=sizeof(ch1)) {
    fprintf(stderr,"CFACV/C) ERROR: you can not sum chapeau objects with different sizes\n");
    exit(-1);
  }
  if (ch1->rmin!=ch2->rmin || ch1->rmax!=ch2->rmax || ch1->dr!=ch2->dr ) {
    fprintf(stderr,"CFACV/C) ERROR: you can not sum chapeau objects with different sizes\n");
    exit(-1);
  }

  for (i=0;i<ch1->N;i++) {
    if ( ch1->mask[i]!=ch2->mask[i] ) {
      fprintf(stderr,"CFACV/C) ERROR: you can not sum chapeau objects with different masks\n");
      exit(-1);
    }
  }    

  gsl_vector_add(ch1->b   ,ch2->b);
  gsl_matrix_add(ch1->A   ,ch2->A);

  for (i=0;i<ch1->m;i++) {
    ch1->hits[i] = ch1->hits[i] + ch2->hits[i];
  }

  // I think that the forward is irrelevant since I am saving the state before
  // the chapeau_output and therefore before the chapeau_update, but I let this
  // here just for the case.

  for (i=0;i<ch1->N;i++) {
    for (j=0;j<3;j++) {
      for (k=0;k<ch1->m;k++) {
        ch1->s[i][j][k] = ch1->s[i][j][k] + ch2->s[i][j][k];
      }
    }
  }
  gsl_vector_add(ch1->lam ,ch2->lam); //also irrelevant

}

chapeau * chapeau_alloc ( int m, double rmin, double rmax, int npart ) {
  chapeau * ch = (chapeau*)malloc(sizeof(chapeau));
  int d,i;

  if (npart<1) {
    fprintf(stderr,"CFACV/C) ERROR: chapeau objects without particles\n");
    exit(-1);
  }         

  if (m<1) {
    fprintf(stderr,"CFACV/C) ERROR: chapeau objects without bins\n");
    exit(-1);
  }         

  ch->m=m;
  ch->N=npart;
  ch->rmin=rmin;
  ch->rmax=rmax;
  ch->dr=(rmax-rmin)/(m-1);
  ch->idr=1.0/ch->dr;
  ch->lam=gsl_vector_calloc(m);
  ch->hits=(long*)calloc(m,sizeof(long));
  
  ch->s=(double***) calloc(npart,sizeof(double**));
  for (i=0;i<npart;i++) {
    ch->s[i]=(double**)calloc(3,sizeof(double*));
    for (d=0;d<3;d++) {
      ch->s[i][d]=(double*)calloc(m,sizeof(double));
    }
  }

  ch->mask=(int*)calloc(npart,sizeof(int));
  
  ch->b=gsl_vector_calloc(m);
  ch->A=gsl_matrix_calloc(m,m);

  fprintf(stderr,"CFACG4/DEBUG) allocated chapeau with %i knots on [%.5f,%.5f] dr %g \n",m,rmin,rmax,ch->dr);
  fflush(stderr);

  ch->alpha=0;

  return ch;
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
/* If m is the left bin of z (i.e. m=(int)((z-ch->rmin)/ch->dr) )
   give the fraction of dr to the nearest right bin:
  
     return           /-x-/
 o    x/dr           z        
 |               o   |        
 |       o       |   |   o    
 |       |       |   |   |    
 |       |       |   v   |    
 +-------+-------+-------+--- 
 0       1      m=2      3    
 rmin  rmin+dr 

*/
  return 1.0-ch->idr*(z-(m*ch->dr+zmin));
}

double left ( chapeau * ch, int m, double z, double zmin ) {
/* If m is the right bin of z (i.e. m=(int)((z-ch->rmin)/ch->dr)+1 )
   give the fraction of dr to the nearest right bin:
  
     return      /-x-/
 o    x/dr           z   
 |               o   |   
 |       o       |   |   o
 |       |       |   |   |
 |       |       |   v   |
 +-------+-------+-------+---
 0       1       2     m=3
 rmin  rmin+dr 

*/
  return ch->idr*(z-((m-1)*ch->dr+zmin));
}

void chapeau_setupoutput ( chapeau * ch, char * filename, int outputFreq, int outputLevel ) {
  char foo[255];
    
  if (ch) {
    char * flag="CxCx";
    ch->ofp=fopen(filename,"w");
    fwrite(flag,sizeof(char),4,ch->ofp);
    fwrite(&outputLevel,sizeof(int),1,ch->ofp);
    fwrite(&ch->m,sizeof(int),1,ch->ofp);
    ch->outputFreq=outputFreq;
    ch->outputLevel=outputLevel;

    strcpy(foo,filename);
    strcat(foo,".restart");
    ch->ofs=fopen(foo,"w");
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
	fwrite(&(ch->hits[i]),sizeof(long),1,ch->ofp);
      }
    }

    chapeau_savestate(ch,timestep);

    fflush(ch->ofp);
  }
}

void chapeau_savestate ( chapeau * ch, int timestep ) {
  if (ch&&ch->ofp&&!(timestep%ch->outputFreq)) {
    int N=ch->N;
    int m=ch->m;

    rewind(ch->ofs);
    fwrite(&timestep,sizeof(int),1,ch->ofs);

    fwrite(ch, sizeof(*ch), 1, ch->ofs);
    fwrite(ch->hits,sizeof(*(ch->hits)),m,ch->ofs);
    fwrite(ch->mask,sizeof(*(ch->mask)),N,ch->ofs);
    fwrite(ch->s,sizeof(***(ch->s)),m*N*3,ch->ofs);

    gsl_matrix_fwrite(ch->ofs,ch->A);
    gsl_vector_fwrite(ch->ofs,ch->b);
    gsl_vector_fwrite(ch->ofs,ch->lam);
    fflush(ch->ofs);
  }
}

chapeau * chapeau_allocloadstate ( char * filename ) {
    int i,m,N,d;
    FILE * ofs=fopen(filename,"r");
    chapeau * ch;
   
    fread(&i,sizeof(int),1,ofs);
    fprintf(stdout,"## Recovering chapeau state of step %i of some previous simulation\n",i);

    ch = (chapeau*)malloc(sizeof(chapeau));
  
    fread(ch, sizeof(*ch), 1, ofs);
    N=ch->N;
    m=ch->m;
   
    ch->hits=(long*)malloc(m*sizeof(long));
    fread(ch->hits,sizeof(*(ch->hits)),m,ofs);

    ch->mask=(int*)malloc(N*sizeof(int));
    fread(ch->mask,sizeof(*(ch->mask)),N,ofs);

    ch->s=(double***) malloc(N*sizeof(double**));
    for (i=0;i<N;i++) {
      ch->s[i]=(double**)malloc(3*sizeof(double*));
      for (d=0;d<3;d++) {
        ch->s[i][d]=(double*)malloc(m*sizeof(double));
        fread(ch->s[i][d],sizeof(*(ch->s[i][d])),m,ofs);
      }
    }

    ch->A=gsl_matrix_calloc(m,m);
    gsl_matrix_fread(ofs,ch->A);

    ch->b=gsl_vector_calloc(m);
    gsl_vector_fread(ofs,ch->b);

    ch->lam=gsl_vector_calloc(m);
    gsl_vector_fread(ofs,ch->lam);

    fread(&ch,sizeof(*ch),1,ofs);
    fclose(ofs);
 
    return ch;
}
    
void chapeau_init_particle_sums ( chapeau * ch ) {
  int d,i,j;
  int m=ch->m;
  for (j=0;j<ch->N;j++) for (d=0;d<3;d++) for (i=0;i<m;i++) ch->s[j][d][i]=0.0;
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
  int i,j,J,I,ii,jj;
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

    // parameters that are allowed to evolve lie between indices for which
    // ch->hits[] is non-zero so extract the proper subspace
    for (i=0;i<ch->m&&(!ch->hits[i]);i++);
    I=i;
    for (i=I;i<ch->m&&ch->hits[i];i++);
    J=i;

    nred=J-I; //minus the last
    Abar=gsl_matrix_alloc(nred,nred);
    bbar=gsl_vector_alloc(nred);
    lambar=gsl_vector_alloc(nred);
    p=gsl_permutation_alloc(nred);
    ii=0;
    for (i=I;i<J;i++) {
      gsl_vector_set(bbar,ii,gsl_vector_get(ch->b,i));
      jj=0;
      for (j=I;j<J;j++) {
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
