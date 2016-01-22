#include "chapeau.h"

chapeau * chapeau_alloc ( int dm, double * rmin, double * rmax, int * N ) {
  chapeau * ch;
  int j,i;

  //TODO: check size of argument vectors and write the errors
  
  // Allocating the object
  ch=(chapeau*)malloc(sizeof(chapeau));
  ch->ofp=0;

  // Allocating the dimension
  ch->dm=dm;
  ch->rmin=(double*)malloc(dm*sizeof(double));
  ch->rmax=(double*)malloc(dm*sizeof(double));
  ch->dr=(double*)malloc(dm*sizeof(double));
  ch->idr=(double*)malloc(dm*sizeof(double));
  ch->N=(int*)malloc(dm*sizeof(int));
  ch->r=(double*)malloc(dm*sizeof(double));
  ch->f=(double*)malloc(dm*sizeof(double));

  if (ch->dm==1) {
    ch->accumulate=accumulate_1D;
  } else {
    ch->accumulate=accumulate_2D;
  }

  ch->m=1;
  for (i=0;i<dm;i++) {
    ch->N[i]=N[i];   
    ch->rmin[i]=rmin[i];
    ch->rmax[i]=rmax[i];
    ch->dr[i]  =(rmax[i]-rmin[i])/(ch->N[i]-1);
    ch->idr[i] =1./ch->dr[i];
    ch->m=ch->m*N[i];
  }

  // Allocating the size
  ch->lam=gsl_vector_calloc(ch->m);
  ch->hits=(int*)calloc(ch->m,sizeof(int));
  ch->b=gsl_vector_calloc(ch->m);
  ch->A=gsl_matrix_calloc(ch->m,ch->m);
  ch->bfull=gsl_vector_calloc(ch->m);
  ch->Afull=gsl_matrix_calloc(ch->m,ch->m);

  return ch;
}

void chapeau_free ( chapeau * ch ) {
  gsl_matrix_free(ch->A);
  gsl_vector_free(ch->b);
  gsl_matrix_free(ch->Afull);
  gsl_vector_free(ch->bfull);
  free(ch->hits);
  free(ch->rmin);
  free(ch->rmax);
  free(ch->dr);
  free(ch->idr);
  free(ch->N);
  free(ch);
}
             
int chapeau_comparesize ( chapeau * ch1,  chapeau * ch2) {
  int i;
  if (ch1->dm          != ch2->dm         ) return 0;
  if (ch1->m           != ch2->m          ) return 0;
  return 1;
}
 
int chapeau_comparegrid ( chapeau * ch1,  chapeau * ch2) {
  int i;
  //if (ch1->nupdate     != ch2->nupdate    ) return 0;
  //if (ch1->outputFreq  != ch2->outputFreq ) return 0;
  //if (ch1->outputLevel != ch2->outputLevel) return 0;
  //if (strcmp(ch1->filename,ch2->filename)) return 0;
  if (ch1->dm          != ch2->dm         ) return 0;
  for (i=0;i<ch2->dm;i++) {
    if (ch1->rmin[i]  != ch2->rmin[i] ) return 0;
    if (ch1->rmax[i]  != ch2->rmax[i] ) return 0;
    if (ch1->dr[i]    != ch2->dr[i]   ) return 0;
    if (ch1->idr[i]   != ch2->idr[i]  ) return 0;
    if (ch1->N[i]  != ch2->N[i] ) return 0;
  }
  return 1;
}

void chapeau_sum ( chapeau * ch1, chapeau * ch2 ) {
  //Overwrite ch1 with ch1+ch2
  int i,j,k;

  if (!chapeau_comparesize(ch1,ch2)) {
    fprintf(stderr,"CFACV/C) ERROR: you can not sum chapeau objects with different sizes\n");
    exit(-1);
  }
  
  //Really?
  //if (!chapeau_comparegrid(ch1,ch2)) {
  //  fprintf(stderr,"CFACV/C) ERROR: you can not sum chapeau objects with different domains\n");
  //  exit(-1);
  //}

  gsl_vector_add(ch1->b   ,ch2->b);
  gsl_matrix_add(ch1->A   ,ch2->A);
  gsl_vector_add(ch1->bfull,ch2->bfull);
  gsl_matrix_add(ch1->Afull,ch2->Afull);

  for (i=0;i<ch1->m;i++) {
    ch1->hits[i] = (ch1->hits[i]||ch2->hits[i]);
  }

  //gsl_vector_add(ch1->lam ,ch2->lam); // irrelevant?

}
        
void chapeau_setupoutput ( chapeau * ch, char * outfile, char * restartfile, int outputFreq, int outputLevel ) {

  if (!ch) exit(-1);

  // output
  ch->ofp=fopen(outfile,"w");
  char * flag="CxCx";
  fwrite(flag,sizeof(char),4,ch->ofp);
  fwrite(&outputLevel,sizeof(int),1,ch->ofp);
  fwrite(&ch->m,sizeof(int),1,ch->ofp);
  ch->outputFreq=outputFreq;
  ch->outputLevel=outputLevel;

  // restart
  strcpy(ch->filename,restartfile);
    
}

void chapeau_output ( chapeau * ch, int timestep ) {
  int outputlevel=ch->outputLevel;
  int i;
  
  if ( !ch || !ch->ofp || timestep%ch->outputFreq ) return;

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

void chapeau_savestate ( chapeau * ch, char * filename ) {
  FILE *ofs;

  ofs=fopen(filename,"w");

  fwrite(ch, sizeof(*ch), 1, ofs);

  fwrite(ch->rmin,sizeof(*(ch->rmin)),ch->dm,ofs);
  fwrite(ch->rmax,sizeof(*(ch->rmax)),ch->dm,ofs);
  fwrite(ch->N,sizeof(*(ch->N)),ch->dm,ofs);
   
  gsl_vector_fwrite(ofs,ch->lam);
  fwrite(ch->hits,sizeof(*(ch->hits)),ch->m,ofs);
  gsl_matrix_fwrite(ofs,ch->A);
  gsl_vector_fwrite(ofs,ch->b);
  gsl_matrix_fwrite(ofs,ch->Afull);
  gsl_vector_fwrite(ofs,ch->bfull);

  fclose(ofs);
  
}

void chapeau_loadstate ( chapeau * ch, char * filename ) {
    chapeau * chaux;
    int i;
    FILE * ofs;
   
    ofs=fopen(filename,"r");

    if (!ch) {
      fprintf(stderr,"CFACV) ERROR in load chapeau file because holding object was not allocated\n");
      exit(-1);
    }

    // With chaux we prevent override of addresses (hits, A, etc) when
    // reading integers (N,M,et), that is only what we want... better way to do
    // this?
    chaux = (chapeau*)malloc(sizeof(chapeau));
    fread(chaux, sizeof(*chaux), 1, ofs);
    if (!chapeau_comparesize(ch,chaux)) {
      fprintf(stdout,"CFACV) ERROR, the chapeau object from the file is diferent\n");
    }
            
    // This variables are alredy set during allocation, might be better to
    // compare them instead
    ch->dm      = chaux->dm;
    ch->m       = chaux->m;
    
    // This variables might be not set
    ch->nupdate     = chaux->nupdate    ;
    strcpy(ch->filename,chaux->filename);
    ch->outputFreq  = chaux->outputFreq ;
    ch->outputLevel = chaux->outputLevel;
    
    // Read directly
    fread(ch->rmin,sizeof(*(ch->rmin)),ch->dm,ofs);
    fread(ch->rmax,sizeof(*(ch->rmax)),ch->dm,ofs);
    fread(ch->N,sizeof(*(ch->N)),ch->dm,ofs);

    gsl_vector_fread(ofs,ch->lam);
    fread(ch->hits,sizeof(*(ch->hits)),ch->m,ofs);
    gsl_matrix_fread(ofs,ch->A);
    gsl_vector_fread(ofs,ch->b);
    gsl_matrix_fread(ofs,ch->Afull);
    gsl_vector_fread(ofs,ch->bfull);

    fclose(ofs);
    free(chaux);
}
    
chapeau * chapeau_allocloadstate ( char * filename ) {
    int i,m,N;
    FILE * ofs;
    chapeau * ch;
   
    fprintf(stdout,"CFACV) Allocating chapeau object from restart file\n");

    ch = (chapeau*)malloc(sizeof(chapeau));
  
    // Reading only the sizes
    ofs=fopen(filename,"r");
    fread(&i,sizeof(int),1,ofs);
    fread(ch, sizeof(*ch), 1, ofs);
    fclose(ofs);
 
    // Allocating the dimension
    ch->rmin=(double*)malloc(ch->dm*sizeof(double));
    ch->rmax=(double*)malloc(ch->dm*sizeof(double));
    ch->dr=(double*)malloc(ch->dm*sizeof(double));
    ch->idr=(double*)malloc(ch->dm*sizeof(double));
    ch->N=(int*)malloc(ch->dm*sizeof(int));

    // Allocating the size
    ch->lam=gsl_vector_calloc(ch->m);
    ch->hits=(int*)calloc(ch->m,sizeof(int));
    ch->b=gsl_vector_calloc(ch->m);
    ch->A=gsl_matrix_calloc(ch->m,ch->m);
    ch->bfull=gsl_vector_calloc(ch->m);
    ch->Afull=gsl_matrix_calloc(ch->m,ch->m);
 
    // Now reading the state
    chapeau_loadstate(ch,filename);

    return ch;
}

void chapeau_loadlambda ( chapeau * ch, char * filename ) {
    int i,m,N;
    chapeau * chaux;
    FILE * ofs;
 
    if (!ch) {
      fprintf(stderr,"CFACV) ERROR in load chapeau file because holding object was not allocated\n");
      exit(-1);
    }
                                           
    chaux=chapeau_allocloadstate(filename);

    // With chaux we prevent override of addresses (hits, A, etc) when
    // reading integers (N,M,et), that is only what we want... better way to do
    // this?
    if (!chapeau_comparesize(ch,chaux)) {
      fprintf(stdout,"CFACV) ERROR, the chapeau object from the file is diferent\n");
    }
            
    // This variables are alredy set during allocation, might be better to
    // compare them instead
    ch->dm      = chaux->dm;
    ch->m       = chaux->m;
                                                 
    for (i=0;i<ch->dm;i++) {
      ch->rmin[i] = chaux->rmin[i];
      ch->rmax[i] = chaux->rmax[i];
      ch->N[i] = chaux->N[i];   
      ch->dr[i] = (ch->rmax[i]-ch->rmin[i])/(ch->m-1);
      ch->idr[i] = 1./ch->dr[i];
    } 
           
    // This variables might be not set
    ch->nupdate     = chaux->nupdate    ;
    strcpy(ch->filename,chaux->filename);
    ch->outputFreq  = chaux->outputFreq ;
    ch->outputLevel = chaux->outputLevel;
    
    // Discarding all this readings by using chaux
    gsl_vector_memcpy(ch->lam,chaux->lam);
    chapeau_free(chaux);
}

void chapeau_update_peaks ( chapeau * ch ) {
  int i,j,J,I;
  int s;
  double ninv;
  double lo,lb;
  gsl_matrix * Abar;
  gsl_vector * bbar;
  gsl_vector * lambar;
  gsl_permutation * p;
  int nred;

  //DB//for (i=0;i<ch->m;i++) fprintf(stderr,"AAA %i %i\n",i,ch->hits[i]); exit(1);

  //// parameters that are allowed to evolve lie between indices for which
  //// ch->hits[] is non-zero so extract the proper subspace
  nred=0;
  for (i=0;i<ch->m;i++) if (ch->hits[i]) nred++;

  // If the nred is small, singular matrix can occur
  //fprintf(stderr,"Warning: No chapeau update: only %i non cero elements\n",nred);
  if (nred<10) return;

  Abar=gsl_matrix_alloc(nred,nred);
  bbar=gsl_vector_alloc(nred);
  lambar=gsl_vector_alloc(nred);
  p=gsl_permutation_alloc(nred);
  
  // Add the new and old statistic to the reduced matrix
  I=0;
  for (i=0;i<ch->m;i++) {
    if (!ch->hits[i]) continue;
    lb=gsl_vector_get(ch->b,i)+gsl_vector_get(ch->bfull,i);
    gsl_vector_set(bbar,I,-lb); //TODO: Trace back the origin of the minus. See fes1D procedure. 
    J=0;
    for (j=0;j<ch->m;j++) {
      if (!ch->hits[j]) continue;
      lb=gsl_matrix_get(ch->A,i,j)+gsl_matrix_get(ch->Afull,i,j);
      gsl_matrix_set(Abar,I,J,lb);
      J++;
    }
    I++;
  } 

  gsl_linalg_LU_decomp(Abar,p,&s);
  gsl_linalg_LU_solve(Abar,p,bbar,lambar);
  
  // update the vector of coefficients
  I=0;
  for (i=0;i<ch->m;i++) {
    if (!ch->hits[i]) continue;

    lb=gsl_vector_get(lambar,I);
    if (lb!=lb) {fprintf(stderr,"PARANOIA) Tripped at 3333. Too many chapeau additions?\n");exit(-1);}
    gsl_vector_set(ch->lam,i,lb);
    
    //lo=gsl_vector_get(ch->lam,i);
    //gsl_vector_set(ch->lam,i,alpha*lo+(1-alpha)*lb);

    I++;
  } 
   
  gsl_matrix_free(Abar);
  gsl_vector_free(bbar);
  gsl_vector_free(lambar);
  gsl_permutation_free(p);

  //chapeau_baselinehits(ch); 
}

void chapeau_set_peaks ( chapeau * ch, char * filename ) {
  FILE * fp = fopen(filename,"r");
  double knots;
  int i=0;
  char ln[255];

  fprintf(stdout,"INFO) Read knots from %s\n",filename);

  for (i=0;i<ch->m;i++) {
    fgets(ln,255,fp);
    sscanf(ln,"%lf",&knots);
    fprintf(stdout,"INFO) %i,%s\n",i++,knots);
    gsl_vector_set(ch->lam,i,knots);
  }

  fflush(stdout);
}

// Process for replica exechange
    
double chapeau_evalf_1simplex ( chapeau * ch, double z ) {
  // Evaluate the free energy of a given CV vector from the linear
  // interpolation of the lambda vector on a chapeau object.

  int m;
  double dm,la,lb;

  // Early return to avoid interpolations beyond the boundaries
  if (ch->dm > 1 ) {
    fprintf(stderr,"CFACV/C) ERROR: chapeau object does not leaves in 1D");
    exit(-1);
  }
  if ( z > ch->rmax[0] ) return 0.;
  if ( z <= ch->rmin[0] ) return 0.;
   
  /*
    o                    z   
    |                o   |   
    |       o        |   |   o
    |       |        |   |   |
    |       |        |   |   |
    min  min+dr ...  |   v   |
    +-------+---//---+-------+---
    0       1   ...  m      m+1
    /---------dm---------/
  */

  dm=ch->idr[0]*(z-ch->rmin[0]);
  m=(int) dm;
  dm=dm-m;

  la=gsl_vector_get(ch->lam,m);
  lb=gsl_vector_get(ch->lam,m+1);
  return la+(lb-la)*dm; 

} 
    
double chapeau_evalf_2simplex ( chapeau * ch, double z1, double z2 ) {
  int m;
  int i,j,k,ni,nj,nk;
  double f,dx,dy,lb;

  // Early return to avoid interpolations beyond the boundaries
  if ( z1 > ch->rmax[0] ) return 0.;
  if ( z1 <= ch->rmin[0] ) return 0.;
  if ( z2 > ch->rmax[1] ) return 0.;
  if ( z2 <= ch->rmin[1] ) return 0.;
  f=0;

  /*
 
  -\  |      -\  |      -\  |    
    -\|        -\|        -\|    
  ----*----------*----------*----
      |-\        |-\        |-\
      |  -\  R5  |  -\      |  -\
      |    -\    |    -\    |    
  -\  |  R6  -\  | R4   -\  |    
    -\|        -\|        -\|    
  ----*----------*----------*----
      |-\        |-\        |-\
      |  -\  R1  |  -\  R3  |  -\
      |    -\    |    -\    |    
  -\  |      -\  |  R2  -\  |    
    -\|        -\|        -\|    
  ----*----------*----------*----
      |-\        |-\        |-\
      |  -\      |  -\      |  -\

  */
 
  // Identifico el cuadrado, esto define dos nodos
  i=(int)( (z1 - ch->rmin[0]) * ch->idr[0] );
  j=(int)( (z2 - ch->rmin[1]) * ch->idr[1] );
  ni=j*(ch->N[0]-1)+i+1;
  nj=j+1*(ch->N[0]-1)+i;
    
  // ahora el triangulo, define el tercer nodo
  dx=z1 - i*ch->idr[0] -1;
  dy=z2 - j*ch->idr[1];
  if(dy/dx<1) {
    nk=j*(ch->N[0]-1)+i;
    /*nj----------* 
       |-\        | 
       |  -\      | 
       |    -\    | 
       |  x   -\  | 
       |        -\| 
      nk---------ni */

    lb=gsl_vector_get(ch->lam,nk);
    f-=(dx*ch->idr[0]+1)*lb;
    f-=(dy*ch->idr[1]+1)*lb;

  } else {
    nk=(j+1)*(ch->N[0]-1)+(i+1);
    /*nj---------nk
       |-\        |  
       |  -\   x  |  
       |    -\    |  
       |      -\  |  
       |        -\|  
       *---------ni */

    lb=gsl_vector_get(ch->lam,nk);
    f+=(dx*ch->idr[0]+1)*lb;
    f+=(dy*ch->idr[1]+1)*lb;

  }
    
  // Sumando contribucion de ni y nj
  lb=gsl_vector_get(ch->lam,ni);
  f+=(dx*ch->idr[0]+1)*lb;
  lb=gsl_vector_get(ch->lam,ni);
  f+=(dy*ch->idr[1]+1)*lb;

  return f;
}  

char * chapeau_serialize ( chapeau * ch ) {
  // return a str containing the serialized chapeau with partial statistics. If
  // this is the central replica, or non replica scheme is used, this serialize
  // the full statistic.
  int i;
  int size;
  char* buf;

  size = (3*ch->m-1)*13+ch->m+1;
  buf  = (char*) malloc(size*sizeof(char));

  // Seems that if I do this:
  //  return (char*)&ch
  //  is not portable, since the way that the cast is made is undefined
  //  (depends of the machine). So it is needed a proper serialization
  
  size=0;
  for (i=0;i<ch->m;i++) {

    // Writing the partial b
    sprintf(buf+size, "%13.5e", gsl_vector_get(ch->b,i)); size+=13;

    // Writing hits
    sprintf(buf+size, "%1i", ch->hits[i]); size+=1;
  }

  // Writing the partial A (simetric and tridiagonal)
  sprintf(buf+size, "%13.5e",gsl_matrix_get(ch->A,0,0)); size+=13;
  for (i=1;i<ch->m;i++) {
    sprintf(buf+size, "%13.5e",gsl_matrix_get(ch->A,i,i)); size+=13;
    sprintf(buf+size, "%13.5e",gsl_matrix_get(ch->A,i,i-1)); size+=13; 
  }

  return buf; 
}     

void chapeau_addserialized ( chapeau *ch, char * str ) {
  // str contains the serialized chapeau that comes from the partial sampling
  // of other replica partial statistic information. This is added to the
  // partial sampling ot the principal replica computing the sum.
  int i,err;
  int size1;
  int size2;
  double aux;
  int iaux;
  char word1[14],word2[2];
  
  // Add null terminators
  word1[13]='\0';
  word2[1]='\0';

  size1 = (3*ch->m-1)*13+ch->m;
  size2 = strlen(str);
  if (size1!=size2) {
    fprintf(stderr,"CFACV/C) ERROR: you can not sum chapeaus with different sizes %d %d\n",size1,size2);
    exit(-1);
  }
 
  // Warning, sscanf, strtod, strtoi, all start in the first nonblank
  // position, then there is not a proper way to read the serialized object
  // using sscanf.
  //
  // TODO, strtod and strtoi more efficient that sscanf?

  size1=0;
  for (i=0;i<ch->m;i++) {

    // Reading the partial b
    memcpy(word1, &str[size1], 13 ); size1+=13;
    err=sscanf(word1,"%13le",&aux); if(!err) {fprintf(stderr,"CFACV/C) Error 1112 on read %s\n",word1);}
    aux+=gsl_vector_get(ch->b,i);
    gsl_vector_set(ch->b,i,aux);

    // Reading hits
    memcpy(word2, &str[size1], 1 ); size1+=1;
    err=sscanf(word2,"%1i",&iaux); if(!err) {fprintf(stderr,"CFACV/C) Error 1111 on read %s\n",word2);}
    ch->hits[i]=(ch->hits[i]||iaux);
  }
 
  // Reading the partial A (simetric and tridiagonal)
  memcpy(word1, &str[size1], 13 ); size1+=13;
  err=sscanf(word1,"%13le",&aux); if(!err) {fprintf(stderr,"CFACV/C) Error 1112 on read %s\n",word1);}
  aux+=gsl_matrix_get(ch->A,0,0);
  gsl_matrix_set(ch->A,0,0,aux);

  for (i=1;i<ch->m;i++) {

    //diagonal
    memcpy(word1, &str[size1], 13 ); size1+=13;
    err=sscanf(word1,"%13le",&aux); if(!err) {fprintf(stderr,"CFACV/C) Error 1113 on read %s\n",word1);}
    aux+=gsl_matrix_get(ch->A,i,i);
    gsl_matrix_set(ch->A,i,i,aux);              

    //offdiagonal
    memcpy(word1, &str[size1], 13 ); size1+=13;
    err=sscanf(word1,"%13le",&aux); if(!err) {fprintf(stderr,"CFACV/C) Error 1114 on read %s\n",word1);}
    aux+=gsl_matrix_get(ch->A,i,i-1);
    gsl_matrix_set(ch->A,i,i-1,aux);
    gsl_matrix_set(ch->A,i-1,i,aux);
  }
   

  //TODO, give error if rmin and rmax are different
  
}

void chapeau_setserialized ( chapeau *ch, char * str ) {
  // str contains the full statistics information that comes from all the
  // replicas contributions. Therefore, this should be stored in Afull and
  // bfull and the partial A and b should be reset to cero.
  int i,err;
  int size1;
  int size2;
  double aux;
  int iaux;
  char word1[14],word2[2];
  
  // Add null terminators
  word1[13]='\0';
  word2[1]='\0';
  
  size1 = (3*ch->m-1)*13+ch->m;
  size2 = strlen(str);
  if (size1!=size2) {
    fprintf(stderr,"CFACV/C) ERROR: you can not set chapeaus with different sizes %d %d\n",size1,size2);
    exit(-1);
  }
  
  // Warning, sscanf, strtod, strtoi, all start in the first nonblank
  // position, then there is not a proper way to read the serialized object
  // using sscanf.
  //
  // TODO, strtod and strtoi more efficient that sscanf?
    
  size1=0;
  for (i=0;i<ch->m;i++) {

    // Read b
    memcpy(word1, &str[size1], 13 ); size1+=13;
    err=sscanf(word1,"%13le",&aux); if(!err) {fprintf(stderr,"CFACV/C) Error 2111 on read %s\n",word1);}
    gsl_vector_set(ch->bfull,i,aux);

    // Read hits
    memcpy(word2, &str[size1], 1 ); size1+=1;
    err=sscanf(word2,"%1i",&iaux); if(!err) {fprintf(stderr,"CFACV/C) Error 2112 on read %s\n",word1);}
    ch->hits[i]=iaux;
  }
     
  // Read A (simetric and tridiagonal)
  memcpy(word1, &str[size1], 13 ); size1+=13; 
  err=sscanf(word1,"%13le",&aux); if(!err) {fprintf(stderr,"CFACV/C) Error 2113 on read %s\n",word1);}
  gsl_matrix_set(ch->Afull,0,0,aux);

  for (i=1;i<ch->m;i++) {

    //diagonal
    memcpy(word1, &str[size1], 13 ); size1+=13;
    err=sscanf(word1,"%13le",&aux); if(!err) {fprintf(stderr,"CFACV/C) Error 2114 on read %s\n",word1);}
    gsl_matrix_set(ch->Afull,i,i,aux);

    //off diagonal
    memcpy(word1, &str[size1], 13 ); size1+=13;
    err=sscanf(word1,"%13le",&aux); if(!err) {fprintf(stderr,"CFACV/C) Error 2115 on read %s\n",word1);}
    gsl_matrix_set(ch->Afull,i,i-1,aux);
    gsl_matrix_set(ch->Afull,i-1,i,aux);
  }
  
  // Reset the variables to hold the partial statistic information
  gsl_vector_set_zero(ch->b);
  gsl_matrix_set_zero(ch->A);

  //TODO, give error if rmin and rmax are different
}

 int accumulate_1D( chapeau * ch ) { 
  int i,j,k,m;
  double F;
  
  // Following the paper E definition A and b should be
  //
  // A=1./nsteps \sum_{steps} \sum_{i} dphi_m(z)/dz dphi_n(z)/dz
  // b=-1./nsteps \sum_{steps} \sum_{i} dphi_m(z)/dz F
  //
  // But here we compute -b*nsteps and A*nsteps and then A and b are scaled
  // by the proper factor when equation A\lambda=b is solved (see the
  // chapeau_update_peaks procedure)
  
  // Early return to avoid interpolations out of the domain
  if ( ch->r[0] > ch->rmax[0] ) return;
  if ( ch->r[0] <= ch->rmin[0] ) return;

  //ch->nsample++;

  /* Interpolation R->m
        m R
        | |
    o---o-x-o---o--*/
  m=(int)((ch->r[0]-ch->rmin[0])*ch->idr[0]);
  ch->hits[m]=1;
  ch->hits[m+1]=1; // Testing

  /* 1/dz  /| 
          / |n=m+1
         /  | 
    o---o-x-o---o--*/ 
  ch->b->data[m+1]+=ch->idr[0]*ch->f[0];
  ch->A->data[(m+1)*ch->A->tda+(m+1)]+=ch->idr[0]*ch->idr[0];

  /*    |\   -1/dz
     n=m| \
        |  \
    o---o-x-o---o--*/ 
  ch->b->data[m]-=ch->idr[0]*ch->f[0];
  ch->A->data[m*ch->A->tda+m]+=ch->idr[0]*ch->idr[0];

  /*    |\ /| 
       m| X |n=m+1
        |/ \| 
    o---o-x-o---o--*/ 
  ch->A->data[(m+1)*ch->A->tda+m]-=ch->idr[0]*ch->idr[0];
  ch->A->data[m*ch->A->tda+m+1]-=ch->idr[0]*ch->idr[0];

  return 0;
}
 
int accumulate_2D( chapeau * ch ) { 
  int m;
  int i,j,k,ni,nj,nk;
  double dx,dy;
  
  // Following the paper E definition A and b should be
  //
  // A=1./nsteps \sum_{steps} \sum_{i} dphi_m(z)/dz dphi_n(z)/dz
  // b=-1./nsteps \sum_{steps} \sum_{i} dphi_m(z)/dz F
  //
  // But here we compute -b*nsteps and A*nsteps and then A and b are scaled
  // by the proper factor when equation A\lambda=b is solved (see the
  // chapeau_update_peaks procedure)
 
  // Early return to avoid interpolations out of the domain
  if ( ch->r[0] > ch->rmax[0] ) return;
  if ( ch->r[1] > ch->rmax[1] ) return;
  if ( ch->r[0] <= ch->rmin[0] ) return;
  if ( ch->r[1] <= ch->rmin[1] ) return;

  /* LA BASE
    
  Li, Z., Qiao, Z., & Tang, T. (2015). The Finite Element Method for 2D
  Problems. In Numerical Solutions of Partial differential equations- An
  introduction to finite difference and finite element methods (pp. 199–245).
  Url: http://www4.ncsu.edu/~zhilin/TEACHING/MA587/

  -\  |      -\  |      -\  |    
    -\|        -\|        -\|    
  ----*----------*----------*----
      |-\        |-\        |-\
      |  -\  R5  |  -\      |  -\
      |    -\    |    -\    |    
  -\  |  R6  -\  | R4   -\  |    
    -\|        -\|        -\|    
  ----*----------*----------*----
      |-\        |-\        |-\
      |  -\  R1  |  -\  R3  |  -\
      |    -\    |    -\    |    
  -\  |      -\  |  R2  -\  |    
    -\|        -\|        -\|    
  ----*----------*----------*----
      |-\        |-\        |-\
      |  -\      |  -\      |  -\
  
  */

  // Identifico el cuadrado, esto define dos nodos
  /*nj----------*
     |          |  
     |          |  
     |     x    |  
     |          |  
     |          |  
     *---------ni */
  dx=(ch->r[0]-ch->rmin[0])*ch->idr[0]; 
  dy=(ch->r[1]-ch->rmin[1])*ch->idr[1]; 

  i=(int)(dx);
  j=(int)(dy);
  ni=j*ch->N[0]+i+1;
  //nj=(j+1)*ch->N[0]+i;
  nj=ni+ch->N[0]-1;
   
  // ahora el triangulo, define el tercer nodo
  if(dy<(j+1)-(dx-i)) {

    //nk=j*ch->N[0]+i;
    nk=ni-1;

    /*nj----------* 
       |-\        | 
       |  -\      | 
       |    -\    | 
       |  x   -\  | 
       |        -\| 
      nk---------ni */
    
    ch->b->data[ni]               +=ch->idr[0]*ch->f[0];
    ch->A->data[ni*ch->A->tda+ni] +=ch->idr[0]*ch->idr[0];

    ch->b->data[nj]               += ch->idr[1]*ch->f[1];
    ch->A->data[nj*ch->A->tda+nj] += ch->idr[1]*ch->idr[1];
     
    ch->b->data[nk]               -= ch->idr[0]*ch->f[0];
    ch->b->data[nk]               -= ch->idr[1]*ch->f[1];
    ch->A->data[nk*ch->A->tda+nk] += ch->idr[0]*ch->idr[0];
    ch->A->data[nk*ch->A->tda+nk] += ch->idr[1]*ch->idr[1];
    
    // Aca podria evitar 2 operaciones si eligiera el triangulo superior
    ch->A->data[nj*ch->A->tda+nk] -= ch->idr[1]*ch->idr[1];
    ch->A->data[nk*ch->A->tda+nj] -= ch->idr[1]*ch->idr[1];
    ch->A->data[ni*ch->A->tda+nk] -= ch->idr[0]*ch->idr[0];
    ch->A->data[nk*ch->A->tda+ni] -= ch->idr[0]*ch->idr[0];

  } else {
    //nk=(j+1)*ch->N[0]+i+1;
    nk=nj+1;
 
    /*nj---------nk
       |-\        |  
       |  -\   x  |  
       |    -\    |  
       |      -\  |  
       |        -\|  
       *---------ni */

    ch->b->data[ni]               -=ch->idr[1]*ch->f[1];
    ch->A->data[ni*ch->A->tda+ni] +=ch->idr[1]*ch->idr[1];

    ch->b->data[nj]               -= ch->idr[0]*ch->f[0];
    ch->A->data[nj*ch->A->tda+nj] += ch->idr[0]*ch->idr[0];
     
    ch->b->data[nk]               += ch->idr[0]*ch->f[0];
    ch->b->data[nk]               += ch->idr[1]*ch->f[1];
    ch->A->data[nk*ch->A->tda+nk] += ch->idr[0]*ch->idr[0];
    ch->A->data[nk*ch->A->tda+nk] += ch->idr[1]*ch->idr[1];
       
    // Aca podria evitar 2 operaciones si eligiera el triangulo superior
    ch->A->data[nk*ch->A->tda+nj] -= ch->idr[0]*ch->idr[0];
    ch->A->data[nj*ch->A->tda+nk] -= ch->idr[0]*ch->idr[0];
    ch->A->data[nk*ch->A->tda+ni] -= ch->idr[1]*ch->idr[1];
    ch->A->data[ni*ch->A->tda+nk] -= ch->idr[1]*ch->idr[1];
        
  }
 
  ch->hits[ni]=1;
  ch->hits[nj]=1;
  ch->hits[nk]=1;

  return 0;
}
 

