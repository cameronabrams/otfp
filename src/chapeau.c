#include "chapeau.h"

void chapeau_sum ( chapeau * ch1, chapeau * ch2 ) {
  //Overwrite ch1 with ch1+ch2
  int i,j,k;

  if (sizeof(ch1)!=sizeof(ch2)) {
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
  gsl_vector_add(ch1->bfull,ch2->bfull);
  gsl_matrix_add(ch1->Afull,ch2->Afull);

  for (i=0;i<ch1->m;i++) {
    ch1->hits[i] = (ch1->hits[i]||ch2->hits[i]);
  }

  //// I think that the next is irrelevant since I am saving the state before
  //// the chapeau_output and therefore before the chapeau_update, but I let this
  //// here just for the case.
  //for (i=0;i<ch1->N;i++) {
  //  for (j=0;j<3;j++) {
  //    for (k=0;k<ch1->m;k++) {
  //      ch1->s[i][j][k] = ch1->s[i][j][k] + ch2->s[i][j][k];
  //    }
  //  }
  //}
  //gsl_vector_add(ch1->lam ,ch2->lam); //also irrelevant

}

chapeau * chapeau_alloc ( int m, double rmin, double rmax, int npart ) {
  chapeau * ch;
  int d,i;

  ch=(chapeau*)malloc(sizeof(chapeau));

  if (npart<1) {
    fprintf(stderr,"CFACV/C) ERROR: chapeau objects without particles\n");
    exit(-1);
  }         

  if (m<1) {
    fprintf(stderr,"CFACV/C) ERROR: chapeau objects without bins\n");
    exit(-1);
  }         

  ch->m=m;
  //ch->mref=0;
  ch->N=npart;
  ch->rmin=rmin;
  ch->rmax=rmax;
  ch->dr=(rmax-rmin)/(m-1);
  ch->idr=1.0/ch->dr;
  ch->lam=gsl_vector_calloc(m);
  ch->nsample=0;
  ch->hits=(int*)calloc(m,sizeof(int));
  
  //ch->s=(double***) calloc(npart,sizeof(double**));
  //for (i=0;i<npart;i++) {
  //  ch->s[i]=(double**)calloc(3,sizeof(double*));
  //  for (d=0;d<3;d++) {
  //    ch->s[i][d]=(double*)calloc(m,sizeof(double));
  //  }
  //}

  ch->mask=(int*)calloc(npart,sizeof(int));
  
  ch->b=gsl_vector_calloc(m);
  ch->A=gsl_matrix_calloc(m,m);
  ch->bfull=gsl_vector_calloc(m);
  ch->Afull=gsl_matrix_calloc(m,m);

  fprintf(stdout,"CFACG4/DEBUG) allocated chapeau with %i knots on [%.5f,%.5f] dr %g \n",m,rmin,rmax,ch->dr);
  fflush(stdout);

  ch->alpha=0;

  return ch;
}

void chapeau_free ( chapeau * ch ) {
  free(ch->hits);
  free(ch->mask);
  gsl_vector_free(ch->b);
  gsl_matrix_free(ch->A);
  gsl_vector_free(ch->bfull);
  gsl_matrix_free(ch->Afull);
  free(ch);
}

int chapeau_quickcompare ( chapeau * ch1,  chapeau * ch2) {
  //this comparision routine does not compare any high rank member
  if (ch1->N           != ch2->N          ) return 0;
  if (ch1->m           != ch2->m          ) return 0;
  //if (ch1->nsample     != ch2->nsample    ) return 0;
  //if (ch1->nupdate     != ch2->nupdate    ) return 0;
  //if (ch1->alpha       != ch2->alpha      ) return 0;
  //if (ch1->outputFreq  != ch2->outputFreq ) return 0;
  //if (ch1->outputLevel != ch2->outputLevel) return 0;
  //if (ch1->rmin        != ch2->rmin       ) return 0;
  //if (ch1->rmax        != ch2->rmax       ) return 0;
  //if (ch1->dr          != ch2->dr         ) return 0;
  //if (ch1->idr         != ch2->idr        ) return 0;
  //if (strcmp(ch1->filename,ch2->filename)) return 0;
  return 1;
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
  if (ch&&ch->ofp&&!(timestep%ch->outputFreq)) {
    int outputlevel=ch->outputLevel;
    int i;
  
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

void chapeau_savestate ( chapeau * ch, int timestep, char * filename ) {
  int N=ch->N;
  int m=ch->m;
  FILE *ofs;

  ofs=fopen(filename,"w");

  fwrite(&timestep,sizeof(int),1,ofs);
  fwrite(ch, sizeof(*ch), 1, ofs);
  fwrite(ch->hits,sizeof(*(ch->hits)),m,ofs);
  fwrite(ch->mask,sizeof(*(ch->mask)),N,ofs);
  //fwrite(ch->s,sizeof(***(ch->s)),m*N*3,ofs);

  gsl_matrix_fwrite(ofs,ch->A);
  gsl_vector_fwrite(ofs,ch->b);
  gsl_matrix_fwrite(ofs,ch->Afull);
  gsl_vector_fwrite(ofs,ch->bfull);
  gsl_vector_fwrite(ofs,ch->lam);

  fclose(ofs);
  
}

chapeau * chapeau_allocloadstate ( char * filename ) {
    int i,m,N;
    FILE * ofs;
    chapeau * ch;
   
    fprintf(stdout,"CFACV) Allocating chapeau object from restart file\n");

    ch = (chapeau*)malloc(sizeof(chapeau));
  
    // Reading only the size
    ofs=fopen(filename,"r");
    fread(&i,sizeof(int),1,ofs);
    fread(ch, sizeof(*ch), 1, ofs);
    N=ch->N;
    m=ch->m;
    fclose(ofs);
   
    // Allocating
    ch->hits=(int*)malloc(m*sizeof(int));
    ch->mask=(int*)malloc(N*sizeof(int));
    ch->A=gsl_matrix_calloc(m,m);
    ch->b=gsl_vector_calloc(m);
    ch->Afull=gsl_matrix_calloc(m,m);
    ch->bfull=gsl_vector_calloc(m);
    ch->lam=gsl_vector_calloc(m);
    //ch->s=(double***) malloc(N*sizeof(double**));
    //for (i=0;i<N;i++) {
    //  ch->s[i]=(double**)malloc(3*sizeof(double*));
    //  for (d=0;d<3;d++) {
    //    ch->s[i][d]=(double*)malloc(m*sizeof(double));
    //    fread(ch->s[i][d],sizeof(*(ch->s[i][d])),m,ofs);
    //  }
    //}
  
    //BUG? fread(&ch,sizeof(*ch),1,ofs);

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

    // With chaux we prevent override of addresses (hits, mask, etc) when
    // reading integers (N,M,et), that is only what we want... better way to do
    // this?
    if (!chapeau_quickcompare(ch,chaux)) {
      fprintf(stdout,"CFACV) ERROR, the chapeau object from the file is diferent\n");
    }

    // This variables are alredy set during allocation, might be better to
    // compare them instead
    ch->nsample     = chaux->nsample    ;
    ch->rmin        = chaux->rmin       ;
    ch->rmax        = chaux->rmax       ;
    ch->dr          = chaux->dr         ;
    ch->idr         = chaux->idr        ;
    
    // This variables might be not set
    ch->nupdate     = chaux->nupdate    ;
    ch->alpha       = chaux->alpha      ;
    strcpy(ch->filename,chaux->filename);
    ch->outputFreq  = chaux->outputFreq ;
    ch->outputLevel = chaux->outputLevel;
    
    // Discarding all this readings by using chaux
    gsl_vector_memcpy(ch->lam,chaux->lam);
    chapeau_free(chaux);
}


void chapeau_loadstate ( chapeau * ch, char * filename ) {
    int i;
    chapeau * chaux = (chapeau*)malloc(sizeof(chapeau));
    FILE * ofs;
   
    ofs=fopen(filename,"r");

    if (!ch) {
      fprintf(stderr,"CFACV) ERROR in load chapeau file because holding object was not allocated\n");
      exit(-1);
    }

    fread(&i,sizeof(int),1,ofs);
    fprintf(stdout,"CFACV) Recovering chapeau state of step %i of some previous simulation\n",i);

    // With chaux we prevent override of addresses (hits, mask, etc) when
    // reading integers (N,M,et), that is only what we want... better way to do
    // this?
    fread(chaux, sizeof(*chaux), 1, ofs);
    if (!chapeau_quickcompare(ch,chaux)) {
      fprintf(stdout,"CFACV) ERROR, the chapeau object from the file is diferent\n");
    }
            
    // This variables are alredy set during allocation, might be better to
    // compare them instead
    ch->nsample     = chaux->nsample    ;
    ch->rmin        = chaux->rmin       ;
    ch->rmax        = chaux->rmax       ;
    ch->dr          = chaux->dr         ;
    ch->idr         = chaux->idr        ;
    
    // This variables might be not set
    ch->nupdate     = chaux->nupdate    ;
    ch->alpha       = chaux->alpha      ;
    strcpy(ch->filename,chaux->filename);
    ch->outputFreq  = chaux->outputFreq ;
    ch->outputLevel = chaux->outputLevel;
     
    fread(ch->hits,sizeof(*(ch->hits)),ch->m,ofs);
    fread(ch->mask,sizeof(*(ch->mask)),ch->N,ofs);

    gsl_matrix_fread(ofs,ch->A);
    gsl_vector_fread(ofs,ch->b);
    gsl_matrix_fread(ofs,ch->Afull);
    gsl_vector_fread(ofs,ch->bfull);
    gsl_vector_fread(ofs,ch->lam);

    fclose(ofs);
    free(chaux);
}
    
void chapeau_update_peaks ( chapeau * ch ) {
  int i,j,J,I;
  int s;
  double ninv;
  double lo,lb,alpha;

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

  //XXX: What is alpha?
  alpha=0.0;//1.0-exp(-1.0e4/nsamples);

  gsl_linalg_LU_decomp(Abar,p,&s);
  gsl_linalg_LU_solve(Abar,p,bbar,lambar);
  
  // update the vector of coefficients
  I=0;
  for (i=0;i<ch->m;i++) {
    if (!ch->hits[i]) continue;

    lb=gsl_vector_get(lambar,I);
    if (lb!=lb) {fprintf(stderr,"PARANOIA) Tripped at 3333. Too many chapeau additions?\n");exit(-1);}
    gsl_vector_set(ch->lam,i,lb);
    
    //XXX: What is alpha?
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
    
double chapeau_evalf ( chapeau * ch, double z ) {
  // Evaluate the free energy of a given CV vector from the linear
  // interpolation of the lambda vector on a chapeau object.
  // For now the interpolation is 1D.
  int m;
  double dm,la,lb;

  // Early return to avoid interpolations beyond the boundaries
  if ( z > ch->rmax ) return 0.;
  if ( z <= ch->rmin ) return 0.;
   
  // o                    z   
  // |                o   |   
  // |       o        |   |   o
  // |       |        |   |   |
  // |       |        |   |   |
  // min  min+dr ...  |   v   |
  // +-------+---//---+-------+---
  // 0       1   ...  m      m+1
  // /---------dm---------/

  dm=ch->idr*(z-ch->rmin);
  m=(int) dm;
  dm=dm-m;

  la=gsl_vector_get(ch->lam,m);
  lb=gsl_vector_get(ch->lam,m+1);
  return la+(lb-la)*dm; 

}  

//void chapeau_setmref ( chapeau * ch, double z ) {
//  // set the mref peak from the z rage and use it as baseline
//  int i;
//  double dm;
//  double la,lb;
//   
//  // Early return to avoid interpolations beyond the boundaries
//  if ( z > ch->rmax ) return;
//  if ( z <= ch->rmin ) return;
//   
//  // o                    z   
//  // |                o   |   
//  // |       o        |   |   o
//  // |       |        |   |   |
//  // |       |        |   |   |
//  // min  min+dr ...  |   v   |
//  // +-------+---//---+-------+---
//  // 0       1   ...  m      m+1
//  // /---------dm---------/
//
//  dm=ch->idr*(z-ch->rmin);
//  ch->mref=(int) dm; 
//
//  la=gsl_vector_get(ch->lam,ch->mref);
//  for (i=0;i<ch->m;i++) {
//    lb=gsl_vector_get(ch->lam,i);
//    gsl_vector_set(ch->lam,i,lb-la);
//  }  
//}  
//     
//void chapeau_baselinehits ( chapeau * ch) {
//  // Add a constant to the lam vector in order to set the baseline in the
//  // element k.
//  int i,m;
//  double dm,la,lb;
//   
//  la=gsl_vector_get(ch->lam,ch->mref);
//  for (i=0;i<ch->m;i++) {
//    if (!ch->hits[i]) continue;
//    lb=gsl_vector_get(ch->lam,i);
//    gsl_vector_set(ch->lam,i,lb-la);
//  }  
//}  
    
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
 
