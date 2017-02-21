#include "chapeau.h"

chapeau * chapeau_alloc ( int dm, double * rmin, double * rmax, int * N, int periodic ) {
  chapeau * ch;
  int i;

  //TODO: check size of argument vectors and write the errors
  
  // Allocating the object
  ch=(chapeau*)malloc(sizeof(chapeau));
  ch->ofp=0;

  // Allocating the dimension
  ch->dm=dm;
  ch->periodic=periodic;
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

  // chapeau 1D => 1 uperdiagonals
  // chapeau 2D => ch->N[0] uperdiagonals
  if (dm==1) {
    ch->ku=1;
  } else{
    ch->ku=ch->N[0];
  }  
  ch->ldad=2*ch->ku+1;

  // Allocating the size
  ch->lam=(double*)calloc(ch->m,sizeof(double));
  ch->hits=(int*)calloc(ch->m,sizeof(int));
  ch->b=(double*)calloc(ch->m,sizeof(double));
  ch->bfull=(double*)calloc(ch->m,sizeof(double));

  ch->A=(double**)calloc(ch->ldad,sizeof(double*));
  ch->Afull=(double**)calloc(ch->ldad,sizeof(double*));
  for (i=0;i<ch->ldad;i++) {
    ch->A[i]=(double*)calloc(ch->m,sizeof(double));
    ch->Afull[i]=(double*)calloc(ch->m,sizeof(double));
  }

  return ch;
}

void chapeau_free ( chapeau * ch ) {
  int i;
  for (i = 0; i < ch->ldad; i++) { 
      free(ch->A[i]);
      free(ch->Afull[i]);
  }
  free(ch->A);
  free(ch->Afull);
  free(ch->b);
  free(ch->bfull);
  free(ch->hits);
  free(ch->rmin);
  free(ch->rmax);
  free(ch->dr);
  free(ch->idr);
  free(ch->N);
  free(ch);
}
             
int chapeau_comparesize ( chapeau * ch1,  chapeau * ch2) {
  if (ch1->dm          != ch2->dm         ) return 0;
  if (ch1->m        != ch2->m       ) return 0;
  if (ch1->ldad        != ch2->ldad       ) return 0;
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
  int i,j;

  if (!chapeau_comparesize(ch1,ch2)) {
    fprintf(stderr,"CFACV/C) ERROR: you can not sum chapeau objects with different sizes\n");
    exit(-1);
  }
  
  //Really?
  //if (!chapeau_comparegrid(ch1,ch2)) {
  //  fprintf(stderr,"CFACV/C) ERROR: you can not sum chapeau objects with different domains\n");
  //  exit(-1);
  //}
  for (i=0;i<ch1->m;i++) {
    ch1->b[i] += ch2->b[i];
    ch1->bfull[i] += ch2->bfull[i];
    ch1->hits[i] = (ch1->hits[i]||ch2->hits[i]);
    for (j=0;j<ch1->ldad;j++) {
      ch1->A[j][i] += ch2->A[j][i];
      ch1->Afull[j][i] += ch2->Afull[j][i];
    }
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
      fwrite(&(ch->lam[i]),sizeof(double),1,ch->ofp);
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
  int i;

  ofs=fopen(filename,"w");

  fwrite(ch, sizeof(*ch), 1, ofs);

  fwrite(ch->rmin,sizeof(*(ch->rmin)),ch->dm,ofs);
  fwrite(ch->rmax,sizeof(*(ch->rmax)),ch->dm,ofs);
  fwrite(ch->N,sizeof(*(ch->N)),ch->dm,ofs);
   
  fwrite(ch->lam,sizeof(*(ch->lam)),ch->m,ofs);
  fwrite(ch->hits,sizeof(*(ch->hits)),ch->m,ofs);

  for (i=0;i<ch->ldad;i++) {
    fwrite(ch->A[i],sizeof(*(ch->A[i])),ch->m,ofs);
  }
  fwrite(ch->b,sizeof(*(ch->b)),ch->m,ofs);

  for (i=0;i<ch->ldad;i++) {
    fwrite(ch->Afull[i],sizeof(*(ch->Afull[i])),ch->m,ofs);
  }
  fwrite(ch->bfull,sizeof(*(ch->bfull)),ch->m,ofs);

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
    fread(chaux, sizeof(*chaux), 1,ofs);
    if (!chapeau_comparesize(ch,chaux)) {
      fprintf(stdout,"CFACV) ERROR, the chapeau object from the file is diferent\n");
    }
            
    // This variables are alredy set during allocation, might be better to
    // compare them instead
    ch->dm      = chaux->dm;
    ch->m    = chaux->m;
    ch->ldad    = chaux->ldad;
    
    // This variables might be not set
    ch->nupdate     = chaux->nupdate    ;
    ch->outputFreq  = chaux->outputFreq ;
    ch->outputLevel = chaux->outputLevel;
    
    // Read directly
    fread(ch->rmin,sizeof(*(ch->rmin)),ch->dm,ofs);
    fread(ch->rmax,sizeof(*(ch->rmax)),ch->dm,ofs);
    fread(ch->N,sizeof(*(ch->N)),ch->dm,ofs);
     
    for (i=0;i<ch->dm;i++) {
      ch->dr[i]  =(ch->rmax[i]-ch->rmin[i])/(ch->N[i]-1);
      ch->idr[i] =1./ch->dr[i];
    }
                     
    fread(ch->lam,sizeof(*(ch->lam)),ch->m,ofs);
    fread(ch->hits,sizeof(*(ch->hits)),ch->m,ofs);

    for (i=0;i<ch->ldad;i++) {
      fread(ch->A[i],sizeof(*(ch->A[i])),ch->m,ofs);
    }
    fread(ch->b,sizeof(*(ch->b)),ch->m,ofs);

    for (i=0;i<ch->ldad;i++) {
      fread(ch->Afull[i],sizeof(*(ch->Afull[i])),ch->m,ofs);
    }
    fread(ch->bfull,sizeof(*(ch->bfull)),ch->m,ofs);
             
    fclose(ofs);
    free(chaux);
}
    
chapeau * chapeau_allocloadstate ( char * filename ) {
    int i;
    FILE * ofs;
    chapeau * ch;
   
    fprintf(stdout,"CFACV) Allocating chapeau object from restart file\n");

    ch = (chapeau*)malloc(sizeof(chapeau));
  
    // Reading only the sizes
    ofs=fopen(filename,"r");
    fread(ch, sizeof(*ch), 1, ofs);
    fclose(ofs);
 
    // Allocating the dimension
    ch->rmin=(double*)malloc(ch->dm*sizeof(double));
    ch->rmax=(double*)malloc(ch->dm*sizeof(double));
    ch->dr=(double*)malloc(ch->dm*sizeof(double));
    ch->idr=(double*)malloc(ch->dm*sizeof(double));
    ch->N=(int*)malloc(ch->dm*sizeof(int));
    ch->r=(double*)malloc(ch->dm*sizeof(double));
    ch->f=(double*)malloc(ch->dm*sizeof(double));

    // Allocating the size. TODO: En este bloque me parece que con malloc basta.
    ch->lam=(double*)malloc(ch->m*sizeof(double));
    ch->hits=(int*)calloc(ch->m,sizeof(int));
    ch->b=(double*)calloc(ch->m,sizeof(double));
    ch->bfull=(double*)calloc(ch->m,sizeof(double));

    ch->A=(double**)calloc(ch->ldad,sizeof(double*));
    ch->Afull=(double**)calloc(ch->ldad,sizeof(double*));
    for (i=0;i<ch->ldad;i++) {
      ch->A[i]=(double*)calloc(ch->m,sizeof(double));
      ch->Afull[i]=(double*)calloc(ch->m,sizeof(double));
    }
                     
    // Now reading the state
    chapeau_loadstate(ch,filename);

    return ch;
}

void chapeau_loadlambda ( chapeau * ch, char * filename ) {
    int i;
    chapeau * chaux;
 
    if (!ch) {
      fprintf(stderr,"OTFP) ERROR in load chapeau file because holding object was not allocated\n");
      exit(-1);
    }
                                           
    chaux=chapeau_allocloadstate(filename);

    // With chaux we prevent override of addresses (hits, A, etc) when
    // reading integers (N,M,et), that is only what we want... better way to do
    // this?
    if (!chapeau_comparesize(ch,chaux)) {
      fprintf(stdout,"OTFP) ERROR, the chapeau object from the file is diferent\n");
    }
            
    // This variables are alredy set during allocation, might be better to
    // compare them instead
    ch->dm      = chaux->dm;
    ch->m    = chaux->m;
    ch->ldad    = chaux->ldad;
                                                 
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
    
    for (i=0;i<ch->m;i++) {
      ch->lam[i]=chaux->lam[i];
    }

    chapeau_free(chaux);
}

void chapeau_solve ( chapeau * ch ) {
  int i,j,J,I,s,k,l;
  double lb;
  double * Abar;
  double * bbar;
  int * pivot;
  int nred;
  int ldad;
  char aux; 

  //// parameters that are allowed to evolve lie between indices for which
  //// ch->hits[] is non-zero so extract the proper subspace
  nred=0;
  for (i=0;i<ch->m;i++)  if (ch->hits[i]) nred++;
      
  // If the nred is small, singular matrix can occur easily
  if (nred<ch->ku+1) {
   // Avoiding to print this for REOTFP simulations
   // fprintf(stderr,"Warning: No chapeau update, nred (%d) < nro diagonals (%d)\n",nred,ch->ku);
   return;
  }

  fprintf(stdout,"OTFP) Inverting matrix A (%dx%d)\n",nred,nred);
  fflush(stdout);

  // Allocating the size. TODO: En este bloque me parece que con malloc basta.
  ldad=ch->ldad+ch->ku;
  bbar=(double*)malloc(nred*sizeof(double));
  Abar=(double*)calloc(nred*ldad,sizeof(double));
  pivot=(int*)malloc(nred*sizeof(int));

  // the band storage scheme for a 
  // m = n = 6, kl = 2, ku = 2 example:
  //
  //   *    *    *    *    +    +      
  //   *    *    *    +    +    +      
  //   *    *   a02  a13  a24  a35
  //   *   a01  a12  a23  a34  a45 (uperdiagonal)
  //  a00  a11  a22  a33  a44  a55 (diagonal)
  //  a10  a21  a32  a43  a54   * 
  //  a20  a31  a42  a53   *    * 
  //
  // storate[(kl+ku+(i-j))][j]=a[i][j]; 
  // note that columns index are the same
  
  // Now, what happens if column 2 is always cero. Since 
  // matrix A is symetric row 2 is also cero. 
  //
  //   *    *    *    *    +    +      
  //   *    *    *    +    +    +      
  //   *    *    0   a13   0   a35
  //   *   a01   0    0   a34  a45 (uperdiagonal)
  //  a00  a11   0   a33  a44  a55 (diagonal)
  //  a10    0   0   a43  a54   * 
  //    0  a31   0   a53   *    * 
  //
  // The reduced matrix does not have neither the row not the column. So the
  // relabel is 1->1, 2 disapear, 3->2 and go on in all columns and rows.
  //
  //   *    *        *        *        +        +                 
  //   *    *        *        +        +        +                 
  //   *    *        0 ( - )  0 ( 0 ) A24(a35)  -                 
  //   *   a01      A12(a13) A23(a34) A34(a45)  -  (uperdiagonal) 
  //  a00  a11      A22(a33) A33(a44) A44(a55)  -  (diagonal)     
  //  a10  A21(a31) A32(a43) A43(a54) A54(a65)  *                 
  //    0   *       A42(a53) A53(a64)  * (a75)  *                 
  //
  //   The net flux is:
  //
  //   *    *    *    *    +    +      
  //   *    *    *    +    +    +      
  //   *    *        L         <-
  //   *   a01            <-   <-
  //  a00  a11       <-   <-   <-
  //  a10            <-   <-   <-
  //        ^        <-   <-   <-
  //                       
  // So it is better to loop in the other space!
  
  // Add the new and old statistic to the reduced matrix
  J=0;
  for (j=0;j<ch->m;j++) {

    if (!ch->hits[j]) continue;
    bbar[J]=-(ch->b[j]+ch->bfull[j]); //TODO: Trace back the origin of the minus. See accumulate procedures. 
    
    // Upper diagonals (and main diagonal) in the row
    I=J;
    for (i=j;i<=j+ch->ku;i++) {
      if (i>=ch->m) break;
      if (!ch->hits[i]) continue;

      Abar[2*ch->ku+I-J+ldad*J]=ch->A[ch->ku+i-j][j]+ch->Afull[ch->ku+i-j][j];
      I++;
    }  

    // Lower diagonals in the row
    I=J-1;
    for (k=1;k<=ch->ku;k++) {
      i=j-k;
      if (i<0) break;
      if (!ch->hits[i]) continue;

      Abar[2*ch->ku+I-J+ldad*J]=ch->A[ch->ku+i-j][j]+ch->Afull[ch->ku+i-j][j];
      I--;
    }
    J++;
  } 

  if (nred!=J) {fprintf(stderr,"Bad matrix nred size: %d != %d\n",nred,J);exit(-1);}

  //// TODO: To call a Fortran routine from C we have to transpose the matrix.
  //// However, this is is no needed if Abar is symmetric, but I should change
  //// Abar to 1 dimension array in the rest of the code
  //for (i=0; i<nred; i++){
  //  for(j=0; j<ldad; j++) AT[j+nred*i]=Abar[j][i];           
  //}                                               

  // find solution using LAPACK routine SGESV.
  J=1;                       
  aux='N';

  // LU factorization of A.
  dgbtrf_(&nred, &nred, &ch->ku, &ch->ku, Abar, &ldad, pivot, &s); 
              
  if (s!=0) {
    fprintf(stdout,"Matrix nearly singular with flag: %d\n",s);

    free(Abar);
    free(bbar);
    free(pivot);   
    
    // chapeau_solve_secure(ch);
    s=s-1;
    I=0;
    for (i=0;i<ch->m;i++) {
      if (!ch->hits[i]) continue;
      if (I==s) {
        fprintf(stdout,"Removing row %d\n",I);
        ch->hits[i]=0;
        break;
      }
      I++;
    }    
    fflush(stdout);

    // chapeau_solve_secure(ch);
    chapeau_solve(ch);

    return;
  }
 
  // Find solution
  dgbtrs_(&aux,&nred, &ch->ku, &ch->ku, &J, Abar, &ldad, pivot, bbar, &nred, &s);
  
  // update the vector of coefficients
  I=0;
  for (i=0;i<ch->m;i++) {
    if (!ch->hits[i]) continue;

    // Insted of solving Ax=b, I rather solve (dr[0]*dr[0]*A)x/dr[0]=(dr[0]*b).
    // So, I have to remember multiply the solution by dr[0].
    lb=bbar[I]*ch->dr[0];
    if (lb!=lb) {fprintf(stderr,"PARANOIA) Tripped at I=%d i=%d: %.5f?\n",I,i,lb);exit(-1);}
    ch->lam[i]=lb;
    
    //lo=gsl_vector_get(ch->lam,i);
    //gsl_vector_set(ch->lam,i,alpha*lo+(1-alpha)*lb);

    I++;
  } 
   
  free(Abar);
  free(bbar);
  free(pivot);   

  //chapeau_baselinehits(ch); 
}

void chapeau_solve_secure ( chapeau * ch ) {
  int i,j,J,I,k,s;
  double lb;
  double * Abar;
  double * Ebar;
  double * Dbar;
  double * nullbar;
  double * U;
  double * VT;
  double * bbar;
  double * work;
  int nred;
  char aux; 

  //DB//for (i=0;i<ch->m;i++) fprintf(stderr,"AAA %i %i\n",i,ch->hits[i]); exit(1);

  //// parameters that are allowed to evolve lie between indices for which
  //// ch->hits[] is non-zero so extract the proper subspace
  nred=0;
  for (i=0;i<ch->m;i++)  if (ch->hits[i]) nred++;
      
  // If the nred is small, singular matrix can occur
  if (nred<ch->ku+1) {
   fprintf(stderr,"Warning: No chapeau update: 0 columns (%d) > diagonals (%d)\n",nred,ch->ku+1);
   return;
  }

  // Allocating the size. TODO: En este bloque me parece que con malloc basta.
  bbar=(double*)malloc(nred*sizeof(double));
  work=(double*)malloc(4*nred*sizeof(double));
  Abar=(double*)calloc(nred*ch->ldad,sizeof(double));
  Dbar=(double*)calloc(nred,sizeof(double));
  Ebar=(double*)calloc(nred-1,sizeof(double));
  U=(double*)calloc(nred*nred,sizeof(double));
  VT=(double*)calloc(nred*nred,sizeof(double));

  // the band storage scheme for a 
  // m = n = 6, kl = 2, ku = 2 example:
  //
  //   *    *   a02  a13  a24  a35
  //   *   a01  a12  a23  a34  a45 (uperdiagonal)
  //  a00  a11  a22  a33  a44  a55 (diagonal)
  //  a10  a21  a32  a43  a54   * 
  //  a20  a31  a42  a53   *    * 
  //
  // storate[(ku+(i-j))][j]=a[i][j]; 
  // note that columns index are the same
  
  // Now, what happens if column 2 is always cero. Since 
  // matrix A is symetric row 2 is also cero. 
  //
  //   *    *    0   a13   0   a35
  //   *   a01   0    0   a34  a45 (uperdiagonal)
  //  a00  a11   0   a33  a44  a55 (diagonal)
  //  a10    0   0   a43  a54   * 
  //    0  a31   0   a53   *    * 
  //
  // The reduced matrix does not have neither the row not the column. So the
  // relabel is 1->1, 2 disapear, 3->2 and go on in all columns and rows.
  //
  //   *    *        0 ( - )  0 ( 0 ) A24(a35)  -                 
  //   *   a01      A12(a13) A23(a34) A34(a45)  -  (uperdiagonal) 
  //  a00  a11      A22(a33) A33(a44) A44(a55)  -  (diagonal)     
  //  a10  A21(a31) A32(a43) A43(a54) A54(a65)  *                 
  //    0   *       A42(a53) A53(a64)  * (a75)  *                 
  //
  //   The net flux is:
  //
  //   *    *        L         <-
  //   *   a01            <-   <-
  //  a00  a11       <-   <-   <-
  //  a10            <-   <-   <-
  //        ^        <-   <-   <-
  //                       
  // So it is better to loop in the other space!
  
  // Add the new and old statistic to the reduced matrix
  J=0;
  for (j=0;j<ch->m;j++) {

    if (!ch->hits[j]) continue;
    bbar[J]=-(ch->b[j]+ch->bfull[j]); //TODO: Trace back the origin of the minus. See accumulate procedures. 
    
    // Upper diagonals (and main diagonal) in the row
    I=J;
    for (i=j;i<=j+ch->ku;i++) {
      if (i>=ch->m) break;
      if (!ch->hits[i]) continue;

      Abar[ch->ku+I-J+ch->ldad*J]=ch->A[ch->ku+i-j][j]+ch->Afull[ch->ku+i-j][j];
      I++;
    }  

    // Lower diagonals in the row
    I=J-1;
    for (k=1;k<=ch->ku;k++) {
      i=j-k;
      if (i<0) break;
      if (!ch->hits[i]) continue;

      Abar[ch->ku+I-J+ch->ldad*J]=ch->A[ch->ku+i-j][j]+ch->Afull[ch->ku+i-j][j];
      I--;
    }
    J++;
  } 
     
  //if (J!=I) {fprintf(stderr,"Matrix not square!!");exit(-1);}
  if (nred!=J) {fprintf(stderr,"Bad matrix nred size: %d != %d\n",nred,J);exit(-1);}

  //// TODO: To call a Fortran routine from C we have to transpose the matrix.
  //// However, this is is no needed if Abar is symmetric, but I should change
  //// Abar to 1 dimension array in the rest of the code
  //for (i=0; i<nred; i++){
  //  for(j=0; j<ch->ldad; j++) AT[j+nred*i]=Abar[j][i];           
  //}                                               

  // Reduction to bidiagonal form by orthogonal transformation A=U B VT.
  I=0;                       
  J=1;                       
  aux='B';
  dgbbrd_(&aux, &nred, &nred, &I, &ch->ku, &ch->ku, Abar, &ch->ldad, Dbar, Ebar, 
          U, &nred, VT, &nred, nullbar, &J, work, &s); 
              
  if (s!=0) {fprintf(stderr,"dgbbrd: %d argument had illegal value\n",s);exit(-1);}

  // SVD of the bidiagonal form
  aux='U';
  dbdsqr_(&aux,&nred, &nred, &nred, &I, Dbar, Ebar, VT, &nred, U, &nred, nullbar, &J, work, &s);
  if (s!=0) {fprintf(stderr,"dbdsqr: error with flag %d\n",s);exit(-1);}


  // update the vector of coefficients
  I=0;
  for (i=0;i<ch->m;i++) {

    if (!ch->hits[i]) continue;
    lb=0;

    // TODO: exchange i and j loops

    for (j=0;j<nred;j++) {
      // Filtro valor singular cercano a cero
      if (Dbar[j]<1.e-10) {
        lb+=0.;
      } else {

        // El elemento Ik de V D^I U^T es
        //   V[i][j]*1./Dbar[j]*UT[j][k] 

        for (k=0;k<nred;k++)  lb+=VT[j+nred*I]*U[k+nred*j]*bbar[k]/Dbar[j];
        // for (k=0;k<nred;k++)  lb+=VT[I+nred*j]*U[j+nred*k]*bbar[k]/Dbar[j];

      }
    }

    if (lb!=lb) {fprintf(stderr,"PARANOIA) Tripped at I=%d i=%d: %.5f?\n",I,i,lb);exit(-1);}
    ch->lam[i]=lb*ch->dr[0];
    I++;

  } 

  free(Abar);
  free(Dbar);
  free(Ebar);
  free(VT);
  free(U);
  free(bbar);
  free(work);   

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
    fprintf(stdout,"INFO) %i,%.5f\n",i++,knots);
    ch->lam[i]=knots;
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

  la=ch->lam[m];
  lb=ch->lam[m+1];
  return la+(lb-la)*dm; 

} 
    
double chapeau_evalf_2simplex ( chapeau * ch, double z1, double z2 ) {
  int i,j,ni,nj,nk;
  double f,dx,dy;

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

    f-=(dx*ch->idr[0]+1)*ch->lam[nk];
    f-=(dy*ch->idr[1]+1)*ch->lam[nk];

  } else {
    nk=(j+1)*(ch->N[0]-1)+(i+1);
    /*nj---------nk
       |-\        |  
       |  -\   x  |  
       |    -\    |  
       |      -\  |  
       |        -\|  
       *---------ni */

    f+=(dx*ch->idr[0]+1)*ch->lam[nk];
    f+=(dy*ch->idr[1]+1)*ch->lam[nk];

  }
    
  // Sumando contribucion de ni y nj
  f+=(dx*ch->idr[0]+1)*ch->lam[ni];
  f+=(dy*ch->idr[1]+1)*ch->lam[nj];

  return f;
}  

char * chapeau_serialize ( chapeau * ch ) {
  // return a str containing the serialized chapeau with partial statistics. If
  // this is the central replica, or non replica scheme is used, this serialize
  // the full statistic.
  int i, j;
  int size;
  char* buf;

  // size = (3*ch->m-1)*13+ch->m+1;
  size  = ch->m*13;           // space for ch->b
  size += ch->m*ch->ldad*13;; // space for ch->A
  size += ch->m;;             // space for ch->hits
  size += 1;;                 // null character?
  buf  = (char*) malloc(size*sizeof(char));

  // Seems that if I do this:
  //  return (char*)&ch
  //  is not portable, since the way that the cast is made is undefined
  //  (depends of the machine). So it is needed a proper serialization
  
  size=0;
  for (i=0;i<ch->m;i++) {

    // Writing the partial b
    sprintf(buf+size, "%13.5e", ch->b[i]); size+=13;

    // Writing hits
    sprintf(buf+size, "%1i", ch->hits[i]); size+=1;
  }

  // TODO take advantage of the simetric property to reduce half of the numbers
  // Writing the partial A (banded)
  for (i=0;i<ch->ldad;i++) {
    for (j=0;j<ch->m;j++) {
      sprintf(buf+size, "%13.5e",ch->A[i][j]); size+=13;
    }
  }

  return buf; 
}     

void chapeau_addserialized ( chapeau *ch, char * str ) {
  // str contains the serialized chapeau that comes from the partial sampling
  // of other replica partial statistic information. This is added to the
  // partial sampling ot the principal replica computing the sum.
  int i,j,err;
  int size1;
  int size2;
  double aux;
  int iaux;
  char word1[14],word2[2];
  
  // Add null terminators
  word1[13]='\0';
  word2[1]='\0';

  size1  = ch->m*13;           // space for ch->b
  size1 += ch->m*ch->ldad*13;; // space for ch->A
  size1 += ch->m;;             // space for ch->hits
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
    ch->b[i]+=aux;

    // Reading hits
    memcpy(word2, &str[size1], 1 ); size1+=1;
    err=sscanf(word2,"%1i",&iaux); if(!err) {fprintf(stderr,"CFACV/C) Error 1111 on read %s\n",word2);}
    ch->hits[i]=(ch->hits[i]||iaux);
  }
 
  // Reading the partial A (simetric and tridiagonal)
  for (i=0;i<ch->ldad;i++) {
    for (j=0;j<ch->m;j++) {
      memcpy(word1, &str[size1], 13 ); size1+=13;
      err=sscanf(word1,"%13le",&aux); if(!err) {fprintf(stderr,"CFACV/C) Error 1112 on read %s\n",word1);}
      ch->A[i][j]+=aux;
    }
  }

  //TODO, give error if rmin and rmax are different
  
}

void chapeau_setserialized ( chapeau *ch, char * str ) {
  // str contains the full statistics information that comes from all the
  // replicas contributions. Therefore, this should be stored in Afull and
  // bfull and the partial A and b should be reset to cero.
  int i,j,err;
  int size1;
  int size2;
  double aux;
  int iaux;
  char word1[14],word2[2];
  
  // Add null terminators
  word1[13]='\0';
  word2[1]='\0';
  
  size1  = ch->m*13;           // space for ch->b
  size1 += ch->m*ch->ldad*13;; // space for ch->A
  size1 += ch->m;;             // space for ch->hits
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
    ch->bfull[i]=aux;

    // Read hits
    memcpy(word2, &str[size1], 1 ); size1+=1;
    err=sscanf(word2,"%1i",&iaux); if(!err) {fprintf(stderr,"CFACV/C) Error 2112 on read %s\n",word1);}
    ch->hits[i]=iaux;
  }
     
  // Reading the partial A (simetric and tridiagonal)
  for (i=0;i<ch->ldad;i++) {
    for (j=0;j<ch->m;j++) {
      memcpy(word1, &str[size1], 13 ); size1+=13;
      err=sscanf(word1,"%13le",&aux); if(!err) {fprintf(stderr,"CFACV/C) Error 1112 on read %s\n",word1);}
      ch->Afull[i][j]+=aux;
    }
  }

  // Reset the variables to hold the partial statistic information
  for (j=0;j<ch->m;j++) ch->b[j]=0.;

  for (i=0;i<ch->ldad;i++) {
    for (j=0;j<ch->m;j++) ch->A[i][j]=0.;
  }

  //TODO, give error if rmin and rmax are different
}

 int accumulate_1D( chapeau * ch ) { 
  int m;
  
  // Following the paper E definition A and b should be
  //
  // A=1./nsteps \sum_{steps} \sum_{i} dphi_m(z)/dz dphi_n(z)/dz
  // b=-1./nsteps \sum_{steps} \sum_{i} dphi_m(z)/dz F
  //
  // But here we compute -b*nsteps and A*nsteps and then A and b are scaled
  // by the proper factor when equation A\lambda=b is solved (see the
  // chapeau_solve procedure)
  
  // Early return to avoid interpolations out of the domain
  if ( ch->r[0] > ch->rmax[0] ) return 0;
  if ( ch->r[0] <= ch->rmin[0] ) return 0;
 
  // Insted of solve Ax=b, I rather solve (dr[0]*dr[0]*A)x/dr[0]=(dr[0]*b).
  // Then, I have to remember multiply the solution by dr[0]. The factor idr[0]
  // in the terms of A and b are now factor 1. 
                       
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
  ch->b[m+1]+=ch->f[0];
  ch->A[ch->ku][m+1]+=1.0;

  /*    |\   -1/dz
     n=m| \
        |  \
    o---o-x-o---o--*/ 
  ch->b[m]-=ch->f[0];
  ch->A[ch->ku][m]+=1.0;

  /*    |\ /| 
       m| X |n=m+1
        |/ \| 
    o---o-x-o---o--*/ 
  ch->A[ch->ku+1][m]-=1.0;
  ch->A[ch->ku-1][m+1]-=1.0;

  return 0;
}
 
int accumulate_2D( chapeau * ch ) { 
  int i,j,ni,nj,nk;
  double dx,dy,ratio,ratio2;
  
  // Following the paper E definition A and b should be
  //
  // A=1./nsteps \sum_{steps} \sum_{i} dphi_m(z)/dz dphi_n(z)/dz
  // b=-1./nsteps \sum_{steps} \sum_{i} dphi_m(z)/dz F
  //
  // But here we compute -b*nsteps and A*nsteps and then A and b are scaled
  // by the proper factor when equation A\lambda=b is solved (see the
  // chapeau_solve procedure)
 
  // Early return to avoid interpolations out of the domain
  if ( ch->r[0] > ch->rmax[0] ) return 0;
  if ( ch->r[1] > ch->rmax[1] ) return 0;
  if ( ch->r[0] <= ch->rmin[0] ) return 0;
  if ( ch->r[1] <= ch->rmin[1] ) return 0;

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

  // The ratio allow me to transform the matrix A and vector b to make them
  // indpendent of the space. So, insted of solve Ax=b, I rather solve
  // (dr[0]*dr[0]*A)x/dr[0]=(dr[0]*b). Then, I have to remember multiply the
  // solution by dr[0]. The factor idr[0] in the terms of A and b are now
  // factor 1 and the terms with idr[1] are now the factor ratio. 
  ratio=ch->dr[0]*ch->idr[1];
  ratio2=ratio*ratio;

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


  // The band storage scheme for a 
  // M = N = 6, KL = 2, KU = 1 example:
  //
  //   *    *    *    +    +    +      
  //   *    *    +    +    +    + 
  //   *   a01  a12  a23  a34  a45 (uperdiagonal)
  //  a00  a11  a22  a33  a44  a55 (diagonal)
  //  a10  a21  a32  a43  a54   * 
  //  a20  a31  a42  a53   *    * 
  //
  // storate[(kl+ku+(i-j))][j]=a[i][j]; 
  //
  // Note that columns index are the same
 
           
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
    
    ch->b[ni]     += ch->f[0];
    ch->A[ch->ku][ni]+= 1.0;         

    ch->b[nj]     += ratio*ch->f[1];
    ch->A[ch->ku][nj] += ratio2;     
     
    ch->b[nk]     -= ch->f[0];
    ch->b[nk]     -= ratio*ch->f[1];
    ch->A[ch->ku][nk]+= 1.0;         
    ch->A[ch->ku][nk]+= ratio2;      
    
    // Aca podria evitar 2 operaciones si eligiera el triangulo superior
    ch->A[ch->ku+(ni-nk)][nk]-= 1.0;   
    ch->A[ch->ku+(nk-ni)][ni]-= 1.0;   
    ch->A[ch->ku+(nj-nk)][nk]-= ratio2;
    ch->A[ch->ku+(nk-nj)][nj]-= ratio2;

    ch->hits[ni]=1;
    ch->hits[nj]=1;
    ch->hits[nk]=1;
 
    if(ch->periodic){
      if(i==ch->N[0]){
        ni=ni-ch->N[0];
        ch->b[ni]     += ch->f[0];
        ch->A[ch->ku][ni]+= 1.0;
        ch->hits[ni]=1;
      }
      if(j==ch->N[1]){
        nj=i;
        ch->b[nj]     += ratio*ch->f[1];
        ch->A[ch->ku][nj] += ratio2;     
        ch->hits[nj]=1;
      }
    }

  } else {
    //nk=(j+1)*ch->N[0]+ch->ku+1;
    nk=nj+1;
 
    /*nj---------nk
       |-\        |  
       |  -\   x  |  
       |    -\    |  
       |      -\  |  
       |        -\|  
       *---------ni */

    ch->b[ni]        -= ratio*ch->f[1];
    ch->A[ch->ku][ni] += ratio2;
                     
    ch->b[nj]        -= ch->f[0];
    ch->A[ch->ku][nj] += 1.0;
                     
    ch->b[nk]        += ch->f[0];
    ch->b[nk]        += ratio*ch->f[1];
    ch->A[ch->ku][nk] += 1.0;
    ch->A[ch->ku][nk] += ratio2;

    // Aca podria evitar 2 operaciones si eligiera el triangulo superior
    ch->A[ch->ku+(nk-nj)][nj] -= 1.0;
    ch->A[ch->ku+(nj-nk)][nk] -= 1.0;
    ch->A[ch->ku+(nk-ni)][ni] -= ratio2;
    ch->A[ch->ku+(ni-nk)][nk] -= ratio2;
 
    ch->hits[ni]=1;
    ch->hits[nj]=1;
    ch->hits[nk]=1;

    if(ch->periodic){
      if(i==ch->N[0]){
        ni=ni-ch->N[0];
        ch->b[ni]        -= ratio*ch->f[1];
        ch->A[ch->ku][ni] += ratio2;
        nk=nk-ch->N[0];
        ch->b[nk]        += ch->f[0];
        ch->b[nk]        += ratio*ch->f[1];
        ch->A[ch->ku][nk] += 1.0;
        ch->A[ch->ku][nk] += ratio2;
        ch->hits[ni]=1;
        ch->hits[nk]=1;
      }
      if(j==ch->N[1]){
        nj=i;
        ch->b[nj]        -= ch->f[0];
        ch->A[ch->ku][nj] += 1.0;
        nk=i+1;
        ch->b[nk]        += ch->f[0];
        ch->b[nk]        += ratio*ch->f[1];
        ch->A[ch->ku][nk] += 1.0;
        ch->A[ch->ku][nk] += ratio2;
        ch->hits[nj]=1;
        ch->hits[nk]=1;
      }
      if (nk==ch->N[0]){
        nk=0;
        ch->b[nk]        += ch->f[0];
        ch->b[nk]        += ratio*ch->f[1];
        ch->A[ch->ku][nk] += 1.0;
        ch->A[ch->ku][nk] += ratio2;
        ch->hits[nk]=1;
      }
    }
     
  }
 
  return 0;
}
 

