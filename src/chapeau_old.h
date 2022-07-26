typedef struct OLD_CHAPEAU {
  // Conjunto de funciones base, sets de chapeaus que
  // constituyen una base para alguna funcion reconstruida por elementos
  // finitos. Los coeficientes tambien se almacenana aca.

  // number of data (not timsteps) to acumulate before update
  int nupdate;
 
  // Dimension
  int dm; 

  // Discretizacion del dominio por dimension
  int * N;       // size in mult-dimensional real space 

  int m;       // m[0] size in 1-dimensional wrap space 

  int ku;      // number of uperdiagonals 
  int ldad;    // dimension of the packed matrix
                 
               
  int periodic; 
  double * rmin;
  double * rmax;
  double * dr;
  double * idr;  

  double * r; // Vector para guardar la posicion actual
  double * f; // Vector para guardar el gradiente de la funcion

  // Para sacar los coeficientes
  FILE * ofp; 
  int outputFreq;
  int outputLevel;

  // Since restart file is open and closed in the subrroutine, might be better
  // to store the name and not the unit? 
  char filename[255]; 

  // Variables that accumulate partial statistics to optimze FEP coeficients.
  // This is used to send between replicas and is a private copy of the
  // information acquired for the self sampling. After all the replicas
  // comuncates, and if this replica is not the center replica that add all the
  // other contributions, this variables are set to cero. If this is the
  // central replica, or non replica scheme is used, this variables contains
  // the full statistics of the sampling.
  double * b;
  double * lam;
  int * hits; 

  //matrix A will have dimensions (dm-1)*3*(ch->N[0]-1)+1 x m
  double ** A;

  // If this replica is not the center replica that add all the other
  // contributions, this variables accumulates the full statistics (including
  // all other replicas) infromation needed to optimze FEP coeficients. It this
  // replica is the center replica, ot non replica scheme is used, this
  // variables are always empty.
  double * bfull;
  double ** Afull;

  //pointer to procedure
  int (*accumulate)(struct OLD_CHAPEAU * self);
                
} old_chapeau;            
 
void old_chapeau_loadstate ( old_chapeau * ch, char * filename ) {
    old_chapeau * chaux;
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
    chaux = (old_chapeau*)malloc(sizeof(old_chapeau));
    fread(chaux, sizeof(*chaux), 1,ofs);
            
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
 
old_chapeau * old_chapeau_allocloadstate ( char * filename ) {
    int i;
    FILE * ofs;
    old_chapeau * ch;
   
    fprintf(stdout,"CFACV) Allocating chapeau object from restart file\n");

    ch = (old_chapeau*)malloc(sizeof(old_chapeau));
  
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
    old_chapeau_loadstate(ch,filename);

    return ch;
}


