#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Chapeau type definition:  This defines the pair potential as 
// a linear expansion of chapeau functions ("peaks").  It also
// contains the variables responsible for evolving the vector
// of coefficients \lambda; these are the single-particle sums
// and the global accumulators
typedef struct CHAPEAU {
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
                 
               
  int * periodic; 
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
  int (*accumulate)(struct CHAPEAU * self);
                
} chapeau;

chapeau * chapeau_alloc ( int dm, double * rmin, double * rmax, int * N, int * periodic );
void chapeau_free ( chapeau * ch );
int chapeau_comparesize ( chapeau * ch1,  chapeau * ch2);
int chapeau_comparegrid ( chapeau * ch1,  chapeau * ch2);

void chapeau_sum ( chapeau * ch1, chapeau * ch2 );


void chapeau_setUpdateInterval ( chapeau * ch, int i );

void chapeau_pair_eval_g ( chapeau * ch, double z, double * u, double * g_r );

// Output system
void chapeau_setupoutput ( chapeau * ch,  char * outfile, char * restartfile, int outputFreq, int outputLevel );
void chapeau_output ( chapeau * ch, int timestep );

// Restart system
void chapeau_savestate ( chapeau * ch, char * filename );
chapeau * chapeau_allocloadstate ( char * filename );
void chapeau_loadstate ( chapeau * ch, char * filename );
void chapeau_loadlambda ( chapeau * ch, char * filename );

void chapeau_init_global_accumulators ( chapeau * ch );

void chapeau_init_particle_sums ( chapeau * ch );

void chapeau_increment_particle_sum ( chapeau * ch, int i, int j, double * Zij, double zij );
void chapeau_increment_global_accumulators ( chapeau * ch, int i, double * F );
void chapeau_solve ( chapeau * ch );
void chapeau_solve_secure ( chapeau * ch );

double chapeau_evalf_1simplex ( chapeau * ch, double z );
char * chapeau_serialize ( chapeau * ch );
void chapeau_addserialized ( chapeau * ch, char * str );
void chapeau_setserialized ( chapeau *ch, char * str );

void chapeau_set_peaks ( chapeau * ch, char * filename );
//void chapeau_baselinehits ( chapeau * ch );
//void chapeau_setmref ( chapeau * ch, double z );


int accumulate_1D ( chapeau * ch );
int accumulate_2D ( chapeau * ch );
