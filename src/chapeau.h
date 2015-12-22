#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_linalg.h>

// Chapeau type definition:  This defines the pair potential as 
// a linear expansion of chapeau functions ("peaks").  It also
// contains the variables responsible for evolving the vector
// of coefficients \lambda; these are the single-particle sums
// and the global accumulators
typedef struct CHAPEAU {
  int m; // number of peaks
  int N; // number of particles

  double rmin, rmax, dr, idr;  // range and increment or argument of
			       // pair potential

  FILE * ofp; // for output chapeu
  int outputFreq;
  int outputLevel;

  // Since restart file is open and closed in the subrroutine, might be better
  // to store the name and not the unit? Besides, if I want to take advantage
  // of save state subroutines using another name?
  char filename[255]; 

  // file name is the same of ofp but with a .restart prefix

  //// single-particle-sums; initialize at every step, every particle
  //double *** s;  // [particle][dimension][peak]
  
  // This mask allows to control each interacion.
  // each mask item have the value 0,-1,1 or 2
  // Interaction is done when mask[i]+mask[j]=0
  // Here int because unsigned is not safe with swig
  int * mask;  // [particle]


  // Variables that accumulate partial statistics to optimze FEP coeficients.
  // This is used to send between replicas and is a private copy of the
  // information acquired for the self sampling. After all the replicas
  // comuncates, and if this replica is not the center replica that add all the
  // other contributions, this variables are set to cero. If this is the
  // central replica, or non replica scheme is used, this variables contains
  // the full statistics of the sampling.
  gsl_vector * b;
  gsl_matrix * A;
  gsl_vector * lam;   
  int * hits; 

  // If this replica is not the center replica that add all the other
  // contributions, this variables accumulates the full statistics (including
  // all other replicas) infromation needed to optimze FEP coeficients. It this
  // replica is the center replica, ot non replica scheme is used, this
  // variables are always empty.
  gsl_vector * bfull;
  gsl_matrix * Afull;

  double alpha;

  // number of data (not timsteps) acumulated
  int nsample;

  // number of data (not timsteps) to acumulate before update
  int nupdate;

} chapeau;

chapeau * chapeau_alloc ( int m, double rmin, double rmax, int npart );

void chapeau_setUpdateInterval ( chapeau * ch, int i );

void chapeau_pair_eval_g ( chapeau * ch, double z, double * u, double * g_r );

// Output system
void chapeau_setupoutput ( chapeau * ch, char * filename, int outputFreq, int outputLevel );
void chapeau_output ( chapeau * ch, int timestep );

// Restart system
void chapeau_savestate ( chapeau * ch, int timestep, char * filename );
chapeau * chapeau_allocloadstate ( char * filename ) ;
void chapeau_loadstate ( chapeau * ch, char * filename ) ;

void chapeau_init_global_accumulators ( chapeau * ch );

void chapeau_init_particle_sums ( chapeau * ch );

void chapeau_increment_particle_sum ( chapeau * ch, int i, int j, double * Zij, double zij );
void chapeau_increment_global_accumulators ( chapeau * ch, int i, double * F );
void chapeau_update_peaks ( chapeau * ch );
void chapeau_sum ( chapeau * ch1, chapeau * ch2 );

double chapeau_evalf ( chapeau * ch, double z );
char * chapeau_serialize ( chapeau * ch );
void chapeau_addserialized ( chapeau * ch, char * str );
void chapeau_setserialized ( chapeau *ch, char * str );

void chapeau_set_peaks ( chapeau * ch, char * filename );
