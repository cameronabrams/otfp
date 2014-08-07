#ifndef CHAPEAU_H
#define CHAPEAU_H
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

  FILE * ofs; // for save current chapeu status
  // file name is the same of ofp but with a .restart prefix

  // single-particle-sums; initialize at every step, every particle
  double *** s;  // [particle][dimension][peak]
  
  // This mask allows to control each interacion.
  // each mask item have the value 0,-1,1 or 2
  // Interaction is done when mask[i]+mask[j]=0
  // Here int because unsigned is not safe with swig
  int * mask;  // [particle]

  gsl_vector * b;
  gsl_matrix * A;
  gsl_vector * lam;    // vector of coefficients -- these are what OTFP optimizes!
  long * hits;  // number of hits in each bin value of r

  double alpha;
  int updateinterval;

} chapeau;

chapeau * chapeau_alloc ( int m, double rmin, double rmax, int npart );

void chapeau_setUpdateInterval ( chapeau * ch, int i );

void chapeau_setPeaks ( chapeau * ch, double * peaks );

void chapeau_pair_eval_g ( chapeau * ch, double z, double * u, double * g_r );

void chapeau_setupoutput ( chapeau * ch, char * filename, int outputFreq, int outputLevel );
void chapeau_output ( chapeau * ch, int timestep );
void chapeau_savestate ( chapeau * ch, int timestep );
chapeau * chapeau_allocloadstate ( char * filename ) ;

void chapeau_init_global_accumulators ( chapeau * ch );

void chapeau_init_particle_sums ( chapeau * ch );

void chapeau_increment_particle_sum ( chapeau * ch, int i, int j, double * Zij, double zij );
void chapeau_increment_global_accumulators ( chapeau * ch, int i, double * F );
void chapeau_update_peaks ( chapeau * ch, int nsamples, int timestep );

#endif
