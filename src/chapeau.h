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

  gsl_vector * lam;    // vector of coefficients -- these are what OTFP optimizes!
  gsl_vector * lambar; // 
  gsl_vector * newlam; //
 
  int * hits;  // number of hits in each bin value of r

  FILE * ofp;
  int outputFreq;
  int outputLevel;

  // single-particle-sums; initialize at every step, every particle
  double *** s;  // [particle][dimension][peak]

  // global accumulators
  gsl_vector * b, * bbar;
  gsl_matrix * A, * Abar;

  gsl_permutation * Permutation;

  double alpha;
  int updateinterval;

} chapeau;

chapeau * chapeau_alloc ( int m, double rmin, double rmax, int npart );

void chapeau_setUpdateInterval ( chapeau * ch, int i );

void chapeau_setPeaks ( chapeau * ch, double * peaks );

void chapeau_pair_eval_g ( chapeau * ch, double z, double * u, double * g_r );

void chapeau_setupoutput ( chapeau * ch, char * filename, int outputFreq, int outputLevel );
void chapeau_output ( chapeau * ch, int timestep );

void chapeau_init_global_accumulators ( chapeau * ch );

void chapeau_init_particle_sums ( chapeau * ch );

void chapeau_increment_particle_sum ( chapeau * ch, int i, int j, double * Zij, double zij );
void chapeau_increment_global_accumulators ( chapeau * ch, int i, double * F );
void chapeau_update_peaks ( chapeau * ch, int nsamples, int timestep );

#endif
