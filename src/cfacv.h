#pragma once

/* Collective Variables module for NAMD via tclForces
 * Cameron F. Abrams
 *
 * Adapted for On-the-fly free-energy parameterization via TAMD
 * Sergio A. Paz
 *
 * This header declares structured data types and functions for
 * handling collective variables and their restraints. 
 * 
 * The central concept is that of a "restraint" which is applied to
 * one or a linear combination of "collective variables".  A restraint
 * computes a force on a collective variable ($\theta$) based on its
 * value and the value of a "target" or "auxiliary" variable ($z$).
 * The auxiliary variable can be constant, or driven according to
 * steered MD or temperature-accelerated MD.  It is the job of the
 * cfacv module to compute these forces after receiving from the NAMD
 * tclforces "calcforces" function the positions of all atoms or
 * atom-groups that participate in each CV owned by each restraint.
 * It is the job of the tclforces module to receive those forces and
 * transmit them to the atoms. This two-way communication happens
 * once per timestep on the master node during NAMD execution.
 *
 * In this particular OTFP impelementation, we are using linear chapeau
 * functions as the basis functions.  We use this to compute effective pair
 * potentials between centers of mass of atomgroups (like water molecules).
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "centers.h"
//#include "wrapcoords.h"
#include "chapeau.h"       // basis functions

#include "cvs.h"           // collective variables
#include "measurements.h"  // random numbers


// Temperature-accelerated MD options structure; will be owned by the
// restraint structure
typedef struct TamdOptStruct {
  double kT;           // fictitious thermal energy (kcal/mol)
  double gamma;        // fictitious friction (1/ps)
  double dt;           // timestep (ps)
  double ginv;         // inverse friction
  double noise;        // noise
} tamdOptStruct;


// Steered MD options structure; will be owned by the restraint
// structure
typedef double (*smdUpdateFunc) (double, double, int );
typedef struct SMDOPTSTRUCT {
  int t0;              // timestep to begin steering
  int t1;              // timestep to end steering (and hold restraint fixed)
  double invinterval;  // 1/(t1-t0)
  double initval;      // initial value of steered variable
  double increment;    // increment of steered variable per timestep
  double target;       // target value of steered variable
  smdUpdateFunc update; // function used to update steered variable
} smdOptStruct;


// The Restraint structure
enum {HARMONIC, HARMCUTO, NULL_RF};

typedef struct restraint {
  double k;               // spring constant
  double z;               // target value
  double val;             // restraint value

  double f;              // force on restraint
  double u;              // potential energy stored by restraint
  int nCV;               // number of collective variables in SYSTEM
  double * cvc;          // the coefficient each collective variable
                         // contributes to the restraint.  This allows
                         // a restraint to be applied on any linear
                         // combination of collective variables.  Most
                         // commonly, however, this will be an array
                         // with only one non-zero value (of 1) that
                         // indicates which CV this restraint is
                         // applied to.

  tamdOptStruct * tamdOpt; // pointer to the tamd options structure
  double tamd_noise;
  double tamd_restraint;

  int evolve;        // Indicate when z must evolve. When tamd_evolve goes
                     // from 0 to 1 z value is initialized

  smdOptStruct * smdOpt;   // pointer to the smd options structure

  int rfityp;                   // type of the restraining function
                                // (Harmonic or Periodic)

  
  // pointer to the evolution type 
  int (*evolveFunc)(struct restraint * self, double f);
   
  // pointer to the restraining energy and force function
  //restrEnergyFunc energyFunc;   
  int (*energyFunc)(struct restraint * self);

  // boundaries
  double min;        
  double max;
  double half_domain;
  double boundk;
  int (*boundFunc)(struct restraint * self);
   
  // integer to chapeau index, this will change
  // in replica exchange simulation
  int chid;
  int chdm;

  FILE * ofp;
  int outputFreq;
  int boutput;

} restraint;


void cfacvBanner ( void );

tamdOptStruct * New_tamdOptStruct ( double g, double kt, double dt, int riftyp);

smdOptStruct * New_smdOptStruct ( double target, int t0, int t1);

restraint * New_restraint ( double k, double z, int nCV, double * cvc, char * rftypstr, double zmin, double zmax, char * boundstr, double boundk, char * outfile, int outputFreq);
    

#ifndef MAXNBOR
#define MAXNBOR 250
#endif

typedef struct DataSpace {
  int N,iN; // number of centers
  int K,iK; // number of restraints
  int ch_num,ch_now; // number of chapeaus


  atomCenterStruct ** ac; // defined in centers.h
  restraint ** restr;
  chapeau ** ch;

  int ncv,icv;
  cv ** cv;
  
  char filename[255];
  int restrsavefreq;

  double ** R; // array of center cartesian coordinates R[i][0/1/2]

  // below are bits for doing TAMD/OTFP

  int doAnalyticalCalc; // OTFP on or off
  
  // for evolution of the parameterization
  double lamfric;
  double lamdt;

  // TODO: This should be controles from the input script
  int evolveAnalyticalParameters;
  int beginaccum;
  int beginsolve;

  int useTAMDforces; // indicates whether tamd forces are used to
		     // update auxiliary variables; other choice is
		     // the use forces available from the
		     // instantaneous parameterization of the
		     // analytical free energy (this would make the
		     // calculation conform to the heterogeneous
		     // multiscale method
  
} DataSpace;



//The random seed
unsigned short * Xi;

// Other subroutines
FILE * my_fopen ( char * name, char * code ) ;
DataSpace * NewDataSpace ( int N, int ncv, int K, long int seed );
int DataSpace_SetupChapeau ( DataSpace * ds, int numrep, int dm, double * min,
    int * nKnots, double * max, int * periodic, int beginaccum, int beginsolve, 
    int useTAMDforces, char * outfile, int outfreq, int outlevel, int nupdate);
chapeau * DataSpace_get_chapeauadress ( DataSpace * ds, int i );
int DataSpace_getN ( DataSpace * ds );
double * DataSpace_centerPos ( DataSpace * ds, int i );
int DataSpace_AddAtomCenter ( DataSpace * ds, int n, int * ind, double * m );
cv * DataSpace_add_cv ( DataSpace * ds, char * typ, int nind, int * ind,
		      double zmin, double zmax,char * boundf, double boundk, char * outfile, int outputFreq );

restraint * DataSpace_AddRestr  ( DataSpace * ds, double k, double targ, int nCV, double * cvc, char * rftypstr, double zmin, double zmax,char * boundf, double boundk,char * outfile, int outputFreq);
int restr_UpdateTamdOpt ( restraint * r, double g, double kt, double dt );
int restr_AddTamdOpt ( restraint * r, double g, double kt, double dt, int chid  , int chdm );
int restr_AddSmdOpt  ( restraint * r, double target, int t0, int t1 );
void restr_output  ( restraint * r );
int DataSpace_ComputeCVs ( DataSpace * ds );
int DataSpace_RestrainingForces ( DataSpace * ds, int first, int timestep );
double DataSpace_RestraintEnergy ( DataSpace * ds );
void DataSpace_ReportAll ( DataSpace * ds );
void DataSpace_ReportCV ( DataSpace * ds, int * active, double * res );
int DataSpace_checkdata ( DataSpace * ds );
int DataSpace_dump ( DataSpace * ds ); 
FILE * my_binfopen ( char * name, char * code, unsigned int outputLevel, DataSpace * ds );
void DataSpace_BinaryReportRestraints ( DataSpace * ds, int step, int outputlevel, FILE * fp );


// Evolve Functions
int cbd ( restraint * r, double f );
int uniformvelocity ( restraint * r, double f );


// Boundaries Functions
int SoftUpperWall ( restraint * r );
int SoftLowerWall ( restraint * r );
int SoftWalls ( restraint * r );
int pbc ( restraint * r );
int Periodic ( restraint * r );
int nada ( restraint * r );

// Potential Function
int HarmonicCart ( restraint * r );
int HarmonicCart_cutoff ( restraint * r );
int HarmonicCart_pbc ( restraint * r );
int HarmonicCart_cutoff_pbc ( restraint * r );

double restr_getz ( restraint * r );
double restr_getu ( restraint * r );
int restr_set_rchid ( restraint * r, DataSpace * ds, int chid);

void ds_saverestrains ( DataSpace * ds, char * filename );
void ds_loadrestrains ( DataSpace * ds, char * filename );


