#pragma once
#ifndef _CFACV_H_
#define _CFACV_H_
 
// #ifndef MAXNBOR
// #define MAXNBOR 250
// #endif

#include "cvs_obj.h"
#include "chapeau_obj.h"
#include "centers.h"
         
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


#endif
 
