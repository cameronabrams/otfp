#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct ATOMCENTERSTRUCT {
  int n;
  int * ind;
  double * m;
  double M;
} atomCenterStruct;


extern atomCenterStruct * New_atomCenterStruct ( int n );

typedef struct CENTERSTRUCT * pcenterStruct;
typedef struct CENTERSTRUCT { 
  int id;
  int maxN;
  int iN;
  int * mList;
  double rg;
  double cm[3];
  pcenterStruct left; /* only used in residue blocking */
  pcenterStruct right;
} centerStruct;

extern centerStruct * New_centerStruct ( int id, int maxN );
extern void centerStuct_addMember ( centerStruct * c, int i );

extern void center_rg ( centerStruct * c, double * x, double * y, double * z );

extern int rgyr_sort ( centerStruct * c, int * bin, double * x, double * y, double * z, 
		int nAtom, int minAtom, double * rg, unsigned int Seed  );
extern centerStruct * Null_centerStruct ( void );
extern int bin_sort ( int * bin, double * x, double * y, double * z, int nAtoms, int nCenters, int nCycles, 
	       unsigned int Seed );



