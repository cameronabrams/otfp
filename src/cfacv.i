%module cfa_cvlibc
%include carrays.i
%array_functions(double, array);
%array_functions(int, arrayint);
%inline %{
double get_double(double *a, int index) {
	return a[index];
}
%}
%{
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "cfacv.h"
%}
extern void cfacvBanner ( void );
extern FILE * my_fopen ( char * name, char * code );
extern DataSpace * NewDataSpace ( int N, int M, int K, long int seed ) ;
extern int DataSpace_SetupPairCalc ( DataSpace * ds, double Ox, double Oy, double Oz, double Lx, double Ly, double Lz, double cutoff, double nlcutoff, int beginEvolve, int useTAMDforces, int reportParamFreq, double spline_min, int nKnots, char * splineoutputfile, int splineoutputfreq, int splineoutputlevel, int lamupdateinterval );
extern int * DataSpace_pairmasks ( DataSpace * ds );
extern int DataSpace_AddAtomCenter ( DataSpace * ds, int n, int * ind, double * m );
extern int DataSpace_AddCV ( DataSpace * ds, char * typ, int nind, int * ind ) ;
extern int DataSpace_AddRestr ( DataSpace * ds, double k, double targ, int nCV, double * cvc, char * rftypstr, double zmin, double zmax  );
extern int DataSpace_AddTamdOpt ( DataSpace * ds, int ir, double g, double kt, double dt );
extern int DataSpace_AddSmdOpt  ( DataSpace * ds, int ir, double target, int t0, int t1 );
extern int DataSpace_getN ( DataSpace * ds );
extern void DataSpace_ReportAll ( DataSpace * ds );
void DataSpace_ReportCV ( DataSpace * ds, int * active, double * res );
extern double * DataSpace_centerPos ( DataSpace * ds, int i );
extern int DataSpace_ComputeCVs ( DataSpace * ds );
extern int DataSpace_RestrainingForces ( DataSpace * ds, int first, int timestep );
extern double DataSpace_RestraintEnergy ( DataSpace * ds );
extern void DataSpace_ReportRestraints ( DataSpace * ds, int step, int outputlevel, FILE * fp );
extern int DataSpace_SetRestraints ( DataSpace * ds, double * rval );
extern int DataSpace_checkdata ( DataSpace * ds );
extern int DataSpace_dump ( DataSpace * ds ); 
extern FILE * my_binfopen ( char * name, char * code, unsigned int outputLevel, DataSpace * ds );
extern void DataSpace_BinaryReportRestraints ( DataSpace * ds, int step, int outputlevel, FILE * fp );
extern int DataSpace_InitKnots ( DataSpace * ds, char * filename );
