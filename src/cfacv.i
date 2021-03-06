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
extern DataSpace * NewDataSpace ( int N, int ncv, int K, long int seed ) ;
extern chapeau * DataSpace_get_chapeauadress ( DataSpace * ds, int i );
extern int DataSpace_AddAtomCenter ( DataSpace * ds, int n, int * ind, double * m );
cv * DataSpace_add_cv ( DataSpace * ds, char * typ, int nind, int * ind, double zmin, double zmax,char * boundf, double boundk, char * outfile, int outputFreq );
restraint * DataSpace_AddRestr  ( DataSpace * ds, double k, double targ, int nCV, double * cvc, char * rftypstr, double zmin, double zmax,char * boundf, double boundk,char * outfile, int outputFreq);
extern int DataSpace_SetupChapeau ( DataSpace * ds, int numrep, int dm, double * min, int * nKnots, double * max, int periodic, int beginaccum, int beginsolve, int useTAMDforces, char * outfile, int outfreq, int outlevel, int nupdate);
extern int restr_UpdateTamdOpt ( restraint * r, double g, double kt, double dt );
extern int restr_AddTamdOpt ( restraint * r, double g, double kt, double dt, int chid , int chdm );
extern int restr_AddSmdOpt  ( restraint * r, double target, int t0, int t1 );
extern int DataSpace_getN ( DataSpace * ds );
extern void DataSpace_ReportAll ( DataSpace * ds );
void DataSpace_ReportCV ( DataSpace * ds, int * active, double * res );
extern double * DataSpace_centerPos ( DataSpace * ds, int i );
extern int DataSpace_ComputeCVs ( DataSpace * ds );
extern int DataSpace_RestrainingForces ( DataSpace * ds, int first, int timestep );
extern double DataSpace_RestraintEnergy ( DataSpace * ds );
extern int DataSpace_checkdata ( DataSpace * ds );
extern int DataSpace_dump ( DataSpace * ds ); 
extern FILE * my_binfopen ( char * name, char * code, unsigned int outputLevel, DataSpace * ds );
extern void DataSpace_BinaryReportRestraints ( DataSpace * ds, int step, int outputlevel, FILE * fp );

extern int set_zsd_circle ( double x,double y, double xy, double s );
extern int set_zsd_ring ( double x,double y, double r1, double r2, double s );

extern double chapeau_evalf_1simplex ( chapeau * ch, double z );
extern char * chapeau_serialize ( chapeau * ch );
extern void chapeau_addserialized ( chapeau * ch, char * str );
extern void chapeau_setserialized ( chapeau *ch, char * str );
extern void chapeau_setserialized ( chapeau *ch, char * str );

extern void chapeau_set_peaks ( chapeau * ch, char * filename );
extern void chapeau_loadstate ( chapeau * ch, char * filename );
extern void chapeau_loadlambda ( chapeau * ch, char * filename );
extern void chapeau_savestate ( chapeau * ch, char * filename );
extern void ds_loadrestrains ( DataSpace * ds, char * filename );
extern void ds_saverestrains ( DataSpace * ds, char * filename );

extern double restr_getz ( restraint * r );
extern double restr_getu ( restraint * r );
extern int restr_set_rchid ( restraint * r, DataSpace * ds, int chid);

//extern void chapeau_setmref ( chapeau * ch, double z );
