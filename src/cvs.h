#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "measurements.h"
    
// headers
extern cv * cv_init ( char * typ, int nC, int * ind, 
    double zmin, double zmax, char * boundstr, double boundk, char * outfile, int outputFreq);

extern int cv_getityp ( char * typ );
extern char * cv_getstyp ( int ityp );
 
extern int calccv_line  ( cv * c, double ** R );
extern int calccv_rmsd  ( cv * c, double ** R );
extern int calccv_charmmd ( cv * c, double ** R );
extern int calccv_d2chap  ( cv * c, double ** R );
extern int set_d2chap    ( cv * c, char * filename );
extern int calccv_gauss2d ( cv * c, double ** R );
extern int calccv_cogx  ( cv * c, double ** R );
extern int calccv_cogy  ( cv * c, double ** R );
extern int calccv_cogz  ( cv * c, double ** R );
extern int calccv_x     ( cv * c, double ** R );
extern int calccv_y     ( cv * c, double ** R );
extern int calccv_z     ( cv * c, double ** R );
extern int calccv_s     ( cv * c, double ** R );
extern int calccv_halfbond  ( cv * c, double ** R );
extern int calccv_bond  ( cv * c, double ** R );
extern int calccv_bonds ( cv * c, double ** R );
extern int calccv_dihed ( cv * c, double ** R );
extern int calccv_angle ( cv * c, double ** R );
extern int calccv_zsd_circle ( cv * c, double ** R );
extern int calccv_zsd_xrange ( cv * c, double ** R );
extern int calccv_zsd_ring ( cv * c, double ** R );
extern int set_zsd_circle ( double x,double y, double xy, double s  ); 
extern int set_zsd_ring ( double x,double y, double r1, double r2, double s  ); 
           
extern double cdf(double x);

// Boundaries Functions
extern int cv_amd ( cv * c );
extern int cv_SoftUpperWall ( cv * c );
extern int cv_SoftLowerWall ( cv * c );
extern int cv_SoftWalls ( cv * c );
extern int cv_nada ( cv * c );
         
extern double * cv_access_ref ( cv * c, int i );
extern double * cv_access_ref2 ( cv * c, int i );

extern void cv_output ( cv * c );

extern int set_line (cv * c);

