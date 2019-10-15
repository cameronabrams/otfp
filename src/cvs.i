%module cvs
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
#include "cvs.h"
%}
double * cv_access_ref ( cv * c, int i );
double * cv_access_ref2( cv * c, int i );
int set_line (cv * c);
