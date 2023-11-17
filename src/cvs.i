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
#include "cvs_obj.h"
#include "chapeau.h"
#include "cvs.h"
%}
%include "chapeau.h"
%include "cvs.h"
