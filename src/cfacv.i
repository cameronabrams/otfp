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
#include <execinfo.h>
#include <signal.h>  
#include <unistd.h>  
#include "cfacv_obj.h"
#include "cfacv.h"
%}

%include "cfacv.h"
