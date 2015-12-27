
set terminal png;
j=0

set title "Replica id"
set output "job0_".j.".png"; j=j+1
plot for [i=0:5] 'output/'.i.'/job0.history' w l t 'p '.i
 


set title "Processor id"
set output "job0_".j.".png"; j=j+1
set terminal png; set output "Replica id"
plot for [i=0:5] 'output/'.i.'/job0_sort.history' w l t 'p '.i
 


set title "Sorted temperature traces"
set output "job0_".j.".png"; j=j+1
plot for [i=0:5] 'output/'.i.'/job0_sort.history' u 1:3 w l ls (i+1) lw 3 t 'p '.i\
   , for [i=0:5] 'output/'.i.'/job0_sort.history' u 1:4 w l ls (i+1)  notit
 


set title "Unsorted temperature traces"
set output "job0_".j.".png"; j=j+1
plot for [i=0:5] 'output/'.i.'/job0.history' u 1:3 w l ls (i+1) lw 3 t 'p '.i\
   , for [i=0:5] 'output/'.i.'/job0.history' u 1:4 w l ls (i+1)  notit
 
  

set title "Sorted energy traces"
set output "job0_".j.".png"; j=j+1
plot for [i=0:5] 'output/'.i.'/job0_sort.history' u 1:5 w l t 'p '.i
 


set title "Unsorted energy traces"
set output "job0_".j.".png"; j=j+1
plot for [i=0:5] 'output/'.i.'/job0.history' u 1:5 w l t 'p '.i
 
                      

set title "Unsorted CV traces"
set output "job0_".j.".png"; j=j+1
plot for [i=0:5] 'output/'.i.'/job0.s' u (($0+1)):1 w l t 'p '.i

 

set title "Sorted CV traces"
set output "job0_".j.".png"; j=j+1
plot for [i=0:5] 'output/'.i.'/job0_sort.s' u (($0+1)):1 w l t 'p '.i

  
