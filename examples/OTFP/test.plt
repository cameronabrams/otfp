
set terminal png;
j=0

set title "CV trace"
set output "job0_".j.".png"; j=j+1
plot 'job0.s' u (($0+1)):1 w l notit

 

set title "Free energy job0"
set output "job0_".j.".png"; j=j+1
plot 'job0_ch0.fes' w l notit


set title "CV trace"
set output "job1_".j.".png"; j=j+1
plot 'job1.s' u (($0+1)):1 w l notit

 

set title "Free energy job1 initial"
set output "job1_".j.".png"; j=j+1
plot 'job1_ch0.fes1' w l notit
 
set title "Free energy job1"
set output "job1_".j.".png"; j=j+1
plot 'job1_ch0.fes' w l notit
 
