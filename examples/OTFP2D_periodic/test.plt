
set terminal png;
j=0

set title "CV trace"
set output "job0_".j.".png"; j=j+1
plot [-180:180][-180:180] 'job0.s' w l notit

 

set title "Free energy job0"
set pm3d map
set output "job0_".j.".png"; j=j+1
splot 'job0_ch0.fes' notit
