#!/usr/bin/gnuplot

IF1='NaCl.namd';IF2='NaCl.LAMBDA';eval INTERM('NaCl_both.tex','(8.5*2)cm,7.5cm')

# IF1 the namd input file
# IF2 the LAMBDA file
eval checkfile(IF1)
eval checkfile(IF2)


dr(rmax,rmin,knots)=real(abs(rmax-rmin))/(knots-1)  
  
rmin=real(system('grep "^set SPLINEMIN" '.IF1." | awk '{print $3}'"))
rmax=real(system('grep "^set CUTOFF" '.IF1." | awk '{print $3}'"))
knots=real(system("awk '{print NF;exit}' ".IF2))

print 'grid: ',rmin,rmax,knots

DR=dr(rmax,rmin,knots)

n=system('wc -l '.IF2.' | cut -d" " -f1 ')
print 'lines: ',n

F2=IF2.'_25'
F3=IF2.'_50'
F4=IF2.'_75'
F5=IF2.'_100'
system 'sed -n "'.int(n*1./4.).'p;'.int(n*1./4.).'q" '.IF2." | fmt -w 11 > ".F2
system 'sed -n "'.int(n*2./4.).'p;'.int(n*2./4.).'q" '.IF2." | fmt -w 11 > ".F3
system 'sed -n "'.int(n*3./4.).'p;'.int(n*3./4.).'q" '.IF2." | fmt -w 11 > ".F4
system 'tail -n 1 '.IF2." | fmt -w 11 > ".F5



xc=0.2  # margen derecho   
yc=0.05  # margen arriba    
xof=0.2  # margen izquierdo
yof=0.2  # margen abajo    
eval mploton(2,1)
               
set multiplot

eval mplot
aux=rmin+knots/2*DR
set arrow from aux,1.4 to aux,1.2 lc rgb 'black'
set xtics rmin,1,rmax-0.5
set xlabel 'r (\AA)'
set ylabel 'g (kcal/mol)' offset 2,0
plot [rmin:rmax][-2.5:2.5] 0 lc rgb 'black' t '',\
         F2 u (rmin+$0*DR):1 w l ls 1 lw 3 lc @ZUB5_1 t '\%'.int(100./4.)  , \
         F3 u (rmin+$0*DR):1 w l ls 1 lw 3 lc @ZUB5_2 t '\%'.int(200./4.)  , \
         F4 u (rmin+$0*DR):1 w l ls 1 lw 3 lc @ZUB5_3 t '\%'.int(300./4.)  , \
         F5 u (rmin+$0*DR):1 w l ls 1 lw 3 lc @ZUB5_4 t '\%'.int(100)

unset arrow                           

eval mplot
set y2tics
set xtics auto
set xlabel 'evolution' offset 2,0
plot [][-2.5:2.5] IF2 u (column(knots/2)) w l ls 1 lw 3 lc @ZUB5_1 t ''
 
unset multiplot
