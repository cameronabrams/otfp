#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

# Run simulations
if [ -n "$(which mpirun 2>/dev/null)" ]; then
  mpirun -n 12 namd2 job0.namd 2>&1 | tee job0.log
  mpirun -n 12 namd2 job1.namd 2>&1 | tee job1.log
else
  namd2 +p12 job0.namd 2>&1 | tee job0.log
  namd2 +p12 job1.namd 2>&1 | tee job1.log
fi

for file in $(ls *bsp); do

  f=${file%.*}

  # Decoding the FES
  ../../src/catbinsp -f $file -ol 1 \
    | awk 'BEGIN{a=""} (NR>1&&$2!=a){a=$2;$1="";$2="";print $0}'\
    > ${f}.LAMBDA

  # Unfolding the FES
  awk -v N=151 \
      -v min=-3.14159265358979 -v max=3.14159265358979 \
      'END{
        aux=(max-min)/(N-1)
        for (j=0;j<N;j++){
          for (i=0;i<N;i++){
            a=i+j*N+1
            print min+aux*j,min+aux*i,$a
          }
          print ""
        }
      }' ${f}.LAMBDA > ${f}.fes

  # Plotting the FES
  gnuplot << HEREGNUPLOT

  set terminal png;
  set size square
  set encoding utf8

  Is0(x)=((x==0.) ? 1/0 : x)

  stats '${f}.fes' u (Is0(\$3))
  fesmin=STATS_min

  set cbrange [0:25]
  set output "${f}_fes.png"
  set ylabel 'φ' norotate
  set xlabel 'ψ' offset 1,0
  set pm3d map
  splot [-pi:pi][-pi:pi]\
      '${f}.fes' u 2:1:(Is0(\$3)-fesmin) notit

HEREGNUPLOT

done

#Final message
echo "Test completed sucesfully, compare the plots obtained with the reference"

