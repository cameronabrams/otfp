#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

# Run simulations
namd2 rmsd.namd 2>&1 | tee rmsd.log
namd2 line.namd 2>&1 | tee line.log

# for file in $(ls *bsp); do
for file in line_ch0.bsp; do

  f=${file%.*}

  # Decoding the FES
  ../../src/catbinsp -f $file -ol 1 \
    | awk 'BEGIN{a=""} (NR>1&&$2!=a){a=$2;$1="";$2="";print $0}'\
    > ${f}.LAMBDA

  # Adding x-range
  awk -v N=201 -v min=-.5 -v max=1.5 \
      'END{
        for (i=1;i<=N;i++){
          print min+(max-min)/(N-1)*i,$i
        }
      }' ${f}.LAMBDA > ${f}.fes
 
  # Plotting
  gnuplot << HEREGNUPLOT

  set terminal png;

  Is0(x)=((x==0.) ? 1/0 : x)

  stats '${f}.fes' u (Is0(\$2))
  fesmin=STATS_min

  set title "Free energy"
  set output "${f}.png"
  plot [-.5:1.5][:5] '${f}.fes' u 1:(Is0(\$2)-fesmin) w l notit

HEREGNUPLOT
done

#Final message
echo "Test completed sucesfully, compare the plots obtained with the reference"

