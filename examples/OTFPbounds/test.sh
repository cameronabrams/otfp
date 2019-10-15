#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

# Run simulation
namd2 job0.namd > job0.log

for file in $(ls *bsp); do

  f=${file%.*}

  # Decoding the FES
  ../../src/catbinsp -f $file -ol 1 \
    | awk 'BEGIN{a=""} (NR>1&&$2!=a){a=$2;$1="";$2="";print $0}'\
    > ${f}.LAMBDA

  # Adding x-range
  awk -v N=201 -v min=2.0 -v max=4.8 \
      'END{
        for (i=1;i<=N;i++){
          print min+(max-min)/(N-1)*i,$i
        }
      }' ${f}.LAMBDA > ${f}.fes
 
  # Plotting
  gnuplot << HEREGNUPLOT

  set terminal png;

  set title "Free energy"
  set output "${f}.png"
  plot [2.6:4.6] '${f}.fes' w l notit

HEREGNUPLOT
    
done


for file in $(ls *dcd); do

  f=${file%.*}

  #Computing traces
  vmd -dispdev text << HERETCL
   
  mol new ../systems/buta.psf type psf first 0 last -1 step 1 filebonds 1 autobonds 0 waitfor all
  mol addfile $file type dcd first 0 last -1 step 1 filebonds 1 autobonds 0 waitfor all

  proc distance { i1 i2 } {
    set sel1 [atomselect 0 "index \$i1"]
    set sel2 [atomselect 0 "index \$i2"]
    set coord1 [lindex [\$sel1 get {x y z}] 0]
    set coord2 [lindex [\$sel2 get {x y z}] 0]
    return [vecdist \$coord1 \$coord2]
  } 

  #open file for writing
  set fil [open "${f}.s" w]
  set all [atomselect top all]
  set nframes [molinfo top get numframes]
  puts \$nframes
  for {set frame 0} {\$frame < \$nframes} {incr frame} {
      animate goto \$frame
      puts \$fil [distance 0 10] 
  }
  close \$fil

  exit            

HERETCL
 
  # Plotting the trace
  gnuplot << HEREGNUPLOT

  set terminal png;

  set title "CV trace"
  set output "${f}.png"
  unset key
  plot '${f}.s' u (\$0+1):1 w l, 3.8, 3.5

HEREGNUPLOT

done 

#Final message
echo "Test completed sucesfully, compare the plots obtained with the reference"

